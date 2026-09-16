
import os, json
import sys
import numpy, math

from scipy import interpolate
from scipy.optimize import nnls
from schainpy.model import data
from schainpy.model.proc.jroproc_base import ProcessingUnit, Operation, MPDecorator
from schainpy.model.data.jrodata import Voltage,hildebrand_sekhon
from schainpy.utils import log
from time import time, mktime, strptime, gmtime, ctime
from scipy.optimize import least_squares
import datetime
import collections.abc
import csv
import ast #new added
from scipy import signal

try:
    from schainpy.model.proc import fitacf_guess
    from schainpy.model.proc import fitacf_fit_short
    from schainpy.model.proc import fitacf_acf2
    from schainpy.model.proc import full_profile_profile
except:
    log.warning('Missing Faraday fortran libs')

import matplotlib.pyplot as plt

class VoltageProc(ProcessingUnit):

    def __init__(self):

        ProcessingUnit.__init__(self)

        self.dataOut = Voltage()
        self.flip = 1
        self.setupReq = False

    def run(self, runNextUnit = 0):

        if self.dataIn.type == 'AMISR':
            self.__updateObjFromAmisrInput()

        if self.dataIn.type == 'Voltage':
            self.dataOut.copy(self.dataIn)
            self.dataOut.runNextUnit = runNextUnit

            #print(self.dataOut.flagNoData)
            #print(self.dataOut.flagNoMoreData)
            # print("data shape: ", self.dataOut.data.shape)
            #exit(1)

        #for i in range(len(self.dataOut.heightList)):
        #    print(self.dataOut.heightList[i])
        #print("len(self.dataOut.heightList)", len(self.dataOut.heightList))
        #print("self.dataOut.heightList[0]", self.dataOut.heightList[0])
        #print("self.dataOut.heightList[-1]", self.dataOut.heightList[-1])
        
        

    def __updateObjFromAmisrInput(self):

        self.dataOut.timeZone = self.dataIn.timeZone
        self.dataOut.dstFlag = self.dataIn.dstFlag
        self.dataOut.errorCount = self.dataIn.errorCount
        self.dataOut.useLocalTime = self.dataIn.useLocalTime

        self.dataOut.flagNoData = self.dataIn.flagNoData
        self.dataOut.data = self.dataIn.data
        self.dataOut.utctime = self.dataIn.utctime
        self.dataOut.channelList = self.dataIn.channelList
        #self.dataOut.timeInterval = self.dataIn.timeInterval
        self.dataOut.heightList = self.dataIn.heightList
        self.dataOut.nProfiles = self.dataIn.nProfiles

        self.dataOut.nCohInt = self.dataIn.nCohInt
        self.dataOut.ippSeconds = self.dataIn.ippSeconds
        self.dataOut.frequency = self.dataIn.frequency

        self.dataOut.azimuth = self.dataIn.azimuth
        self.dataOut.zenith = self.dataIn.zenith

        self.dataOut.beam.codeList = self.dataIn.beam.codeList
        self.dataOut.beam.azimuthList = self.dataIn.beam.azimuthList
        self.dataOut.beam.zenithList = self.dataIn.beam.zenithList


    # print(self.dataOut.data)

class selectChannels(Operation):

    def run(self, dataOut, channelList):

        channelIndexList = []
        self.dataOut = dataOut
        for channel in channelList:
            if channel not in self.dataOut.channelList:
                raise ValueError("Channel %d is not in %s" %(channel, str(self.dataOut.channelList)))

            index = self.dataOut.channelList.index(channel)
            channelIndexList.append(index)
        self.selectChannelsByIndex(channelIndexList)

        # print("data shape : ", self.dataOut.data.shape)


        return self.dataOut


    def selectChannelsByIndex(self, channelIndexList):
        """
        Selecciona un bloque de datos en base a canales segun el channelIndexList

        Input:
            channelIndexList    :    lista sencilla de canales a seleccionar por ej. [2,3,7]

        Affected:
            self.dataOut.data
            self.dataOut.channelIndexList
            self.dataOut.nChannels
            self.dataOut.m_ProcessingHeader.totalSpectra
            self.dataOut.systemHeaderObj.numChannels
            self.dataOut.m_ProcessingHeader.blockSize

        Return:
            None
        """

        for channelIndex in channelIndexList:
            if channelIndex not in self.dataOut.channelIndexList:
                raise ValueError("The value %d in channelIndexList is not valid" %channelIndex)

        if self.dataOut.type == 'Voltage':
            if self.dataOut.flagDataAsBlock:
                """
                Si la data es obtenida por bloques, dimension = [nChannels, nProfiles, nHeis]
                """
                data = self.dataOut.data[channelIndexList,:,:]
            else:
                data = self.dataOut.data[channelIndexList,:]

            self.dataOut.data = data
            # self.dataOut.channelList = [self.dataOut.channelList[i] for i in channelIndexList]
            self.dataOut.channelList = range(len(channelIndexList))

        elif self.dataOut.type == 'Spectra':
            data_spc = self.dataOut.data_spc[channelIndexList, :]
            data_dc = self.dataOut.data_dc[channelIndexList, :]

            self.dataOut.data_spc = data_spc
            self.dataOut.data_dc = data_dc

            # self.dataOut.channelList = [self.dataOut.channelList[i] for i in channelIndexList]
            self.dataOut.channelList = range(len(channelIndexList))
            self.__selectPairsByChannel(channelIndexList)

        return 1

    def __selectPairsByChannel(self, channelList=None):

        if channelList == None:
            return

        pairsIndexListSelected = []
        for pairIndex in self.dataOut.pairsIndexList:
            # First pair
            if self.dataOut.pairsList[pairIndex][0] not in channelList:
                continue
            # Second pair
            if self.dataOut.pairsList[pairIndex][1] not in channelList:
                continue

            pairsIndexListSelected.append(pairIndex)

        if not pairsIndexListSelected:
            self.dataOut.data_cspc = None
            self.dataOut.pairsList = []
            return

        self.dataOut.data_cspc = self.dataOut.data_cspc[pairsIndexListSelected]
        self.dataOut.pairsList = [self.dataOut.pairsList[i]
                                  for i in pairsIndexListSelected]

        return

class CombineChannels(Operation):
    '''Digital hybrid implementation'''

    def run(self, dataout, sum_list=[], sub_list=[]):
        '''
        Input:
            sum_list    : list of pairs [[0,2],[1,3]]
            sub_list    : list of pairs [[2,4],[6,8]]
        '''
        #print(dataout.data[0, :, :])
        tmp = []

        if sub_list:
            for i, j in sub_list:
                if dataout.flagDataAsBlock:
                    tmp.append( dataout.data[i, :, :] - dataout.data[j, :, :])
                else:
                    tmp.append(dataout.data[i,:] - dataout.data[j,:,:])
        
        if sum_list:
            for i, j in sum_list:
                if dataout.flagDataAsBlock:
                    tmp.append(dataout.data[i, :, :] + dataout.data[j, :, :])
                    #tmp.append(numpy.sum(dataout.data[i, :, :],dataout.data[j, :, :]))
                else:
                    tmp.append(dataout.data[i,:] + dataout.data[j,:,:])
                    

        

        dataout.data = numpy.array(tmp)
        dataout.channelList = range(len(tmp))

        return dataout


class CombineChannels_V2(Operation):
    '''Digital hybrid implementation'''

    def run(self, dataout, ch_list=[], comb_list=[]):
        '''
        Input:
            ch_list    : list of pairs [[0,2],[1,3]] or single-index [[4], [5]]
            comb_list  : list of operations ['sum', 'sub', 'none']
        '''
        comb_list = [s.strip().lower() for s in comb_list]
        tmp = []
        is_block = dataout.flagDataAsBlock

        for (channels, comb) in zip(ch_list, comb_list):
            i = channels[0]

            if comb == 'none':
                data = dataout.data[i, :, :] if is_block else dataout.data[i, :]
            
            elif comb in ('sum', 'sub'):
                j = channels[1]
                a = dataout.data[i, :, :] if is_block else dataout.data[i, :]
                b = dataout.data[j, :, :] if is_block else dataout.data[j, :]
                data = a + b if comb == 'sum' else a - b

            else:
                raise ValueError(f"Unknown combination operation: '{comb}'")

            tmp.append(data)

        dataout.data = numpy.array(tmp)
        dataout.channelList = list(range(len(tmp)))
        return dataout

class saturatedBlock(Operation):
    '''
    Created by C. Portilla
    hardcoded for the moment to Hybrid experiment
    '''

    def ca_cfar_alpha(self, N, pfa):
        """
        Analytical CA-CFAR scale factor for exponential noise:
        alpha = N * (pfa^{-1/N} - 1)
        Works as a decent initial approximation for SO-CFAR.
        """
        return N * (pfa**(-1.0 / N) - 1.0)


    def so_cfar_1d(self, power, T=12, G=2, pfa=1e-6, alpha=None, pad_mode='edge'):
        """
        Smallest-Of CFAR (SO-CFAR) for 1-D power profile.
        Parameters
        ----------
        power : 1D np.array (float)
            Power values (non-negative). NaN are allowed (treated as missing).
        T : int
            Number of training cells PER SIDE (left and right).
        G : int
            Number of guard cells per side.
        pfa : float
            Desired probability of false alarm (used to compute alpha if not provided).
        alpha : float or None
            If provided, use this scaling factor for threshold. Otherwise compute
            CA-CFAR alpha and use it as approximation.
        pad_mode : {'edge','constant'}
            Padding behavior for edges; 'edge' repeats edge values, 'constant' pads with zero.
        Returns
        -------
        detections : boolean array same shape as power
            True where detections occurred (satellite candidates).
        threshold : float array same shape as power
            Computed threshold for each CUT (NaN where not computable).
        ref_estimate : float array
            reference noise estimate (min of left/right averages).
        notes:
            - cells near edges where full training window does not fit are returned as False
            - input may contain NaN; training windows containing NaN are treated by ignoring NaNs
        Complexity: O(N) using cumulative sums.
        """
        power = numpy.asarray(power, dtype=float)
        N = len(power)
        if alpha is None:
            alpha = self.ca_cfar_alpha(2 * T, pfa)  # use CA alpha as initial approximation

        # Minimum number of samples needed around a CUT
        pad = T + G
        # pad array to simplify indexing; choose edge padding to avoid introducing zeros in quiet backgrounds
        if pad_mode == 'edge':
            p = numpy.pad(power, pad_width=pad, mode='edge')
        else:
            p = numpy.pad(power, pad_width=pad, mode='constant', constant_values=0.0)

        # cumulative sum but treating NaN: we'll compute sums and counts ignoring NaN
        vals = numpy.nan_to_num(p, nan=0.0)
        counts = (~numpy.isnan(p)).astype(float)

        csum = numpy.concatenate(([0.0], numpy.cumsum(vals)))  # length len(p)+1
        ccount = numpy.concatenate(([0.0], numpy.cumsum(counts)))

        # helper to get sum over interval [a,b) in padded indices
        def window_sum(a, b):
            return csum[b] - csum[a]
        def window_count(a, b):
            return ccount[b] - ccount[a]

        detections = numpy.zeros(N, dtype=bool)
        threshold = numpy.full(N, numpy.nan, dtype=float)
        ref_est = numpy.full(N, numpy.nan, dtype=float)

        # prepare padded index of original element i -> i_p = i + pad
        idx = numpy.arange(N) + pad

        for k, i_p in enumerate(idx):
            # left window [i_p - G - T, i_p - G)
            l0 = int(i_p - G - T)
            l1 = int(i_p - G)
            # right window [i_p + G + 1, i_p + G + 1 + T)
            r0 = int(i_p + G + 1)
            r1 = int(r0 + T)

            # safety check (should be fine due to padding)
            if l0 < 0 or r1 > len(p):
                continue

            sum_l = window_sum(l0, l1)
            cnt_l = window_count(l0, l1)
            sum_r = window_sum(r0, r1)
            cnt_r = window_count(r0, r1)

            # compute averages ignoring NaNs (if too few samples, skip)
            if cnt_l < max(1, 0.5 * T) or cnt_r < max(1, 0.5 * T):
                # not enough valid training cells -> leave as no detection (or you could choose to adapt)
                continue

            avg_l = sum_l / cnt_l
            avg_r = sum_r / cnt_r

            ref = min(avg_l, avg_r)  #0.5*(avg_l+avg_r) #min(avg_l, avg_r) #
            ref_est[k] = ref
            thr = alpha * ref
            threshold[k] = thr

            vcut = power[k]
            if numpy.isnan(vcut):
                # cannot test a NaN CUT; treat as not detection (or optionally mark)
                continue
            # Decision
            if vcut > thr:
                detections[k] = True


        return detections, threshold, ref_est

    def run(self, dataOut, limit = 70, mode = 0):
        
        #self.setup()
        #print("TEST is equal Ch 0 ", dataOut.data[0,0:15,:])
        #print("TEST is equal Ch 1 ", dataOut.data[0,0:15,:])
        z = numpy.abs(dataOut.data)
        
        # Working place
        # data z  (4, 150, 334) -> (4,128,334)
        z_new = z [:,:128,:]
        dataOut_extra = dataOut.data[:, 128:, :]
        z_new = z_new.reshape(4, 8, 16, 334)
        dataOut_z_new = z_new.reshape(4, 8, 16, 334)

        idx = []
        jars_fix = 1 # to solve, the samplpes are not enough correct

        # Identify the saturated profile groups
        for i in range(8):
            hardTarget = False

            if mode == 0:
                sample_trh = slice(40,200 - jars_fix)
                time_text = datetime.datetime.utcfromtimestamp(dataOut.utctime)
                if time_text.hour > 23 or time_text.hour < 5: sample_trh = slice(45,200 - jars_fix) #10

                hardTarget = (z_new[1,i,:,sample_trh] > limit).any()

            elif mode == 1:
                # Array to save (and plot) satellites
                dataOut.sat_indices = numpy.zeros([2,200 - jars_fix])

                Ch_list = [0,1]
                detections = []
                aux_hardTarget = []
                sample_trh = slice(25,200 - jars_fix)
                for Ch in Ch_list:
                    self.z_new = numpy.average(z_new[Ch,i,:,sample_trh], axis=0)
                    # Typical values T=12, G=8, pfa=4e-2 Hybrid experiment
                    detection,_,_ = self.so_cfar_1d(self.z_new, T=12, G=8, pfa=4e-2, pad_mode='edge')
                    detections.append(detection)
                    aux_hardTarget.append(sum(detection) >= 2)

                    sat_indices = numpy.where(detection,0.125, 0) # Convert the detection in a fraction of block
                    #if detection.any() == True: print("aaBs", sat_indices)
                    dataOut.sat_indices[Ch, sample_trh] += sat_indices

                detections = numpy.logical_or.reduce(detections)
                #hardTarget = detections.any() # At least 1 target in both Ch
                #hardTarget = sum(detections) >= 2 # At least 2 targets in both Ch
                hardTarget = numpy.array(aux_hardTarget).any()  # At least 2 targets in any Ch

            else:
                return log.warning("Any filter BLock mode selected")

            if hardTarget: # Save index and print
                idx.append(i)
                print(f"Debris detected at {i}", (z_new[1,i,:,sample_trh] > limit).any())
        

        # Deals what to do with saturated profile groups
        if len(idx) < 7:#len(idx) != 8:
            candidates = numpy.setdiff1d(numpy.arange(8), idx)
            result = numpy.array([candidates[numpy.abs(candidates - i).argmin()] for i in idx])
            for n,i in enumerate(idx):
                #print(i,result[n]) # use to print each profile is being changed indexes
                #z_new[:,i,:,23:200] = z_new[:,result[n],:,23:200]
                get_block = lambda x: slice(x*16, (x+1)*16)
                dataOut.data[:,get_block(i),:] = dataOut.data[:,get_block(result[n]),:]
        else:
            print("All profiles saturated")

        # Temporal correction of JARS2 sampling issue - AllISR 2025
        '''print(dataOut.data[0,10,197], dataOut.data[0,10,198], dataOut.data[0,10,199], dataOut.data[0,10,200], dataOut.data[0,10,201])
        print("Jars correction")'''
        dataOut.data[:,:,199] = dataOut.data[:,:,199- jars_fix]

        return dataOut


class saturatedBlock_Sat(saturatedBlock):
    '''
    Created by C. Portilla
    hardcoded for the moment to Hybrid experiment
    '''

    def run(self, dataOut, limit = 70, mode = 0):

        #self.setup()
        #print("TEST is equal Ch 0 ", dataOut.data[0,0:15,:])
        #print("TEST is equal Ch 1 ", dataOut.data[0,0:15,:])
        z = numpy.abs(dataOut.data)
        
        # Working place
        # data z  (4, 150, 334) -> (4,128,334)
        z_new = z [:,:128,:]
        dataOut_extra = dataOut.data[:, 128:, :]
        z_new = z_new.reshape(4, 8, 16, 334)
        dataOut_z_new = z_new.reshape(4, 8, 16, 334)

        idx = []
        jars_fix = 1 # to solve, the samplpes are not enough correct
        sample_trh = slice(40,200 - jars_fix)
        time_text = datetime.datetime.utcfromtimestamp(dataOut.utctime)
        if time_text.hour > 23 or time_text.hour < 5: sample_trh = slice(45,200 - jars_fix) #10

        # Identify the saturated profile groups
        dataOut.sat_indices = numpy.zeros([2,200 - jars_fix])
        for i in range(8):
            hardTarget = False

            if mode == 0:
                hardTarget = (z_new[1,i,:,sample_trh] > limit).any()

            elif mode == 1:
                Ch_list = [0,1]
                detections = []
                aux_hardTarget = []
                sample_trh = slice(25,200 - jars_fix)
                for Ch in Ch_list:
                    self.z_new = numpy.average(z_new[Ch,i,:,sample_trh], axis=0)
                    # Typical values T=12, G=8, pfa=4e-2 Hybrid experiment
                    detection,_,_ = self.so_cfar_1d(self.z_new, T=12, G=8, pfa=4e-2, pad_mode='edge')
                    #print("aa", detection)
                    detections.append(detection)
                    aux_hardTarget.append(sum(detection) >= 2)

                    sat_indices = numpy.where(detection,0.125, 0)
                    #if detection.any() == True: print("aaBs", sat_indices)
                    dataOut.sat_indices[Ch, sample_trh] += sat_indices

                #detections = numpy.array(sorted(list(set(numpy.concatenate(detections))))) # join satellites indexes
                detections = numpy.logical_or.reduce(detections)
                
                hardTarget = numpy.array(aux_hardTarget).any()  # At least 2 targets in any Ch

            else:
                return log.warning("Any filter BLock mode selected")

            if hardTarget:
                idx.append(i)
                print(f"Debris detected at {i}", (z_new[1,i,:,sample_trh] > limit).any()) # use to print each profile is being changed
                #print(z_new[1,i,:,sample_trh])
        

        # Deals what to do with saturated profile groups
        if len(idx) < 7:#len(idx) != 8:
            candidates = numpy.setdiff1d(numpy.arange(8), idx)
            result = numpy.array([candidates[numpy.abs(candidates - i).argmin()] for i in idx])
            for n,i in enumerate(idx):
                #print(i,result[n]) # use to print each profile is being changed indexes
                #z_new[:,i,:,23:200] = z_new[:,result[n],:,23:200]
                get_block = lambda x: slice(x*16, (x+1)*16)
                dataOut.data[:,get_block(i),:] = dataOut.data[:,get_block(result[n]),:]
        else:
            print("All profiles saturated")

        # Temporal correction of JARS2 sampling issue - AllISR 2025
        '''print(dataOut.data[0,10,197], dataOut.data[0,10,198], dataOut.data[0,10,199], dataOut.data[0,10,200], dataOut.data[0,10,201])
        print("Jars correction")'''
        dataOut.data[:,:,199] = dataOut.data[:,:,199- jars_fix]

        return dataOut


class LagsReshape150(Operation):
    '''
    Written by R. Flores
    '''
    """Operation to reshape input data into (Channels,Profiles(with same lag),Heights,Lags) and heights reconstruction.

    Parameters:
    -----------


    Example
    --------

    op = proc_unit.addOperation(name='LagsReshape')


    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

        self.buffer=None
        self.buffer_HR=None
        self.buffer_HRonelag=None

    def LagDistribution(self,dataOut):

        dataOut.datapure=numpy.copy(dataOut.data[:,0:dataOut.NSCAN,:])
        self.buffer = numpy.zeros((dataOut.nChannels,
                                   int(dataOut.NSCAN/dataOut.DPL),
                                   dataOut.nHeights,dataOut.DPL),
                                  dtype='complex')

        for j in range(int(self.buffer.shape[1]/2)):
            for i in range(dataOut.DPL):
                if j+1==int(self.buffer.shape[1]/2) and i+1==dataOut.DPL:
                    self.buffer[:,2*j:,:,i]=dataOut.datapure[:,2*i+int(2*j*dataOut.DPL):,:]
                else:
                    self.buffer[:,2*j:2*(j+1),:,i]=dataOut.datapure[:,2*i+int(2*j*dataOut.DPL):2*(i+1)+int(2*j*dataOut.DPL),:]

        return self.buffer

    def HeightReconstruction(self,dataOut):

        self.buffer_HR = numpy.zeros((int(dataOut.NSCAN/dataOut.DPL),
                                   dataOut.nHeights,dataOut.DPL),
                                  dtype='complex')

        for i in range(int(dataOut.DPL)): #Only channel B
            if i==0:
                self.buffer_HR[:,:,i]=dataOut.datalags[1,:,:,i]
            else:
                self.buffer_HR[:,:,i]=self.HRonelag(dataOut,i)

        return self.buffer_HR


    def HRonelag(self,dataOut,whichlag):
        self.buffer_HRonelag = numpy.zeros((int(dataOut.NSCAN/dataOut.DPL),
                                   dataOut.nHeights),
                                  dtype='complex')

        for i in range(self.buffer_HRonelag.shape[0]):
            for j in range(dataOut.nHeights):
                if j+int(2*whichlag)<dataOut.nHeights:
                    self.buffer_HRonelag[i,j]=dataOut.datalags[1,i,j+2*whichlag,whichlag]
                else:
                    if whichlag!=10:
                        self.buffer_HRonelag[i,j]=dataOut.datalags[1,i,(j+2*whichlag)%dataOut.nHeights,whichlag+1]
                    else:
                        if i+2<self.buffer_HRonelag.shape[0]:
                            self.buffer_HRonelag[i,j]=dataOut.datalags[1,i+2,(j+2*whichlag)%dataOut.nHeights,0]
                        else: #i+1==self.buffer_HRonelag.shape[0]:
                            self.buffer_HRonelag[i,j]=dataOut.datalags[1,i,(j+2*whichlag)%dataOut.nHeights,whichlag]

        return self.buffer_HRonelag



    def run(self,dataOut,DPL=11,NSCAN=132):

        dataOut.DPL=DPL
        dataOut.NSCAN=NSCAN
        dataOut.paramInterval=0#int(dataOut.nint*dataOut.header[7][0]*2 )
        dataOut.lat=-11.95
        dataOut.lon=-76.87
        dataOut.datalags=None

        dataOut.datalags=numpy.copy(self.LagDistribution(dataOut))
        dataOut.datalags[1,:,:,:]=self.HeightReconstruction(dataOut)
        #print(dataOut.datalags[0,:,])
        dataOut.data = numpy.reshape(dataOut.datalags,(2,132,dataOut.nHeights))

        return dataOut

class LagsReshapeHP(Operation):
    '''
    Written by R. Flores
    '''
    """Operation to reshape input data into (Channels,Profiles(with same lag),Heights,Lags) and heights reconstruction.

    Parameters:
    -----------


    Example
    --------

    op = proc_unit.addOperation(name='LagsReshape')


    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

        self.buffer=None
        self.buffer_HR=None
        self.buffer_HRonelag=None

    def LagDistribution(self,dataOut):

        dataOut.datapure=numpy.copy(dataOut.data[:,0:dataOut.NSCAN,:])
        self.buffer = numpy.zeros((dataOut.nChannels,
                                   int(dataOut.NSCAN/dataOut.DPL),
                                   dataOut.nHeights,dataOut.DPL),
                                  dtype='complex')

        for j in range(int(self.buffer.shape[1]/2)):
            for i in range(dataOut.DPL):
                if j+1==int(self.buffer.shape[1]/2) and i+1==dataOut.DPL:
                    self.buffer[:,2*j:,:,i]=dataOut.datapure[:,2*i+int(2*j*dataOut.DPL):,:]
                else:
                    self.buffer[:,2*j:2*(j+1),:,i]=dataOut.datapure[:,2*i+int(2*j*dataOut.DPL):2*(i+1)+int(2*j*dataOut.DPL),:]

        return self.buffer

    def HeightReconstruction(self,dataOut):

        self.buffer_HR = numpy.zeros((int(dataOut.NSCAN/dataOut.DPL),
                                   dataOut.nHeights,dataOut.DPL),
                                  dtype='complex')

        for i in range(int(dataOut.DPL)): #Only channel B
            if i==0:
                self.buffer_HR[:,:,i]=dataOut.datalags[1,:,:,i]
            else:
                self.buffer_HR[:,:,i]=self.HRonelag(dataOut,i)

        return self.buffer_HR


    def HRonelag(self,dataOut,whichlag):
        self.buffer_HRonelag = numpy.zeros((int(dataOut.NSCAN/dataOut.DPL),
                                   dataOut.nHeights),
                                  dtype='complex')

        for i in range(self.buffer_HRonelag.shape[0]):
            for j in range(dataOut.nHeights):
                if j+int(2*whichlag)<dataOut.nHeights:
                    self.buffer_HRonelag[i,j]=dataOut.datalags[1,i,j+2*whichlag,whichlag]
                else:
                    if whichlag!=10:
                        self.buffer_HRonelag[i,j]=dataOut.datalags[1,i,(j+2*whichlag)%dataOut.nHeights,whichlag+1]
                    else:
                        if i+2<self.buffer_HRonelag.shape[0]:
                            self.buffer_HRonelag[i,j]=dataOut.datalags[1,i+2,(j+2*whichlag)%dataOut.nHeights,0]
                        else: #i+1==self.buffer_HRonelag.shape[0]:
                            self.buffer_HRonelag[i,j]=dataOut.datalags[1,i,(j+2*whichlag)%dataOut.nHeights,whichlag]

        return self.buffer_HRonelag



    def run(self,dataOut,DPL=11,NSCAN=132):

        dataOut.DPL=DPL
        dataOut.NSCAN=NSCAN
        dataOut.paramInterval=0#int(dataOut.nint*dataOut.header[7][0]*2 )
        dataOut.lat=-11.95
        dataOut.lon=-76.87
        dataOut.datalags=None

        dataOut.datalags=numpy.copy(self.LagDistribution(dataOut))
        dataOut.datalags[1,:,:,:]=self.HeightReconstruction(dataOut)

        return dataOut

class LagsReshapeDP(Operation):
    '''
    Written by R. Flores
    '''
    """Operation to reshape input data into (Channels,Profiles(with same lag),Heights,Lags) and heights reconstruction.

    Parameters:
    -----------


    Example
    --------

    op = proc_unit.addOperation(name='LagsReshape')


    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

        self.buffer=None

    def LagDistribution(self,dataOut):

        self.buffer = numpy.zeros((dataOut.nChannels,
                                   int(2*2*dataOut.NSCAN/dataOut.NLAG),
                                   dataOut.NDP,dataOut.DPL),
                                  dtype='complex')

        indProfile = numpy.arange(0,dataOut.NSCAN,1)//8

        #dataOut.nNoiseProfiles = dataOut.nProfiles-dataOut.NSCAN

        for i in range(2):
            if i==0:
                aux = 0
            else:
                aux =16
            for j in range(dataOut.NDP):
                for k in range(int(dataOut.NSCAN)):

                    n=dataOut.lagind[k%dataOut.NLAG]

                    data_ChA=dataOut.data[0,k,dataOut.NRANGE+j+i*dataOut.NDT]#-dataOut.dc[0]

                    if dataOut.NRANGE+j+i*dataOut.NDT+2*n<dataOut.read_samples:

                        data_ChB=dataOut.data[1,k,dataOut.NRANGE+j+i*dataOut.NDT+2*n]#-dataOut.dc[1]
                        #print(data_ChB)
                        #exit(1)
                        #print("*1*")

                    else:
                        #print(i,j,n)
                        #exit(1)

                        if k+1<int(dataOut.NSCAN):
                            data_ChB=dataOut.data[1,k+1,(dataOut.NRANGE+j+i*dataOut.NDT+2*n)%dataOut.NDP]
                            #print(data_ChB)
                            #print("*2*")
                            #exit(1)
                        if k+1==int(dataOut.NSCAN):
                            data_ChB=dataOut.data[1,k,(dataOut.NRANGE+j+i*dataOut.NDT+2*n)%dataOut.NDP]
                            #print("*3*")
                    #if n == 7 and j == 65:
                        #print(k)
                        #print(data_ChB)
                                #exit(1)
                    if n == 8 or n == 9 or n == 10:
                        self.buffer[0,int((aux+indProfile[k]-1)/2),j,n] = data_ChA
                        self.buffer[1,int((aux+indProfile[k]-1)/2),j,n] = data_ChB
                    elif n == 1 or n == 2 or n == 7:
                        self.buffer[0,int((aux+indProfile[k])/2),j,n] = data_ChA
                        self.buffer[1,int((aux+indProfile[k])/2),j,n] = data_ChB
                    else:
                        self.buffer[0,aux+indProfile[k],j,n] = data_ChA
                        self.buffer[1,aux+indProfile[k],j,n] = data_ChB

        #FindMe
        pa1 = 20
        pa2 = 10

        #print(self.buffer[0,:,pa1,pa2])
        #print(self.buffer[1,:,pa1,pa2])
        '''
        print(sum(self.buffer[0,:,pa1,pa2]))
        print(sum(self.buffer[1,:,pa1,pa2]))
        #exit(1)
        '''

        '''
        for pa1 in range(67):
            print(sum(self.buffer[0,:,pa1,pa2]))
            print(sum(self.buffer[1,:,pa1,pa2]))
            '''

        '''
        import matplotlib.pyplot as plt
        fft = numpy.fft.fft(self.buffer[0,:,pa1,pa2])
        fft2 = fft*numpy.conjugate(fft)
        fft2 = fft2.real
        fft2 = numpy.fft.fftshift(fft2)
        '''
        #print("before",fft2)
        #plt.plot(fft2)
        #plt.show()
        #import time
        #time.sleep(5)
        #plt.close('all')
        #exit(1)
        return self.buffer



    def run(self,dataOut,DPL=11,NSCAN=128,lagind=(0,1,2,3,4,5,6,7,0,3,4,5,6,8,9,10),lagfirst=(1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1), NLAG = 16, NRANGE = 200):

        dataOut.DPL=DPL
        dataOut.NSCAN=NSCAN
        dataOut.NLAG = NLAG
        deltaHeight   =  dataOut.heightList[1] - dataOut.heightList[0]
        dataOut.NRANGE = NRANGE
        dataOut.read_samples=int(dataOut.nHeights)
        #print(dataOut.read_samples)
        #print(dataOut.nHeights)
        #exit(1)
        dataOut.NDP = dataOut.NDT = int((dataOut.nHeights-dataOut.NRANGE)/2)
        dataOut.heightList = numpy.arange(dataOut.NDP) *deltaHeight# + dataOut.heightList[0]
        #dataOut.NDP = dataOut.NDT = int(dataOut.nHeights/2)#int((dataOut.nHeights-dataOut.NRANGE)/2)
        #print(dataOut.NDP)
        #print(dataOut.heightList)
        dataOut.paramInterval=0#int(dataOut.nint*dataOut.header[7][0]*2 )
        dataOut.lat=-11.95
        dataOut.lon=-76.87
        dataOut.datalags=None
        dataOut.lagind=lagind
        dataOut.lagfirst=lagfirst


        #print(dataOut.data[1,:12,:15])
        #exit(1)
        #print(numpy.shape(dataOut.data))
        dataOut.datalags = numpy.copy(self.LagDistribution(dataOut))
        #print(numpy.shape(dataOut.datalags))
        #exit(1)
        #print("AFTER RESHAPE DP")

        dataOut.data = dataOut.data[:,:,200:]
        #print(numpy.shape(dataOut.data))
        #exit(1)


        return dataOut

class LagsReshapeDP_V2(Operation):
    '''
    Written by R. Flores
    '''
    """Operation to reshape input data into (Channels,Profiles(with same lag),Heights,Lags) and heights reconstruction.

    Parameters:
    -----------


    Example
    --------

    op = proc_unit.addOperation(name='LagsReshape')


    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

        self.buffer=None
        self.data_buffer = []

    def setup(self,dataOut,DPL,NSCAN,NLAG,NRANGE,lagind,lagfirst):
        dataOut.DPL=DPL
        dataOut.NSCAN=NSCAN
        dataOut.NLAG = NLAG
        deltaHeight   =  dataOut.heightList[1] - dataOut.heightList[0]
        dataOut.NRANGE = NRANGE
        dataOut.read_samples=int(dataOut.nHeights)
        #print(dataOut.read_samples)
        #print(dataOut.nHeights)
        #exit(1)
        dataOut.NDP = dataOut.NDT = int((dataOut.nHeights-dataOut.NRANGE)/2)
        dataOut.heightList = numpy.arange(dataOut.NDP) *deltaHeight# + dataOut.heightList[0]
        #dataOut.NDP = dataOut.NDT = int(dataOut.nHeights/2)#int((dataOut.nHeights-dataOut.NRANGE)/2)
        #print(dataOut.NDP)
        #print(dataOut.heightList)
        dataOut.paramInterval=0#int(dataOut.nint*dataOut.header[7][0]*2 )
        dataOut.lat=-11.95
        dataOut.lon=-76.87
        dataOut.datalags=None
        dataOut.lagind=lagind
        dataOut.lagfirst=lagfirst


    def LagDistribution(self,dataOut):

        self.buffer = numpy.zeros((dataOut.nChannels,
                                   int(2*2*dataOut.NSCAN/dataOut.NLAG),
                                   dataOut.NDP,dataOut.DPL),
                                  dtype='complex')

        indProfile = numpy.arange(0,dataOut.NSCAN,1)//8

        #dataOut.nNoiseProfiles = dataOut.nProfiles-dataOut.NSCAN

        for i in range(2):
            if i==0:
                aux = 0
            else:
                aux =16
            for j in range(dataOut.NDP):
                for k in range(int(dataOut.NSCAN)):

                    n=dataOut.lagind[k%dataOut.NLAG]

                    data_ChA=dataOut.data[0,k,dataOut.NRANGE+j+i*dataOut.NDT]#-dataOut.dc[0]

                    if dataOut.NRANGE+j+i*dataOut.NDT+2*n<dataOut.read_samples:

                        data_ChB=dataOut.data[1,k,dataOut.NRANGE+j+i*dataOut.NDT+2*n]#-dataOut.dc[1]
                        #print(data_ChB)
                        #exit(1)
                        #print("*1*")

                    else:
                        #print(i,j,n)
                        #exit(1)

                        if k+1<int(dataOut.NSCAN):
                            data_ChB=dataOut.data[1,k+1,(dataOut.NRANGE+j+i*dataOut.NDT+2*n)%dataOut.NDP]
                            #print(data_ChB)
                            #print("*2*")
                            #exit(1)
                        if k+1==int(dataOut.NSCAN):
                            data_ChB=dataOut.data[1,k,(dataOut.NRANGE+j+i*dataOut.NDT+2*n)%dataOut.NDP]
                            #print("*3*")
                    #if n == 7 and j == 65:
                        #print(k)
                        #print(data_ChB)
                                #exit(1)
                    if n == 8 or n == 9 or n == 10:
                        self.buffer[0,int((aux+indProfile[k]-1)/2),j,n] = data_ChA
                        self.buffer[1,int((aux+indProfile[k]-1)/2),j,n] = data_ChB
                    elif n == 1 or n == 2 or n == 7:
                        self.buffer[0,int((aux+indProfile[k])/2),j,n] = data_ChA
                        self.buffer[1,int((aux+indProfile[k])/2),j,n] = data_ChB
                    else:
                        self.buffer[0,aux+indProfile[k],j,n] = data_ChA
                        self.buffer[1,aux+indProfile[k],j,n] = data_ChB

        #FindMe
        pa1 = 20
        pa2 = 10

        #print(self.buffer[0,:,pa1,pa2])
        #print(self.buffer[1,:,pa1,pa2])
        '''
        print(sum(self.buffer[0,:,pa1,pa2]))
        print(sum(self.buffer[1,:,pa1,pa2]))
        #exit(1)
        '''

        '''
        for pa1 in range(67):
            print(sum(self.buffer[0,:,pa1,pa2]))
            print(sum(self.buffer[1,:,pa1,pa2]))
            '''

        '''
        import matplotlib.pyplot as plt
        fft = numpy.fft.fft(self.buffer[0,:,pa1,pa2])
        fft2 = fft*numpy.conjugate(fft)
        fft2 = fft2.real
        fft2 = numpy.fft.fftshift(fft2)
        '''
        #print("before",fft2)
        #plt.plot(fft2)
        #plt.show()
        #import time
        #time.sleep(5)
        #plt.close('all')
        #exit(1)
        return self.buffer



    def run(self,dataOut,DPL=11,NSCAN=128,lagind=(0,1,2,3,4,5,6,7,0,3,4,5,6,8,9,10),lagfirst=(1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1), NLAG = 16, NRANGE = 200):

        if not self.isConfig:
            self.setup(dataOut,DPL,NSCAN,NLAG,NRANGE,lagind,lagfirst)
            self.isConfig = True

        #print(dataOut.data[1,:12,:15])
        #exit(1)
        #print(numpy.shape(dataOut.data))
        #print(dataOut.profileIndex)

        if not dataOut.flagDataAsBlock:

            dataOut.flagNoData = True
            #print("nProfiles: ",dataOut.nProfiles)
            #if dataOut.profileIndex == 140:
            #print("id: ",dataOut.profileIndex)
            if dataOut.profileIndex == dataOut.nProfiles-1:
                #print("here")
                #print(dataOut.data.shape)
                self.data_buffer.append(dataOut.data)
                dataOut.data = numpy.transpose(numpy.array(self.data_buffer),(1,0,2))
                #print(dataOut.data.shape)
                #print(numpy.sum(dataOut.data))
                #print(dataOut.data[1,100,:])
                #exit(1)
                dataOut.datalags = numpy.copy(self.LagDistribution(dataOut))
                #print(numpy.shape(dataOut.datalags))
                #exit(1)
                #print("AFTER RESHAPE DP")

                dataOut.data = dataOut.data[:,:,200:]
                self.data_buffer = []
                dataOut.flagDataAsBlock = True
                dataOut.flagNoData = False

                deltaHeight   =  dataOut.heightList[1] - dataOut.heightList[0]
                dataOut.heightList = numpy.arange(dataOut.NDP) *deltaHeight# + dataOut.heightList[0]
                #exit(1)
                #print(numpy.sum(dataOut.datalags))
                #exit(1)

            else:
                self.data_buffer.append(dataOut.data)
            #print(numpy.shape(dataOut.data))
            #exit(1)
        else:
            #print(dataOut.data.shape)
            #print(numpy.sum(dataOut.data))
            #print(dataOut.data[1,100,:])
            #exit(1)
            dataOut.datalags = numpy.copy(self.LagDistribution(dataOut))
            #print(dataOut.datalags.shape)
            dataOut.data = dataOut.data[:,:,200:]
            deltaHeight   =  dataOut.heightList[1] - dataOut.heightList[0]
            dataOut.heightList = numpy.arange(dataOut.NDP) * deltaHeight# + dataOut.heightList[0]
            #print(dataOut.nHeights)
            #print(numpy.sum(dataOut.datalags))
            #exit(1)

        return dataOut

class LagsReshapeLP(Operation):
    '''
    Written by R. Flores
    '''
    """Operation to reshape input data into (Channels,Profiles(with same lag),Heights,Lags) and heights reconstruction.

    Parameters:
    -----------


    Example
    --------

    op = proc_unit.addOperation(name='LagsReshape')


    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

        self.buffer=None


    def LagDistributionLP(self,dataOut):

        buffer=dataOut.data
        self.buffer = numpy.zeros((dataOut.nChannels,
                                   dataOut.NSCAN,
                                   dataOut.NRANGE,dataOut.NLAG),
                                  dtype='complex64')


            #self.dataOut.nptsfft2=150
        self.cnorm=float((dataOut.nProfiles-dataOut.NSCAN)/dataOut.NSCAN)
        self.lagp0=numpy.zeros((dataOut.NLAG,dataOut.NSCAN,dataOut.NRANGE),'complex64')
        self.lagp1=numpy.zeros((dataOut.NLAG,dataOut.NSCAN,dataOut.NRANGE),'complex64')
        self.lagp2=numpy.zeros((dataOut.NLAG,dataOut.NSCAN,dataOut.NRANGE),'complex64')
        #self.lagp3=numpy.zeros((dataOut.NLAG,dataOut.nProfiles-dataOut.NSCAN,dataOut.NRANGE),'complex64')
        self.lagp3=numpy.zeros((dataOut.NLAG,dataOut.NSCAN,dataOut.NRANGE),'complex64')
        #self.lagp4=numpy.zeros((dataOut.NLAG,dataOut.NRANGE,dataOut.NAVG),'complex64')

        for i in range(dataOut.nChannels):
            #buffer_dc=dataOut.dc[i]
            for j in range(dataOut.NRANGE):

                range_for_n=numpy.min((dataOut.NRANGE-j,dataOut.NLAG))

                buffer_aux=numpy.conj(buffer[i,:dataOut.nProfiles,j])#-buffer_dc)
                for n in range(range_for_n):

                    c=(buffer_aux)*(buffer[i,:dataOut.nProfiles,j+n])#-buffer_dc)

                    if i==0:
                        self.lagp0[n,:,j]=c[:dataOut.NSCAN]
                        self.lagp3[n,:dataOut.nProfiles-dataOut.NSCAN,j]=c[dataOut.NSCAN:]/self.cnorm
                    elif i==1:
                        self.lagp1[n,:,j]=c[:dataOut.NSCAN]
                    elif i==2:
                        self.lagp2[n,:,j]=c[:dataOut.NSCAN]

        '''
        self.lagp0=numpy.conj(self.lagp0)
        self.lagp1=numpy.conj(self.lagp1)
        self.lagp2=numpy.conj(self.lagp2)
        self.lagp3=numpy.conj(self.lagp3)
        '''
        self.buffer[0,:,:,:] = numpy.transpose(numpy.conj(self.lagp0),(1,2,0))
        self.buffer[1,:,:,:] = numpy.transpose(numpy.conj(self.lagp1),(1,2,0))
        self.buffer[2,:,:,:] = numpy.transpose(numpy.conj(self.lagp2),(1,2,0))
        self.buffer[3,:,:,:] = numpy.transpose(numpy.conj(self.lagp3),(1,2,0))

        #print(self.buffer[3,:,100,15])
        print("Sum: ", sum(self.buffer[3,:,199,2]))
        #print(self.cnorm)
        exit(1)


        return self.buffer

    def LagDistributionLP_V2(self,dataOut):

        buffer=dataOut.data
        self.buffer = numpy.zeros((dataOut.NLAG,
                                   dataOut.NSCAN,
                                   dataOut.NRANGE,dataOut.nChannels),
                                  dtype='complex128')


            #self.dataOut.nptsfft2=150
        self.cnorm=float((dataOut.nProfiles-dataOut.NSCAN)/dataOut.NSCAN)
        #dataOut.nNoiseProfiles = dataOut.nProfiles-dataOut.NSCAN
        self.lagp0=numpy.zeros((dataOut.NLAG,dataOut.NSCAN,dataOut.NRANGE),'complex128')
        self.lagp1=numpy.zeros((dataOut.NLAG,dataOut.NSCAN,dataOut.NRANGE),'complex128')
        self.lagp2=numpy.zeros((dataOut.NLAG,dataOut.NSCAN,dataOut.NRANGE),'complex128')
        #self.lagp3=numpy.zeros((dataOut.NLAG,dataOut.nProfiles-dataOut.NSCAN,dataOut.NRANGE),'complex64')
        self.lagp3=numpy.zeros((dataOut.NLAG,dataOut.NSCAN,dataOut.NRANGE),'complex128')
        #self.lagp4=numpy.zeros((dataOut.NLAG,dataOut.NRANGE,dataOut.NAVG),'complex64')

        for i in range(dataOut.nChannels):
            #buffer_dc=dataOut.dc[i]
            for j in range(dataOut.NRANGE):

                range_for_n=numpy.min((dataOut.NRANGE-j,dataOut.NLAG))

                #buffer_aux=buffer[i,:dataOut.nProfiles,j]#-buffer_dc)
                for n in range(range_for_n):

                    c=buffer[i,:dataOut.nProfiles,j+n]#-buffer_dc)

                    if i==0:
                        self.lagp0[n,:,j]=c[:dataOut.NSCAN]
                        self.lagp3[n,:dataOut.nProfiles-dataOut.NSCAN,j]=c[dataOut.NSCAN:]#/self.cnorm
                    elif i==1:
                        self.lagp1[n,:,j]=c[:dataOut.NSCAN]
                    elif i==2:
                        self.lagp2[n,:,j]=c[:dataOut.NSCAN]

        '''
        self.lagp0=numpy.conj(self.lagp0)
        self.lagp1=numpy.conj(self.lagp1)
        self.lagp2=numpy.conj(self.lagp2)
        self.lagp3=numpy.conj(self.lagp3)
        '''
        self.buffer[:,:,:,0] = self.lagp0
        self.buffer[:,:,:,1] = self.lagp1
        self.buffer[:,:,:,2] = self.lagp2
        self.buffer[:,:,:,3] = self.lagp3

        #print(self.buffer[3,:,100,15])
        #print(sum(self.buffer[3,:,199,2]))
        #print(self.cnorm)
        #exit(1)


        return self.buffer


    def run(self,dataOut,DPL=11,NSCAN=128,lagind=(0,1,2,3,4,5,6,7,0,3,4,5,6,8,9,10),lagfirst=(1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1),NLAG=16,NRANGE=200):

        dataOut.DPL=DPL
        dataOut.NSCAN=NSCAN
        dataOut.NLAG = NLAG
        dataOut.NRANGE = dataOut.nHeights
        #print(dataOut.NRANGE)
        #exit(1)
        dataOut.paramInterval=0#int(dataOut.nint*dataOut.header[7][0]*2 )
        dataOut.lat=-11.95
        dataOut.lon=-76.87
        dataOut.datalags=None
        dataOut.lagind=lagind
        dataOut.lagfirst=lagfirst

        #self.LagDistributionHP(dataOut)
        dataOut.datalags = numpy.copy(self.LagDistributionLP_V2(dataOut))


        #dataOut.final_noise = dataOut.getNoise(Profmin_index=128, Profmax_index=150)[:2]
        #print(dataOut.final_noise.shape)

        dataOut.channelList = range(16)
        #print("AFTER RESHAPE LP")
        #print(dataOut.datalags[15])
        #print("DONE")
        #exit(1)

        return dataOut

class LagsReshapeHP2(Operation):
    '''
    Written by R. Flores
    '''
    """Operation to reshape input data into (Channels,Profiles(with same lag),Heights,Lags) and heights reconstruction.

    Parameters:
    -----------


    Example
    --------

    op = proc_unit.addOperation(name='LagsReshape')


    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

        self.buffer=None

    def LagDistribution(self,dataOut):

        nChannelsDP = 2

        self.buffer = numpy.zeros((nChannelsDP,
                                   int(2*2*dataOut.NSCAN/dataOut.NLAG),
                                   dataOut.NDP,dataOut.DPL),
                                  dtype='complex')

        indProfile = numpy.arange(0,dataOut.NSCAN,1)//8

        for i in range(2):
            if i==0:
                aux = 0
            else:
                aux =16
            for j in range(dataOut.NDP):
                for k in range(int(dataOut.NSCAN)):

                    n=dataOut.lagind[k%dataOut.NLAG]

                    data_ChA=dataOut.data[0,k,dataOut.NRANGE+j+i*dataOut.NDT]#-dataOut.dc[0]

                    if dataOut.NRANGE+j+i*dataOut.NDT+2*n<dataOut.read_samples:

                        data_ChB=dataOut.data[1,k,dataOut.NRANGE+j+i*dataOut.NDT+2*n]#-dataOut.dc[1]

                    else:

                        if k+1<int(dataOut.NSCAN):
                            data_ChB=dataOut.data[1,k+1,(dataOut.NRANGE+j+i*dataOut.NDT+2*n)%dataOut.NDP]
                            #print(data_ChB)
                            #exit(1)
                        if k+1==int(dataOut.NSCAN):
                            data_ChB=dataOut.data[1,k,(dataOut.NRANGE+j+i*dataOut.NDT+2*n)%dataOut.NDP]

                    if n == 8 or n == 9 or n == 10:
                        self.buffer[0,int((aux+indProfile[k]-1)/2),j,n] = data_ChA
                        self.buffer[1,int((aux+indProfile[k]-1)/2),j,n] = data_ChB
                    elif n == 1 or n == 2 or n == 7:
                        self.buffer[0,int((aux+indProfile[k])/2),j,n] = data_ChA
                        self.buffer[1,int((aux+indProfile[k])/2),j,n] = data_ChB
                    else:
                        self.buffer[0,aux+indProfile[k],j,n] = data_ChA
                        self.buffer[1,aux+indProfile[k],j,n] = data_ChB

        #print(self.buffer[0,:,65,8])
        #print(sum(self.buffer[0,:,65,8]))
        #print(sum(self.buffer[1,:,65,8]))
        #exit(1)
        print(sum(self.buffer[0,:,65,7]))
        print(sum(self.buffer[1,:,65,7]))

        exit(1)
        return self.buffer

    def LagDistributionHP(self,dataOut):

        buffer=dataOut.data
        self.buffer = numpy.zeros((dataOut.nChannels,
                                   dataOut.NSCAN,
                                   dataOut.NRANGE,dataOut.NLAG),
                                  dtype='complex64')


            #self.dataOut.nptsfft2=150
        self.cnorm=float((dataOut.nProfiles-dataOut.NSCAN)/dataOut.NSCAN)
        self.lagp0=numpy.zeros((dataOut.NLAG,dataOut.NSCAN,dataOut.NRANGE),'complex64')
        self.lagp1=numpy.zeros((dataOut.NLAG,dataOut.NSCAN,dataOut.NRANGE),'complex64')
        self.lagp2=numpy.zeros((dataOut.NLAG,dataOut.NSCAN,dataOut.NRANGE),'complex64')
        #self.lagp3=numpy.zeros((dataOut.NLAG,dataOut.nProfiles-dataOut.NSCAN,dataOut.NRANGE),'complex64')
        self.lagp3=numpy.zeros((dataOut.NLAG,dataOut.NSCAN,dataOut.NRANGE),'complex64')
        #self.lagp4=numpy.zeros((dataOut.NLAG,dataOut.NRANGE,dataOut.NAVG),'complex64')

        for i in range(dataOut.nChannels):
            #buffer_dc=dataOut.dc[i]
            for j in range(dataOut.NRANGE):

                range_for_n=numpy.min((dataOut.NRANGE-j,dataOut.NLAG))

                buffer_aux=numpy.conj(buffer[i,:dataOut.nProfiles,j])#-buffer_dc)
                for n in range(range_for_n):

                    c=(buffer_aux)*(buffer[i,:dataOut.nProfiles,j+n])#-buffer_dc)

                    if i==0:
                        self.lagp0[n,:,j]=c[:dataOut.NSCAN]
                        self.lagp3[n,:dataOut.nProfiles-dataOut.NSCAN,j]=c[dataOut.NSCAN:]#/self.cnorm
                    elif i==1:
                        self.lagp1[n,:,j]=c[:dataOut.NSCAN]
                    elif i==2:
                        self.lagp2[n,:,j]=c[:dataOut.NSCAN]

        '''
        self.lagp0=numpy.conj(self.lagp0)
        self.lagp1=numpy.conj(self.lagp1)
        self.lagp2=numpy.conj(self.lagp2)
        self.lagp3=numpy.conj(self.lagp3)
        '''
        self.buffer[0,:,:,:] = numpy.transpose(numpy.conj(self.lagp0),(1,2,0))
        self.buffer[1,:,:,:] = numpy.transpose(numpy.conj(self.lagp1),(1,2,0))
        self.buffer[2,:,:,:] = numpy.transpose(numpy.conj(self.lagp2),(1,2,0))
        self.buffer[3,:,:,:] = numpy.transpose(numpy.conj(self.lagp3),(1,2,0))
        '''
        print(self.buffer[3,:,100,15])
        print(sum(self.buffer[3,:,100,15]))
        print(self.cnorm)
        exit(1)
        '''

        return self.buffer


    def run(self,dataOut,DPL=11,NSCAN=128,lagind=(0,1,2,3,4,5,6,7,0,3,4,5,6,8,9,10),lagfirst=(1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1),NLAG=16,NRANGE=200):

        dataOut.DPL=DPL
        dataOut.NSCAN=NSCAN
        dataOut.NLAG = NLAG
        dataOut.NRANGE = NRANGE
        dataOut.NDP = dataOut.NDT = int((dataOut.nHeights-dataOut.NRANGE)/2)
        dataOut.paramInterval=0#int(dataOut.nint*dataOut.header[7][0]*2 )
        dataOut.lat=-11.95
        dataOut.lon=-76.87
        dataOut.datalags=None
        dataOut.lagind=lagind
        dataOut.lagfirst=lagfirst
        dataOut.read_samples=len(dataOut.heightList)

        dataOut.datalagsDP=numpy.copy(self.LagDistribution(dataOut))
        #self.LagDistributionHP(dataOut)
        dataOut.datalagsLP = numpy.copy(self.LagDistributionHP(dataOut))

        B = numpy.zeros((2,96,67,11),dtype='complex')
        #C = numpy.zeros((2,128,133,11),dtype='complex')
        D = numpy.zeros((2,128,67,5),dtype='complex')

        dataOut.datalagsDP = numpy.append(dataOut.datalagsDP,B,axis=1)
        #dataOut.datalagsDP = numpy.append(dataOut.datalagsDP,C,axis=2)
        dataOut.datalagsDP = numpy.append(dataOut.datalagsDP,D,axis=3)

        ##First 2 "channels" DP, 4 other "channels" LP
        #dataOut.datalags = numpy.append(dataOut.datalagsDP,dataOut.datalagsLP,axis=0)
        dataOut.datalags = numpy.append(dataOut.datalagsDP,dataOut.datalagsLP[:2,:],axis=2)
        dataOut.datalags = numpy.append(dataOut.datalags,dataOut.datalagsLP[2:4,:],axis=2)
        #print(dataOut.datalags.shape)
        #exit(1)
        '''
        print(dataOut.datalags.shape)
        print(dataOut.datalags[0,:,65,7])
        print(sum(dataOut.datalags[0,:,65,7]))
        print(sum(dataOut.datalags[1,:,65,7]))
        exit(1)
        '''

        return dataOut

class selectHeights(Operation):

    def run(self, dataOut, minHei=None, maxHei=None, minIndex=None, maxIndex=None):
        """
        Selecciona un bloque de datos en base a un grupo de valores de alturas segun el rango
        minHei <= height <= maxHei

        Input:
            minHei    :    valor minimo de altura a considerar
            maxHei    :    valor maximo de altura a considerar

        Affected:
            Indirectamente son cambiados varios valores a travez del metodo selectHeightsByIndex

        Return:
            1 si el metodo se ejecuto con exito caso contrario devuelve 0
        """

        # print(self.dataOut.heightList)
        # print(self.dataOut.heightList[0])
        # print(self.dataOut.heightList[-1])

        self.dataOut = dataOut

        if minHei and maxHei:
        #if 1:
            if minHei == None:
               minHei = self.dataOut.heightList[0]

            if maxHei == None:
               maxHei = self.dataOut.heightList[-1]

            if (minHei < self.dataOut.heightList[0]):
                minHei = self.dataOut.heightList[0]

            if (maxHei > self.dataOut.heightList[-1]):
                maxHei = self.dataOut.heightList[-1]

            minIndex = 0
            maxIndex = 0
            heights = self.dataOut.heightList

            inda = numpy.where(heights >= minHei)
            indb = numpy.where(heights <= maxHei)

            try:
                minIndex = inda[0][0]
            except:
                minIndex = 0

            try:
                maxIndex = indb[0][-1]
            except:
                maxIndex = len(heights)

        self.selectHeightsByIndex(minIndex, maxIndex)
        #print(self.dataOut.nHeights)


        return self.dataOut

    def selectHeightsByIndex(self, minIndex, maxIndex):
        """
        Selecciona un bloque de datos en base a un grupo indices de alturas segun el rango
        minIndex <= index <= maxIndex

        Input:
            minIndex    :    valor de indice minimo de altura a considerar
            maxIndex    :    valor de indice maximo de altura a considerar

        Affected:
            self.dataOut.data
            self.dataOut.heightList

        Return:
            1 si el metodo se ejecuto con exito caso contrario devuelve 0
        """

        if self.dataOut.type == 'Voltage':
            if (minIndex < 0) or (minIndex > maxIndex):
                raise ValueError("Height index range (%d,%d) is not valid" % (minIndex, maxIndex))

            if (maxIndex >= self.dataOut.nHeights):
                maxIndex = self.dataOut.nHeights

            #voltage
            if self.dataOut.flagDataAsBlock:
                """
                Si la data es obtenida por bloques, dimension = [nChannels, nProfiles, nHeis]
                """
                data = self.dataOut.data[:,:, minIndex:maxIndex]
            else:
                data = self.dataOut.data[:, minIndex:maxIndex]

            #         firstHeight = self.dataOut.heightList[minIndex]

            self.dataOut.data = data
            self.dataOut.heightList = self.dataOut.heightList[minIndex:maxIndex]

            if self.dataOut.nHeights <= 1:
                raise ValueError("selectHeights: Too few heights. Current number of heights is %d" %(self.dataOut.nHeights))
        elif self.dataOut.type == 'Spectra':
            if (minIndex < 0) or (minIndex > maxIndex):
                raise ValueError("Error selecting heights: Index range (%d,%d) is not valid" % (
                    minIndex, maxIndex))

            if (maxIndex >= self.dataOut.nHeights):
                maxIndex = self.dataOut.nHeights - 1

            # Spectra
            data_spc = self.dataOut.data_spc[:, :, minIndex:maxIndex + 1]
            if hasattr(self.dataOut, "ByLags"):
                if self.dataOut.ByLags:
                    self.dataOut.dataLag_spc = self.dataOut.dataLag_spc[:, :, minIndex:maxIndex + 1]

            data_cspc = None
            if self.dataOut.data_cspc is not None:
                data_cspc = self.dataOut.data_cspc[:, :, minIndex:maxIndex + 1]
                if hasattr(self.dataOut, "ByLags"):
                    if self.dataOut.ByLags:
                        self.dataOut.dataLag_cspc = self.dataOut.dataLag_cspc[:, :, minIndex:maxIndex + 1]

            data_dc = None
            if self.dataOut.data_dc is not None:
                data_dc = self.dataOut.data_dc[:, minIndex:maxIndex + 1]
                if hasattr(self.dataOut, "ByLags"):
                    if self.dataOut.ByLags:
                        self.dataOut.dataLag_dc = self.dataOut.dataLag_dc[:, minIndex:maxIndex + 1]

            self.dataOut.data_spc = data_spc
            self.dataOut.data_cspc = data_cspc
            self.dataOut.data_dc = data_dc

            self.dataOut.heightList = self.dataOut.heightList[minIndex:maxIndex + 1]

        return 1


class filterByHeights(Operation):

    def run(self, dataOut, window):
        # len(heightList) = 990
        # heightList = [0, 1.5, 3.0, 4.5, ..., 1480.5, 1482.0, 1483.5] aprox. 1500 km
        deltaHeight = dataOut.heightList[1] - dataOut.heightList[0]
        
        if window == None:
            window = (dataOut.radarControllerHeaderObj.txA/dataOut.radarControllerHeaderObj.nBaud) / deltaHeight

        newdelta = deltaHeight * window
        r = dataOut.nHeights % window
        # nHeights = 990
        newheights = (dataOut.nHeights-r)/window
        if newheights <= 1:
            raise ValueError("filterByHeights: Too few heights. Current number of heights is %d and window is %d" %(dataOut.nHeights, window))

        if dataOut.flagDataAsBlock:
            """
            Si la data es obtenida por bloques, dimension = [nChannels, nProfiles, nHeis]
            """
            buffer = dataOut.data[:, :, 0:int(dataOut.nHeights-r)]
            buffer = buffer.reshape(dataOut.nChannels, dataOut.nProfiles, int(dataOut.nHeights/window), window)
            buffer = numpy.sum(buffer,3)

        else:
            buffer = dataOut.data[:,0:int(dataOut.nHeights-r)]
            buffer = buffer.reshape(dataOut.nChannels,int(dataOut.nHeights/window),int(window))
            buffer = numpy.sum(buffer,2)

        dataOut.data = buffer#/window
        dataOut.heightList = dataOut.heightList[0] + numpy.arange( newheights )*newdelta
        dataOut.windowOfFilter = window

        return dataOut

class setOffset(Operation):

    def run(self, dataOut, offset=None):

        if not offset:
            offset = 0.0

        newHeiRange = dataOut.heightList - offset

        dataOut.heightList = newHeiRange

        return dataOut

class setH0(Operation):

    def run(self, dataOut, h0, deltaHeight = None):

        if not deltaHeight:
            deltaHeight = dataOut.heightList[1] - dataOut.heightList[0]

        nHeights = dataOut.nHeights

        newHeiRange = h0 + numpy.arange(nHeights)*deltaHeight

        dataOut.heightList = newHeiRange

        return dataOut


class deFlip(Operation):
    def __init__(self):

        self.flip = 1

    def run(self, dataOut, channelList = []):

        data = dataOut.data.copy()

        if channelList==1:  #PARCHE #Lista de un solo canal produce error
            channelList=[1]

        

        dataOut.FlipChannels=channelList
        if dataOut.flagDataAsBlock:
            flip = self.flip
            profileList = list(range(dataOut.nProfiles))

            if not channelList:
                for thisProfile in profileList:
                    data[:,thisProfile,:] = data[:,thisProfile,:]*flip
                    flip *= -1.0
            else:
                for thisChannel in channelList:
                    if thisChannel not in dataOut.channelList:
                        continue

                    for thisProfile in profileList:
                        data[thisChannel,thisProfile,:] = data[thisChannel,thisProfile,:]*flip
                        flip *= -1.0

            self.flip = flip

        else:
            if not channelList:
                data[:,:] = data[:,:]*self.flip
            else:
                #channelList=[1]
                #print(channelList)
                #exit(1)
                for thisChannel in channelList:
                    if thisChannel not in dataOut.channelList:
                        continue

                    data[thisChannel,:] = data[thisChannel,:]*self.flip

            self.flip *= -1.

        dataOut.data = data

        #print(dataOut.data.shape)
        #print(dataOut.data[1,:,:])
        #t = [i for i in range(len(dataOut.data[1,0,:]))]

        #for i in range(300):
        #    plt.figure(figsize=(8, 4))      # Nueva figura
        #    plt.plot(t, numpy.real(dataOut.data[1, i, :]))
        #    plt.plot(t, numpy.imag(dataOut.data[1, i, :]))
        #    plt.savefig(f"/opt/results/faraday/plots/plot_1_abs_with_flip_{i:03d}.png",
        #                dpi=300,
        #                bbox_inches="tight")
        #    plt.close()

        # print("dataOut after flip", dataOut.data.shape)
        return dataOut

class deFlipHP(Operation):
    '''
    Written by R. Flores
    '''
    def __init__(self):

        self.flip = 1

    def run(self, dataOut, byHeights = False, channelList = [], HeiRangeList = None):

        data = dataOut.data.copy()

        firstHeight = HeiRangeList[0]
        lastHeight = HeiRangeList[1]+1

        #if channelList==1:  #PARCHE #Lista de un solo canal produce error
            #channelList=[1]

        dataOut.FlipChannels=channelList
        if dataOut.flagDataAsBlock:
            flip = self.flip
            profileList = list(range(dataOut.nProfiles))

            if not channelList:
                for thisProfile in profileList:
                    data[:,thisProfile,:] = data[:,thisProfile,:]*flip
                    flip *= -1.0
            else:
                for thisChannel in channelList:
                    if thisChannel not in dataOut.channelList:
                        continue
                    if not byHeights:
                        for thisProfile in profileList:
                            data[thisChannel,thisProfile,:] = data[thisChannel,thisProfile,:]*flip
                            flip *= -1.0

                    else:
                        firstHeight = HeiRangeList[0]
                        lastHeight = HeiRangeList[1]+1
                        flip = -1.0
                        data[thisChannel,:,firstHeight:lastHeight] = data[thisChannel,:,firstHeight:lastHeight]*flip


            self.flip = flip

        else:
            if not channelList:
                data[:,:] = data[:,:]*self.flip
            else:
                #channelList=[1]

                for thisChannel in channelList:
                    if thisChannel not in dataOut.channelList:
                        continue

                    if not byHeights:
                        data[thisChannel,:] = data[thisChannel,:]*flip

                    else:
                        firstHeight = HeiRangeList[0]
                        lastHeight = HeiRangeList[1]+1
                        flip = -1.0
                        data[thisChannel,firstHeight:lastHeight] = data[thisChannel,firstHeight:lastHeight]*flip

                    #data[thisChannel,:] = data[thisChannel,:]*self.flip

            self.flip *= -1.

        #print(dataOut.data[0,:12,1066+2])
        #print(dataOut.data[1,:12,1066+2])
        dataOut.data =data
        #print(dataOut.data[0,:12,1066+2])
        #print(dataOut.data[1,:12,1066+2])
        #exit(1)

        return dataOut

class setAttribute(Operation):
    '''
    Set an arbitrary attribute(s) to dataOut
    '''

    def __init__(self):

        Operation.__init__(self)
        self._ready = False

    def run(self, dataOut, **kwargs):

        for key, value in kwargs.items():

            setattr(dataOut, key, value)

        return dataOut


@MPDecorator
class printAttribute(Operation):
    '''
    Print an arbitrary attribute of dataOut
    '''

    def __init__(self):

        Operation.__init__(self)

    def run(self, dataOut, attributes):

        if isinstance(attributes, str):
            attributes = [attributes]
        for attr in attributes:
            if hasattr(dataOut, attr):
                log.log(getattr(dataOut, attr), attr)


class interpolateHeights(Operation):

    def run(self, dataOut, topLim, botLim):
        #69 al 72 para julia
        #82-84 para meteoros
        if len(numpy.shape(dataOut.data))==2:
            sampInterp = (dataOut.data[:,botLim-1] + dataOut.data[:,topLim+1])/2
            sampInterp = numpy.transpose(numpy.tile(sampInterp,(topLim-botLim + 1,1)))
            #dataOut.data[:,botLim:limSup+1] = sampInterp
            dataOut.data[:,botLim:topLim+1] = sampInterp
        else:
            nHeights = dataOut.data.shape[2]
            x = numpy.hstack((numpy.arange(botLim),numpy.arange(topLim+1,nHeights)))
            y = dataOut.data[:,:,list(range(botLim))+list(range(topLim+1,nHeights))]
            f = interpolate.interp1d(x, y, axis = 2)
            xnew = numpy.arange(botLim,topLim+1)
            ynew = f(xnew)
            dataOut.data[:,:,botLim:topLim+1]  = ynew

        return dataOut


class LagsReshape(Operation):
    '''
    Written by R. Flores
    '''
    """Operation to reshape input data into (Channels,Profiles(with same lag),Heights,Lags) and heights reconstruction.

    Parameters:
    -----------


    Example
    --------

    op = proc_unit.addOperation(name='LagsReshape')


    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

        self.buffer=None
        self.buffer_HR=None
        self.buffer_HRonelag=None

    def LagDistribution(self,dataOut):
        # Considera los primeros NSCAN perfiles
        dataOut.datapure=numpy.copy(dataOut.data[:,0:dataOut.NSCAN,:])
        self.buffer = numpy.zeros((dataOut.nChannels,
                                   int(dataOut.NSCAN/dataOut.DPL),
                                   dataOut.nHeights,dataOut.DPL),
                                  dtype='complex')

        for j in range(int(self.buffer.shape[1]/2)):
            for i in range(dataOut.DPL):
                if j+1==int(self.buffer.shape[1]/2) and i+1==dataOut.DPL:
                    self.buffer[:,2*j:,:,i]=dataOut.datapure[:,2*i+int(2*j*dataOut.DPL):,:]
                else:
                    self.buffer[:,2*j:2*(j+1),:,i]=dataOut.datapure[:,2*i+int(2*j*dataOut.DPL):2*(i+1)+int(2*j*dataOut.DPL),:]

        return self.buffer

    def HeightReconstruction(self,dataOut):

        self.buffer_HR = numpy.zeros((int(dataOut.NSCAN/dataOut.DPL),
                                   dataOut.nHeights,dataOut.DPL),
                                  dtype='complex')

        for i in range(int(dataOut.DPL)): #Only channel B
            if i==0:
                self.buffer_HR[:,:,i]=dataOut.datalags[1,:,:,i]
            else:
                self.buffer_HR[:,:,i]=self.HRonelag(dataOut,i)

        return self.buffer_HR


    def HRonelag(self,dataOut,whichlag):
        self.buffer_HRonelag = numpy.zeros((int(dataOut.NSCAN/dataOut.DPL),
                                   dataOut.nHeights),
                                  dtype='complex')
        TxLagRate = dataOut.TxLagRate
        for i in range(self.buffer_HRonelag.shape[0]): #perfil
            for j in range(dataOut.nHeights): # height
                if j+int(TxLagRate*whichlag)<dataOut.nHeights:
                    self.buffer_HRonelag[i,j]=dataOut.datalags[1,i,j+TxLagRate*whichlag,whichlag]
                else:
                    if whichlag!=10:
                        self.buffer_HRonelag[i,j]=dataOut.datalags[1,i,(j+TxLagRate*whichlag)%dataOut.nHeights,whichlag+1]
                    else:
                        if i+2<self.buffer_HRonelag.shape[0]:
                            self.buffer_HRonelag[i,j]=dataOut.datalags[1,i+2,(j+TxLagRate*whichlag)%dataOut.nHeights,0]
                        else: #i+1==self.buffer_HRonelag.shape[0]:
                            self.buffer_HRonelag[i,j]=dataOut.datalags[1,i,(j+TxLagRate*whichlag)%dataOut.nHeights,whichlag] #1, 198,64 = 1,198, 0, 10

        return self.buffer_HRonelag



    def run(self,dataOut,DPL=11,NSCAN=132, TxLagRate=2):
        # TxLagRate = 2, 2*newDeltaHeis = 30 km (para el canal 1 todos los perfiles y lags estan al inicio)

        #print(dataOut.flagNoData)
        # print(dataOut.flagNoMoreFiles)
        dataOut.TxLagRate = TxLagRate
        '''PA = dataOut.data[0,:,:]
        PB = dataOut.data[1,:,:]
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 11, figsize=(18, 6), sharex=True, sharey=True)

        for i in range(11):
            axes[0,i].plot(PA[i, :], dataOut.heightList, label=f'PA {i+1}')
            axes[0, i].set_title(f'Lag {i+1}')
            #axes[0, i].set_xscale("log")  # Log scale for y-axis
            #axes[0, i].set_xlim([0,1e+7])
            axes[1,i].plot(PB[i, :], dataOut.heightList, label=f'PB {i+1}')
            #axes[1, i].set_xscale("log")  # Log scale for y-axis
            #axes[1, i].set_xlim([0,1e+7])
            
        
        plt.tight_layout()
        plt.show()'''


        dataOut.DPL=DPL
        dataOut.NSCAN=NSCAN

        #print()
        dataOut.paramInterval=0#int(dataOut.nint*dataOut.header[7][0]*2 )
        dataOut.lat=-11.95
        dataOut.lon=-76.87
        dataOut.datalags=None

        dataOut.datalags=numpy.copy(self.LagDistribution(dataOut))
        dataOut.datalags[1,:,:,:]=self.HeightReconstruction(dataOut)

        '''
        PA = dataOut.datalags[0,5,:,:]
        PB = dataOut.datalags[1,5,:,:]
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 11, figsize=(18, 6), sharex=True, sharey=True)

        for i in range(11):
            axes[0,i].plot(PA[:, i], dataOut.heightList, label=f'PA {i+1}')
            axes[0, i].set_title(f'Lag {i+1}')
            #axes[0, i].set_xscale("log")  # Log scale for y-axis
            #axes[0, i].set_xlim([0,1e+7])
            axes[1,i].plot(PB[:, i], dataOut.heightList, label=f'PB {i+1}')
            #axes[1, i].set_xscale("log")  # Log scale for y-axis
            #axes[1, i].set_xlim([0,1e+7])
            
        
        plt.tight_layout()
        plt.show()#'''

        '''outdir = "/home/david/Documents/faraday-experiment/results/FaradayDP/Faraday_2025_01/plots_lags_reshape"
        os.makedirs(outdir, exist_ok=True)'''

        # Figura 1: dataOut.data
        '''plt.figure(figsize=(12, 5))
        plt.imshow(numpy.abs(dataOut.data[0, :, :]).T,
                aspect='auto',
                origin='lower')
        plt.xlabel("Perfil")
        plt.ylabel("Altura (bin)")
        plt.title("Canal 0 - dataOut.data")
        plt.colorbar(label="|Señal|")
        plt.savefig(f"{outdir}/dataOut_data_ch0.png",
                    dpi=300,
                    bbox_inches="tight")
        plt.close()'''

        '''ch = 1
        # Figura 2: dataOut.datalags
        for i in range(11):

            plt.figure(figsize=(8, 5))
            plt.imshow(numpy.abs(dataOut.datalags[ch, :, :, i]),
                    aspect='auto',
                    origin='lower')
            plt.xlabel("Altura (bin)")
            plt.ylabel("Perfiles")
            plt.title(f"Canal {ch} - Lag {i} - dataOut.datalags")
            plt.colorbar(label="|Señal|")
            plt.savefig(f"{outdir}/dataOut_datalags_ch{ch}_lag{i}.png",
                        dpi=300,
                        bbox_inches="tight")
            plt.close()

        print("Figuras guardadas correctamente.")'''

        return dataOut


class CrossProdDP(Operation):
    '''
    Written by R. Flores
    '''
    """Operation to calculate cross products of the Double Pulse Experiment.

    Parameters:
    -----------
    NLAG : int
        Number of lags Long Pulse.
    NRANGE : int
        Number of samples for Long Pulse.
    NCAL : int
        .*
    DPL : int
        Number of lags Double Pulse.
    NDN : int
        .*
    NDT : int
        Number of heights for Double Pulse.*
    NDP : int
        Number of heights for Double Pulse.*
    NSCAN : int
        Number of profiles when the transmitter is on.
    flags_array : intlist
        .*
    NAVG : int
        Number of blocks to be "averaged".
    nkill : int
        Number of blocks not to be considered when averaging.

    Example
    --------

    op = proc_unit.addOperation(name='CrossProdDP', optype='other')
    op.addParameter(name='NLAG', value='16', format='int')
    op.addParameter(name='NRANGE', value='0', format='int')
    op.addParameter(name='NCAL', value='0', format='int')
    op.addParameter(name='DPL', value='11', format='int')
    op.addParameter(name='NDN', value='0', format='int')
    op.addParameter(name='NDT', value='66', format='int')
    op.addParameter(name='NDP', value='66', format='int')
    op.addParameter(name='NSCAN', value='132', format='int')
    op.addParameter(name='flags_array', value='(0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300)', format='intlist')
    op.addParameter(name='NAVG', value='16', format='int')
    op.addParameter(name='nkill', value='6', format='int')

    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)
        self.bcounter=0
        self.aux=1
        self.lag_products_LP_median_estimates_aux=0

    def set_header_output(self,dataOut):

        dataOut.read_samples=len(dataOut.heightList)#int(dataOut.systemHeaderObj.nSamples/dataOut.windowOfFilter)
        padding=numpy.zeros(1,'int32')
        hsize=numpy.zeros(1,'int32')
        bufsize=numpy.zeros(1,'int32')
        nr=numpy.zeros(1,'int32')
        ngates=numpy.zeros(1,'int32') ###  ###  ### 2
        time1=numpy.zeros(1,'uint64') # pos 3
        time2=numpy.zeros(1,'uint64') # pos 4
        lcounter=numpy.zeros(1,'int32')
        groups=numpy.zeros(1,'int32')
        system=numpy.zeros(4,'int8') # pos 7
        h0=numpy.zeros(1,'float32')
        dh=numpy.zeros(1,'float32')
        ipp=numpy.zeros(1,'float32')
        process=numpy.zeros(1,'int32')
        tx=numpy.zeros(1,'int32')
        ngates1=numpy.zeros(1,'int32')  ###  ###  ### 13
        time0=numpy.zeros(1,'uint64') # pos 14
        nlags=numpy.zeros(1,'int32')
        nlags1=numpy.zeros(1,'int32')
        txb=numpy.zeros(1,'float32')   ###  ###  ### 17
        time3=numpy.zeros(1,'uint64') # pos 18
        time4=numpy.zeros(1,'uint64') # pos 19
        h0_=numpy.zeros(1,'float32')
        dh_=numpy.zeros(1,'float32')
        ipp_=numpy.zeros(1,'float32')
        txa_=numpy.zeros(1,'float32')
        pad=numpy.zeros(100,'int32')
        nbytes=numpy.zeros(1,'int32')
        limits=numpy.zeros(1,'int32')
        ngroups=numpy.zeros(1,'int32') ###  ###  ### 27

        dataOut.header=[hsize,bufsize,nr,ngates,time1,time2,
                lcounter,groups,system,h0,dh,ipp,
                process,tx,ngates1,padding,time0,nlags,
                nlags1,padding,txb,time3,time4,h0_,dh_,
                ipp_,txa_,pad,nbytes,limits,padding,ngroups]


        #dataOut.header[1][0]=81864
        dataOut.FirstHeight=int(dataOut.heightList[0])
        dataOut.MAXNRANGENDT=max(dataOut.NRANGE,dataOut.NDT)
        dataOut.header[3][0]=max(dataOut.NRANGE,dataOut.NDT)
        dataOut.header[7][0]=dataOut.NAVG
        dataOut.header[9][0]=int(dataOut.heightList[0])
        dataOut.header[10][0]=dataOut.DH
        dataOut.header[17][0]=dataOut.DPL
        dataOut.header[18][0]=dataOut.NLAG
        #self.header[5][0]=0
        dataOut.header[15][0]=dataOut.NDP
        dataOut.header[2][0]=dataOut.NR


    def get_products_cabxys(self,dataOut):

        if self.aux==1:
            self.set_header_output(dataOut)
            self.aux=0

        dataOut.lags_array=[x / dataOut.DH for x in dataOut.flags_array]
        self.cax=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.cay=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.cbx=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.cby=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.cax2=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.cay2=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.cbx2=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.cby2=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.caxbx=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.caxby=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.caybx=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.cayby=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.caxay=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.cbxby=numpy.zeros((dataOut.NDP,dataOut.DPL,2))

        for i in range(2):
            for j in range(dataOut.NDP):
                for k in range(int(dataOut.NSCAN/2)):
                    n=k%dataOut.DPL
                    ax=dataOut.data[0,2*k+i,j].real
                    ay=dataOut.data[0,2*k+i,j].imag
                    if j+dataOut.lags_array[n]<dataOut.NDP:
                        bx=dataOut.data[1,2*k+i,j+int(dataOut.lags_array[n])].real
                        by=dataOut.data[1,2*k+i,j+int(dataOut.lags_array[n])].imag
                    else:
                        if k+1<int(dataOut.NSCAN/2):
                            bx=dataOut.data[1,2*(k+1)+i,(dataOut.NRANGE+dataOut.NCAL+j+int(dataOut.lags_array[n]))%dataOut.NDP].real
                            by=dataOut.data[1,2*(k+1)+i,(dataOut.NRANGE+dataOut.NCAL+j+int(dataOut.lags_array[n]))%dataOut.NDP].imag

                        if k+1==int(dataOut.NSCAN/2):
                            bx=dataOut.data[1,2*k+i,(dataOut.NRANGE+dataOut.NCAL+j+int(dataOut.lags_array[n]))%dataOut.NDP].real
                            by=dataOut.data[1,2*k+i,(dataOut.NRANGE+dataOut.NCAL+j+int(dataOut.lags_array[n]))%dataOut.NDP].imag

                    if(k<dataOut.DPL):
                        self.cax[j][n][i]=ax
                        self.cay[j][n][i]=ay
                        self.cbx[j][n][i]=bx
                        self.cby[j][n][i]=by
                        self.cax2[j][n][i]=ax*ax
                        self.cay2[j][n][i]=ay*ay
                        self.cbx2[j][n][i]=bx*bx
                        self.cby2[j][n][i]=by*by
                        self.caxbx[j][n][i]=ax*bx
                        self.caxby[j][n][i]=ax*by
                        self.caybx[j][n][i]=ay*bx
                        self.cayby[j][n][i]=ay*by
                        self.caxay[j][n][i]=ax*ay
                        self.cbxby[j][n][i]=bx*by
                    else:
                        self.cax[j][n][i]+=ax
                        self.cay[j][n][i]+=ay
                        self.cbx[j][n][i]+=bx
                        self.cby[j][n][i]+=by
                        self.cax2[j][n][i]+=ax*ax
                        self.cay2[j][n][i]+=ay*ay
                        self.cbx2[j][n][i]+=bx*bx
                        self.cby2[j][n][i]+=by*by
                        self.caxbx[j][n][i]+=ax*bx
                        self.caxby[j][n][i]+=ax*by
                        self.caybx[j][n][i]+=ay*bx
                        self.cayby[j][n][i]+=ay*by
                        self.caxay[j][n][i]+=ax*ay
                        self.cbxby[j][n][i]+=bx*by


    def medi(self,data_navg,NAVG,nkill):
        sorts=sorted(data_navg)
        rsorts=numpy.arange(NAVG)
        result=0.0
        for k in range(NAVG):
            if k>=nkill/2 and k<NAVG-nkill/2:
                result+=sorts[k]*float(NAVG)/(float(NAVG-nkill))
        return result


    def get_dc(self,dataOut):
        if self.bcounter==0:
            dataOut.dc=numpy.zeros(dataOut.NR,dtype='complex64')
    def cabxys_navg(self,dataOut):


        dataOut.header[5][0]=dataOut.TimeBlockSeconds

        dataOut.LastAVGDate=dataOut.TimeBlockSeconds

        if self.bcounter==0:
            dataOut.FirstAVGDate=dataOut.TimeBlockSeconds
            dataOut.header[4][0]=dataOut.header[5][0]#firsttimeofNAVG
            if dataOut.CurrentBlock==1:
                dataOut.FirstBlockDate=dataOut.TimeBlockSeconds
                dataOut.header[16][0]=dataOut.header[5][0]#FirsTimeOfTotalBlocks

            self.cax_navg=[]
            self.cay_navg=[]
            self.cbx_navg=[]
            self.cby_navg=[]
            self.cax2_navg=[]
            self.cay2_navg=[]
            self.cbx2_navg=[]
            self.cby2_navg=[]
            self.caxbx_navg=[]
            self.caxby_navg=[]
            self.caybx_navg=[]
            self.cayby_navg=[]
            self.caxay_navg=[]
            self.cbxby_navg=[]

            dataOut.noisevector=numpy.zeros((dataOut.MAXNRANGENDT,dataOut.NR,dataOut.NAVG),'float32')

            dataOut.noisevector_=numpy.zeros((dataOut.read_samples,dataOut.NR,dataOut.NAVG),'float32')

        self.noisevectorizer(dataOut.NSCAN,dataOut.nProfiles,dataOut.NR,dataOut.MAXNRANGENDT,dataOut.noisevector,dataOut.data,dataOut.dc)   #30/03/2020

        self.cax_navg.append(self.cax)
        self.cay_navg.append(self.cay)
        self.cbx_navg.append(self.cbx)
        self.cby_navg.append(self.cby)
        self.cax2_navg.append(self.cax2)
        self.cay2_navg.append(self.cay2)
        self.cbx2_navg.append(self.cbx2)
        self.cby2_navg.append(self.cby2)
        self.caxbx_navg.append(self.caxbx)
        self.caxby_navg.append(self.caxby)
        self.caybx_navg.append(self.caybx)
        self.cayby_navg.append(self.cayby)
        self.caxay_navg.append(self.caxay)
        self.cbxby_navg.append(self.cbxby)
        self.bcounter+=1

    def noise_estimation4x_DP(self,dataOut):
        if self.bcounter==dataOut.NAVG:
            dataOut.noise_final=numpy.zeros(dataOut.NR,'float32')
            snoise=numpy.zeros((dataOut.NR,dataOut.NAVG),'float32')
            nvector1=numpy.zeros((dataOut.NR,dataOut.NAVG,dataOut.MAXNRANGENDT),'float32')
            for i in range(dataOut.NR):
                dataOut.noise_final[i]=0.0
                for k in range(dataOut.NAVG):
                    snoise[i][k]=0.0
                    for j in range(dataOut.MAXNRANGENDT):
                        nvector1[i][k][j]= dataOut.noisevector[j][i][k];
                    snoise[i][k]=self.noise_hs4x(dataOut.MAXNRANGENDT, nvector1[i][k])
                dataOut.noise_final[i]=self.noise_hs4x(dataOut.NAVG, snoise[i])

    def kabxys(self,dataOut):

        if self.bcounter==dataOut.NAVG:

            dataOut.flagNoData =  False

            self.kax=numpy.zeros((dataOut.NDP,dataOut.DPL,2),'float32')
            self.kay=numpy.zeros((dataOut.NDP,dataOut.DPL,2),'float32')
            self.kbx=numpy.zeros((dataOut.NDP,dataOut.DPL,2),'float32')
            self.kby=numpy.zeros((dataOut.NDP,dataOut.DPL,2),'float32')
            self.kax2=numpy.zeros((dataOut.NDP,dataOut.DPL,2),'float32')
            self.kay2=numpy.zeros((dataOut.NDP,dataOut.DPL,2),'float32')
            self.kbx2=numpy.zeros((dataOut.NDP,dataOut.DPL,2),'float32')
            self.kby2=numpy.zeros((dataOut.NDP,dataOut.DPL,2),'float32')
            self.kaxbx=numpy.zeros((dataOut.NDP,dataOut.DPL,2),'float32')
            self.kaxby=numpy.zeros((dataOut.NDP,dataOut.DPL,2),'float32')
            self.kaybx=numpy.zeros((dataOut.NDP,dataOut.DPL,2),'float32')
            self.kayby=numpy.zeros((dataOut.NDP,dataOut.DPL,2),'float32')
            self.kaxay=numpy.zeros((dataOut.NDP,dataOut.DPL,2),'float32')
            self.kbxby=numpy.zeros((dataOut.NDP,dataOut.DPL,2),'float32')
            auxx = numpy.array(self.cax2_navg)
            auxx2 = numpy.array(self.cay2_navg)
            auxx3 = numpy.array(self.cbx2_navg)
            auxx4 = numpy.array(self.cby2_navg)

            pa1 = 20
            pa2 = 10
            #aux1 = dataOut.kax2[pa1,pa2,0]+dataOut.kax2[pa1,pa2,1]+dataOut.kay2[pa1,pa2,0]+dataOut.kay2[pa1,pa2,1]
            #aux2 = dataOut.kbx2[pa1,pa2,0]+dataOut.kbx2[pa1,pa2,1]+dataOut.kby2[pa1,pa2,0]+dataOut.kby2[pa1,pa2,1]

            #print(aux1)#*numpy.conjugate(aux1))
            #print(aux2)#*numpy.conjugate(aux2))
            #print(auxx.shape)
            '''
            print(numpy.sum(auxx2[:,pa1,pa2,0]+auxx2[:,pa1,pa2,1]))#auxx[:,pa1,pa2,0]+auxx[:,pa1,pa2,1]))#+auxx2[:,pa1,pa2,0]+auxx2[:,pa1,pa2,1]))
            print(numpy.sum(auxx3[:,pa1,pa2,0]+auxx3[:,pa1,pa2,1]+auxx4[:,pa1,pa2,0]+auxx4[:,pa1,pa2,1]))
            #print(self.cax_navg[:,53,0,1])
            #exit(1)
            '''

            for i in range(self.cax_navg[0].shape[0]):
                        for j in range(self.cax_navg[0].shape[1]):
                            for k in range(self.cax_navg[0].shape[2]):
                                data_navg=[item[i,j,k] for item in self.cax_navg]
                                self.kax[i,j,k]=self.medi(data_navg,dataOut.NAVG,dataOut.nkill)
                                data_navg=[item[i,j,k] for item in self.cay_navg]
                                self.kay[i,j,k]=self.medi(data_navg,dataOut.NAVG,dataOut.nkill)
                                data_navg=[item[i,j,k] for item in self.cbx_navg]
                                self.kbx[i,j,k]=self.medi(data_navg,dataOut.NAVG,dataOut.nkill)
                                data_navg=[item[i,j,k] for item in self.cby_navg]
                                self.kby[i,j,k]=self.medi(data_navg,dataOut.NAVG,dataOut.nkill)
                                data_navg=[item[i,j,k] for item in self.cax2_navg]
                                self.kax2[i,j,k]=self.medi(data_navg,dataOut.NAVG,dataOut.nkill)
                                data_navg=[item[i,j,k] for item in self.cay2_navg]
                                self.kay2[i,j,k]=self.medi(data_navg,dataOut.NAVG,dataOut.nkill)
                                data_navg=[item[i,j,k] for item in self.cbx2_navg]
                                self.kbx2[i,j,k]=self.medi(data_navg,dataOut.NAVG,dataOut.nkill)
                                data_navg=[item[i,j,k] for item in self.cby2_navg]
                                self.kby2[i,j,k]=self.medi(data_navg,dataOut.NAVG,dataOut.nkill)
                                data_navg=[item[i,j,k] for item in self.caxbx_navg]
                                self.kaxbx[i,j,k]=self.medi(data_navg,dataOut.NAVG,dataOut.nkill)
                                data_navg=[item[i,j,k] for item in self.caxby_navg]
                                self.kaxby[i,j,k]=self.medi(data_navg,dataOut.NAVG,dataOut.nkill)
                                data_navg=[item[i,j,k] for item in self.caybx_navg]
                                self.kaybx[i,j,k]=self.medi(data_navg,dataOut.NAVG,dataOut.nkill)
                                data_navg=[item[i,j,k] for item in self.cayby_navg]
                                self.kayby[i,j,k]=self.medi(data_navg,dataOut.NAVG,dataOut.nkill)
                                data_navg=[item[i,j,k] for item in self.caxay_navg]
                                self.kaxay[i,j,k]=self.medi(data_navg,dataOut.NAVG,dataOut.nkill)
                                data_navg=[item[i,j,k] for item in self.cbxby_navg]
                                self.kbxby[i,j,k]=self.medi(data_navg,dataOut.NAVG,dataOut.nkill)


            dataOut.kax=self.kax
            dataOut.kay=self.kay
            dataOut.kbx=self.kbx
            dataOut.kby=self.kby
            dataOut.kax2=self.kax2
            dataOut.kay2=self.kay2
            dataOut.kbx2=self.kbx2
            dataOut.kby2=self.kby2
            dataOut.kaxbx=self.kaxbx
            dataOut.kaxby=self.kaxby
            dataOut.kaybx=self.kaybx
            dataOut.kayby=self.kayby
            dataOut.kaxay=self.kaxay
            dataOut.kbxby=self.kbxby

            #FindMe
            pa1 = 20
            pa2 = 0
            '''
            aux1 = dataOut.kax2[pa1,pa2,0]+dataOut.kax2[pa1,pa2,1]#+dataOut.kay2[pa1,pa2,0]+dataOut.kay2[pa1,pa2,1]
            aux2 = dataOut.kbx2[pa1,pa2,0]+dataOut.kbx2[pa1,pa2,1]+dataOut.kby2[pa1,pa2,0]+dataOut.kby2[pa1,pa2,1]
            aux3 = dataOut.kay2[pa1,pa2,0]+dataOut.kay2[pa1,pa2,1]
            aux4 = dataOut.kax2[pa1,pa2,0]+dataOut.kax2[pa1,pa2,1]+dataOut.kay2[pa1,pa2,0]+dataOut.kay2[pa1,pa2,1]

            print(aux1)#*numpy.conjugate(aux1))
            print(aux3)
            print("sum",aux4)
            print(aux2)#*numpy.conjugate(aux2))
            '''
            '''
            aux1 = dataOut.kaxbx[pa1,pa2,0]+dataOut.kaxbx[pa1,pa2,1]+dataOut.kayby[pa1,pa2,0]+dataOut.kayby[pa1,pa2,1]
            aux2 = dataOut.kaybx[pa1,pa2,0]+dataOut.kaybx[pa1,pa2,1]-dataOut.kaxby[pa1,pa2,0]-dataOut.kaxby[pa1,pa2,1]
            print(aux1)
            print(aux2)
            exit(1)
            '''

            #print(dataOut.kax[53,0,0])
            #exit(1)

            self.bcounter=0

            dataOut.crossprods=numpy.zeros((3,4,numpy.shape(dataOut.kax)[0],numpy.shape(dataOut.kax)[1],numpy.shape(dataOut.kax)[2]))

            dataOut.crossprods[0]=[dataOut.kax,dataOut.kay,dataOut.kbx,dataOut.kby]
            dataOut.crossprods[1]=[dataOut.kax2,dataOut.kay2,dataOut.kbx2,dataOut.kby2]
            dataOut.crossprods[2]=[dataOut.kaxay,dataOut.kbxby,dataOut.kaxbx,dataOut.kaxby]
            dataOut.data_for_RTI_DP=numpy.zeros((3,dataOut.NDP))
            dataOut.data_for_RTI_DP[0],dataOut.data_for_RTI_DP[1],dataOut.data_for_RTI_DP[2]=self.RTI_COLUMN(dataOut.kax2,dataOut.kay2,dataOut.kbx2,dataOut.kby2,dataOut.kaxbx,dataOut.kayby,dataOut.kaybx,dataOut.kaxby, dataOut.NDP)



    def RTI_COLUMN(self,kax2,kay2,kbx2,kby2,kaxbx,kayby,kaybx,kaxby, NDP):
        x00=numpy.zeros(NDP,dtype='float32')
        x01=numpy.zeros(NDP,dtype='float32')
        x02=numpy.zeros(NDP,dtype='float32')
        for j in range(2):# first couple lags
            for k in range(2): #flip
                for i in range(NDP): #
                    fx=numpy.sqrt((kaxbx[i,j,k]+kayby[i,j,k])**2+(kaybx[i,j,k]-kaxby[i,j,k])**2)
                    x00[i]=x00[i]+(kax2[i,j,k]+kay2[i,j,k])
                    x01[i]=x01[i]+(kbx2[i,j,k]+kby2[i,j,k])
                    x02[i]=x02[i]+fx

                    x00[i]=10.0*numpy.log10(x00[i]/512.)
                    x01[i]=10.0*numpy.log10(x01[i]/512.)
                    x02[i]=10.0*numpy.log10(x02[i])
        return x02,x00,x01


    def noisevectorizer(self,NSCAN,nProfiles,NR,MAXNRANGENDT,noisevector,data,dc):

        rnormalizer= 1./(float(nProfiles - NSCAN))
        #rnormalizer= float(NSCAN)/((float(nProfiles - NSCAN))*float(MAXNRANGENDT))
        for i in range(NR):
            for j in range(MAXNRANGENDT):
                for k in range(NSCAN,nProfiles):
                    #TODO:integrate just 2nd quartile gates
                    if k==NSCAN:
                        noisevector[j][i][self.bcounter]=(abs(data[i][k][j]-dc[i])**2)*rnormalizer
                    else:
                        noisevector[j][i][self.bcounter]+=(abs(data[i][k][j]-dc[i])**2)*rnormalizer




    def  noise_hs4x(self, ndatax, datax):
        divider=10#divider was originally 10
        noise=0.0
        data=numpy.zeros(ndatax,'float32')
        ndata1=int(ndatax/4)
        ndata2=int(2.5*(ndatax/4.))
        ndata=int(ndata2-ndata1)
        sorts=sorted(datax)

        for k in range(ndata2): # select just second quartile
            data[k]=sorts[k+ndata1]
        nums_min= int(ndata/divider)
        if(int(ndata/divider)> 2):
            nums_min= int(ndata/divider)
        else:
            nums_min=2
        sump=0.0
        sumq=0.0
        j=0
        cont=1
        while ( (cont==1) and (j<ndata)):
            sump+=data[j]
            sumq+= data[j]*data[j]
            j=j+1
            if (j> nums_min):
                rtest= float(j/(j-1)) +1.0/ndata
                if( (sumq*j) > (rtest*sump*sump ) ):
                    j=j-1
                    sump-= data[j]
                    sumq-=data[j]*data[j]
                    cont= 0
        noise= (sump/j)

        return noise



    def run(self, dataOut, NLAG=16, NRANGE=0, NCAL=0, DPL=11,
        NDN=0, NDT=66, NDP=66, NSCAN=132,
        flags_array=(0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300), NAVG=16, nkill=6, **kwargs):

        dataOut.NLAG=NLAG
        dataOut.NR=len(dataOut.channelList)
        dataOut.NRANGE=NRANGE
        dataOut.NCAL=NCAL
        dataOut.DPL=DPL
        dataOut.NDN=NDN
        dataOut.NDT=NDT
        dataOut.NDP=NDP
        dataOut.NSCAN=NSCAN
        dataOut.DH=dataOut.heightList[1]-dataOut.heightList[0]
        dataOut.H0=int(dataOut.heightList[0])
        dataOut.flags_array=flags_array
        dataOut.NAVG=NAVG
        dataOut.nkill=nkill
        dataOut.flagNoData =  True

        self.get_dc(dataOut)
        self.get_products_cabxys(dataOut)
        self.cabxys_navg(dataOut)
        self.noise_estimation4x_DP(dataOut)
        self.kabxys(dataOut)

        return dataOut



class IntegrationDP(Operation):
    '''
    Written by R. Flores
    '''
    """Operation to integrate the Double Pulse data.

    Parameters:
    -----------
    nint : int
        Number of integrations.

    Example
    --------

    op = proc_unit.addOperation(name='IntegrationDP', optype='other')
    op.addParameter(name='nint', value='30', format='int')

    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

        self.counter=0
        self.aux=0
        self.init_time=None

    def integration_for_double_pulse(self,dataOut):

        if self.aux==1:

            dataOut.TimeBlockSeconds_for_dp_power=dataOut.utctime
            dataOut.bd_time=gmtime(dataOut.TimeBlockSeconds_for_dp_power)
            dataOut.year=dataOut.bd_time.tm_year+(dataOut.bd_time.tm_yday-1)/364.0
            dataOut.ut_Faraday=dataOut.bd_time.tm_hour+dataOut.bd_time.tm_min/60.0+dataOut.bd_time.tm_sec/3600.0
            self.aux=0

        if self.counter==0:

            tmpx=numpy.zeros((dataOut.NDP,dataOut.DPL,2),'float32')
            dataOut.kabxys_integrated=[tmpx,tmpx,tmpx,tmpx,tmpx,tmpx,tmpx,tmpx,tmpx,tmpx,tmpx,tmpx,tmpx,tmpx]
            self.init_time=dataOut.utctime

        if self.counter < dataOut.nint:

            dataOut.final_cross_products=[dataOut.kax,dataOut.kay,dataOut.kbx,dataOut.kby,dataOut.kax2,dataOut.kay2,dataOut.kbx2,dataOut.kby2,dataOut.kaxbx,dataOut.kaxby,dataOut.kaybx,dataOut.kayby,dataOut.kaxay,dataOut.kbxby]

            for ind in range(len(dataOut.kabxys_integrated)): #final cross products
                dataOut.kabxys_integrated[ind]=dataOut.kabxys_integrated[ind]+dataOut.final_cross_products[ind]

            self.counter+=1

            if self.counter==dataOut.nint-1:
                self.aux=1

            if self.counter==dataOut.nint:
                dataOut.flagNoData=False
                pa1 = 20
                pa2 = 10
                '''
                print(32*(dataOut.kabxys_integrated[4][pa1,pa2,0]+dataOut.kabxys_integrated[5][pa1,pa2,0]+dataOut.kabxys_integrated[4][pa1,pa2,1]+dataOut.kabxys_integrated[5][pa1,pa2,1]))
                print(32*(dataOut.kabxys_integrated[6][pa1,pa2,0]+dataOut.kabxys_integrated[7][pa1,pa2,0]+dataOut.kabxys_integrated[6][pa1,pa2,1]+dataOut.kabxys_integrated[7][pa1,pa2,1]))

                exit(1)
                '''
                dataOut.utctime=self.init_time
                self.counter=0
                '''
                print(dataOut.kabxys_integrated[8][53,6,0]+dataOut.kabxys_integrated[11][53,6,0])
                print(dataOut.kabxys_integrated[8][53,9,0]+dataOut.kabxys_integrated[11][53,9,0])
                exit(1)
                '''


    def run(self,dataOut,nint=20):

        dataOut.flagNoData=True
        dataOut.nint=nint
        dataOut.paramInterval=0#int(dataOut.nint*dataOut.header[7][0]*2 )
        dataOut.lat=-11.95
        dataOut.lon=-76.87

        self.integration_for_double_pulse(dataOut)


        return dataOut


class SumFlips(Operation):
    '''
    Written by R. Flores
    '''
    """Operation to sum the flip and unflip part of certain cross products of the Double Pulse.

    Parameters:
    -----------
    None

    Example
    --------

    op = proc_unit.addOperation(name='SumFlips', optype='other')

    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)


    def rint2DP(self,dataOut):

        dataOut.rnint2=numpy.zeros(dataOut.DPL,'float32')

        for l in range(dataOut.DPL):

            dataOut.rnint2[l]=1.0/(dataOut.nint*dataOut.NAVG*12.0)


    def SumLags(self,dataOut):

        for l in range(dataOut.DPL):
            '''
            if l == 10:
                print(32*(dataOut.kabxys_integrated[4][20,10,0]+dataOut.kabxys_integrated[5][20,10,0]+dataOut.kabxys_integrated[4][20,10,1]+dataOut.kabxys_integrated[5][20,10,1]))
                print(32*(dataOut.kabxys_integrated[6][20,10,0]+dataOut.kabxys_integrated[7][20,10,0]+dataOut.kabxys_integrated[6][20,10,1]+dataOut.kabxys_integrated[7][20,10,1]))
                '''
            dataOut.kabxys_integrated[4][:,l,0]=(dataOut.kabxys_integrated[4][:,l,0]+dataOut.kabxys_integrated[4][:,l,1])*dataOut.rnint2[l]
            dataOut.kabxys_integrated[5][:,l,0]=(dataOut.kabxys_integrated[5][:,l,0]+dataOut.kabxys_integrated[5][:,l,1])*dataOut.rnint2[l]
            dataOut.kabxys_integrated[6][:,l,0]=(dataOut.kabxys_integrated[6][:,l,0]+dataOut.kabxys_integrated[6][:,l,1])*dataOut.rnint2[l]
            dataOut.kabxys_integrated[7][:,l,0]=(dataOut.kabxys_integrated[7][:,l,0]+dataOut.kabxys_integrated[7][:,l,1])*dataOut.rnint2[l]

            dataOut.kabxys_integrated[8][:,l,0]=(dataOut.kabxys_integrated[8][:,l,0]-dataOut.kabxys_integrated[8][:,l,1])*dataOut.rnint2[l]
            dataOut.kabxys_integrated[9][:,l,0]=(dataOut.kabxys_integrated[9][:,l,0]-dataOut.kabxys_integrated[9][:,l,1])*dataOut.rnint2[l]
            dataOut.kabxys_integrated[10][:,l,0]=(dataOut.kabxys_integrated[10][:,l,0]-dataOut.kabxys_integrated[10][:,l,1])*dataOut.rnint2[l]
            dataOut.kabxys_integrated[11][:,l,0]=(dataOut.kabxys_integrated[11][:,l,0]-dataOut.kabxys_integrated[11][:,l,1])*dataOut.rnint2[l]
            '''
            if l == 10:
                print(32*(dataOut.kabxys_integrated[4][20,10,0]+dataOut.kabxys_integrated[5][20,10,0]))
                print(32*(dataOut.kabxys_integrated[6][20,10,0]+dataOut.kabxys_integrated[7][20,10,0]))
                exit(1)
                '''
    def run(self,dataOut):

        self.rint2DP(dataOut)
        self.SumLags(dataOut)


        return dataOut


class FlagBadHeights(Operation):
    '''
    Written by R. Flores
    '''
    """Operation to flag bad heights (bad data) of the Double Pulse.

    Parameters:
    -----------
    None

    Example
    --------

    op = proc_unit.addOperation(name='FlagBadHeights', optype='other')

    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

    def run(self,dataOut):

        dataOut.ibad=numpy.zeros((dataOut.NDP,dataOut.DPL),'int32')

        for j in range(dataOut.NDP):
            for l in range(dataOut.DPL):
                ip1=j+dataOut.NDP*(0+2*l)

                if( (dataOut.kabxys_integrated[5][j,l,0] <= 0.) or (dataOut.kabxys_integrated[4][j,l,0] <= 0.) or (dataOut.kabxys_integrated[7][j,l,0] <= 0.) or (dataOut.kabxys_integrated[6][j,l,0] <= 0.)):
                    dataOut.ibad[j][l]=1
                else:
                    dataOut.ibad[j][l]=0
        #print("dataOut.ibad",dataOut.ibad)
        return dataOut

class FlagBadHeightsSpectra(Operation):
    '''
    Written by R. Flores
    '''
    """Operation to flag bad heights (bad data) of the Double Pulse.

    Parameters:
    -----------
    None

    Example
    --------

    op = proc_unit.addOperation(name='FlagBadHeightsSpectra', optype='other')

    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

    def run(self,dataOut):

        dataOut.ibad=numpy.zeros((dataOut.NDP,dataOut.DPL),'int32')

        for j in range(dataOut.NDP):
            for l in range(dataOut.DPL):
                ip1=j+dataOut.NDP*(0+2*l)

                if( (dataOut.kabxys_integrated[4][j,l,0] <= 0.) or (dataOut.kabxys_integrated[6][j,l,0] <= 0.)):
                    dataOut.ibad[j][l]=1
                else:
                    dataOut.ibad[j][l]=0

        return dataOut

class CleanCohEchoes(Operation):
    '''
    Written by R. Flores
    '''
    """Operation to clean coherent echoes.

    Parameters:
    -----------
    None

    Example
    --------

    op = proc_unit.addOperation(name='CleanCohEchoes')

    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

    def remove_coh(self,pow):
        #print("pow inside: ",pow)
        #print(pow.shape)
        q75,q25 = numpy.percentile(pow,[75,25],axis=0)
        #print(q75,q25)
        intr_qr = q75-q25

        max = q75+(1.5*intr_qr)
        min = q25-(1.5*intr_qr)

        pow[pow > max] = numpy.nan

        #print("Max: ",max)
        #print("Min: ",min)

        return pow

    def mad_based_outlier_V0(self, points, thresh=3.5):
        #print("points: ",points)
        if len(points.shape) == 1:
            points = points[:,None]
        median = numpy.nanmedian(points, axis=0)
        diff = numpy.nansum((points - median)**2, axis=-1)
        diff = numpy.sqrt(diff)
        med_abs_deviation = numpy.nanmedian(diff)

        modified_z_score = 0.6745 * diff / med_abs_deviation
        #print(modified_z_score)
        return modified_z_score > thresh

    def mad_based_outlier(self, points, thresh=3.5):

        median = numpy.nanmedian(points)
        diff = (points - median)**2
        diff = numpy.sqrt(diff)
        med_abs_deviation = numpy.nanmedian(diff)

        modified_z_score = 0.6745 * diff / med_abs_deviation

        return modified_z_score > thresh

    def removeSpreadF_V0(self,dataOut):
        for i in range(11):
            print("BEFORE Chb: ",i,dataOut.kabxys_integrated[6][:,i,0])
        #exit(1)

        #Removing echoes greater than 35 dB
        maxdB = 35 #DEBERÍA SER NOISE+ALGO!!!!!!!!!!!!!!!!!!!!!!
        #print(dataOut.kabxys_integrated[6][:,0,0])
        data = numpy.copy(10*numpy.log10(dataOut.kabxys_integrated[6][:,0,0])) #Lag0 ChB
        #print(data)
        for i in range(12,data.shape[0]):
            #for j in range(data.shape[1]):
            if data[i]>maxdB:
                dataOut.kabxys_integrated[4][i-2:i+3,:,0] = numpy.nan #Debido a que estos ecos son intensos, se
                dataOut.kabxys_integrated[6][i-2:i+3,:,0] = numpy.nan #remueve además dos muestras antes y después
                #dataOut.kabxys_integrated[4][i-1,:,0] = numpy.nan
                #dataOut.kabxys_integrated[6][i-1,:,0] = numpy.nan
                #dataOut.kabxys_integrated[4][i+1,:,0] = numpy.nan
                #dataOut.kabxys_integrated[6][i+1,:,0] = numpy.nan
                dataOut.flagSpreadF = True
                print("Removing Threshold",i)
                #print("i: ",i)

        #print("BEFORE Chb: ",dataOut.kabxys_integrated[6][:,0,0])
        #exit(1)

        #Removing outliers from the profile
        nlag = 9
        minHei = 180
        #maxHei = 600
        maxHei = 525
        inda = numpy.where(dataOut.heightList >= minHei)
        indb = numpy.where(dataOut.heightList <= maxHei)
        minIndex = inda[0][0]
        maxIndex = indb[0][-1]
        #l0 = 0
        #print("BEFORE Cha: ",dataOut.kabxys_integrated[4][:,l0,0])
        #print("BEFORE Chb: ",dataOut.kabxys_integrated[6][:,l0,0])
        #exit(1)
        #'''
        l0 = 0
        #print("BEFORE Cha: ",dataOut.kabxys_integrated[4][:,l0,0])
        #print("BEFORE Chb: ",dataOut.kabxys_integrated[6][:,l0,0])

        import matplotlib.pyplot as plt
        for i in range(l0,l0+11):
            plt.plot(dataOut.kabxys_integrated[6][:,i,0],dataOut.heightList,label='{}'.format(i))
        #plt.xlim(1.e5,1.e8)
        plt.legend()
        plt.xlim(0,2000)
        plt.show()
        #'''
        #dataOut.kabxys_integrated[4][minIndex:,:,0] = self.remove_coh(dataOut.kabxys_integrated[4][minIndex:,:,0
        outliers_IDs = []
        '''
        for lag in range(11):
            outliers = self.mad_based_outlier(dataOut.kabxys_integrated[4][minIndex:,lag,0], thresh=3.)
            #print("Outliers: ",outliers)
            #indexes.append(outliers.nonzero())
            #numpy.concatenate((outliers))
            #dataOut.kabxys_integrated[4][minIndex:,lag,0][outliers == True] = numpy.nan
            outliers_IDs=numpy.append(outliers_IDs,outliers.nonzero())
            '''
        for lag in range(11):
            #outliers = self.mad_based_outlier(dataOut.kabxys_integrated[6][minIndex:maxIndex,lag,0], thresh=2.)
            outliers = self.mad_based_outlier(dataOut.kabxys_integrated[6][minIndex:maxIndex,lag,0])
            outliers_IDs=numpy.append(outliers_IDs,outliers.nonzero())
        #print(outliers_IDs)
        #exit(1)
        if outliers_IDs != []:
            outliers_IDs=numpy.array(outliers_IDs)
            outliers_IDs=outliers_IDs.ravel()
            outliers_IDs=outliers_IDs.astype(numpy.dtype('int64'))

            (uniq, freq) = (numpy.unique(outliers_IDs, return_counts=True))
            aux_arr = numpy.column_stack((uniq,freq))
            #print("repetitions: ",aux_arr)

        #if aux_arr != []:
            final_index = []
            for i in range(aux_arr.shape[0]):
                if aux_arr[i,1] >= 10:
                    final_index.append(aux_arr[i,0])

            if final_index != [] and len(final_index) > 1:
                final_index += minIndex
                #print("final_index: ",final_index)
                following_index = final_index[-1]+1 #Remove following index to ensure we remove remaining SpreadF
                previous_index = final_index[0]-1 #Remove previous index to ensure we remove remaning SpreadF
                final_index = numpy.concatenate(([previous_index],final_index,[following_index]))
                final_index = numpy.unique(final_index) #If there was only one outlier
                #print("final_index: ",final_index)
                #exit(1)
                dataOut.kabxys_integrated[4][final_index,:,0] = numpy.nan
                dataOut.kabxys_integrated[6][final_index,:,0] = numpy.nan

                dataOut.flagSpreadF = True

        #print(final_index+minIndex)
        #print(outliers_IDs)
        #exit(1)
        #print("flagSpreadF",dataOut.flagSpreadF)

        '''
        for lag in range(11):
            #print("Lag: ",lag)
            outliers = self.mad_based_outlier(dataOut.kabxys_integrated[6][minIndex:maxIndex,lag,0], thresh=2.)
            dataOut.kabxys_integrated[6][minIndex:maxIndex,lag,0][outliers == True] = numpy.nan
            '''
        #dataOut.kabxys_integrated[4][minIndex:,:,0] = self.remove_coh(dataOut.kabxys_integrated[4][minIndex:,:,0])
        '''
        import matplotlib.pyplot as plt
        for i in range(11):
            plt.plot(dataOut.kabxys_integrated[6][:,i,0],dataOut.heightList,label='{}'.format(i))
        plt.xlim(0,2000)
        plt.legend()
        plt.grid()
        plt.show()
        '''
        '''
        for nlag in range(11):
            print("BEFORE",dataOut.kabxys_integrated[6][:,nlag,0])
        #exit(1)
        '''
        #dataOut.kabxys_integrated[6][minIndex:,:,0] = self.remove_coh(dataOut.kabxys_integrated[6][minIndex:,:,0])


        '''
        for nlag in range(11):
            print("AFTER",dataOut.kabxys_integrated[6][:,nlag,0])
        exit(1)
        '''
        #print("AFTER",dataOut.kabxys_integrated[4][33,:,0])
        #print("AFTER",dataOut.kabxys_integrated[6][33,:,0])
        #exit(1)

    def removeSpreadF(self,dataOut):
        #for i in range(11):
            #print("BEFORE Chb: ",i,dataOut.kabxys_integrated[6][:,i,0])
        #exit(1)

        #for i in range(12,data.shape[0]):
            #for j in range(data.shape[1]):
            #if data[i]>maxdB:
                #dataOut.kabxys_integrated[4][i-2:i+3,:,0] = numpy.nan #Debido a que estos ecos son intensos, se
                #dataOut.kabxys_integrated[6][i-2:i+3,:,0] = numpy.nan #remueven además dos muestras antes y después
                #dataOut.flagSpreadF = True
                #print("Removing Threshold",i)
                #print("i: ",i)

        #print("BEFORE Chb: ",dataOut.kabxys_integrated[6][:,0,0])
        #exit(1)

        #Removing outliers from the profile
        nlag = 9
        minHei = 180
        #maxHei = 600
        maxHei = 525
        inda = numpy.where(dataOut.heightList >= minHei)
        indb = numpy.where(dataOut.heightList <= maxHei)
        minIndex = inda[0][0]
        maxIndex = indb[0][-1]
        #l0 = 0
        #print("BEFORE Cha: ",dataOut.kabxys_integrated[4][:,l0,0])
        #print("BEFORE Chb: ",dataOut.kabxys_integrated[6][:,l0,0])
        #exit(1)
        '''
        l0 = 0
        #print("BEFORE Cha: ",dataOut.kabxys_integrated[4][:,l0,0])
        #print("BEFORE Chb: ",dataOut.kabxys_integrated[6][:,l0,0])

        import matplotlib.pyplot as plt
        for i in range(l0,l0+11):
            plt.plot(dataOut.kabxys_integrated[6][:,i,0],dataOut.heightList,label='{}'.format(i))
        #plt.xlim(1.e5,1.e8)
        plt.legend()
        plt.xlim(0,2000)
        plt.show()
        '''
        #dataOut.kabxys_integrated[4][minIndex:,:,0] = self.remove_coh(dataOut.kabxys_integrated[4][minIndex:,:,0
        outliers_IDs = []
        '''
        for lag in range(11):
            outliers = self.mad_based_outlier(dataOut.kabxys_integrated[4][minIndex:,lag,0], thresh=3.)
            #print("Outliers: ",outliers)
            #indexes.append(outliers.nonzero())
            #numpy.concatenate((outliers))
            #dataOut.kabxys_integrated[4][minIndex:,lag,0][outliers == True] = numpy.nan
            outliers_IDs=numpy.append(outliers_IDs,outliers.nonzero())
            '''
        '''
        for lag in range(11):
            #outliers = self.mad_based_outlier(dataOut.kabxys_integrated[6][minIndex:maxIndex,lag,0], thresh=2.)
            outliers = self.mad_based_outlier(dataOut.kabxys_integrated[6][minIndex:maxIndex,lag,0])
            outliers_IDs=numpy.append(outliers_IDs,outliers.nonzero())
            '''

        for i in range(15):
          minIndex = 12+i#12
          #maxIndex = 22+i#35
          if gmtime(dataOut.utctime).tm_hour >= 23. or gmtime(dataOut.utctime).tm_hour < 3.:
            maxIndex = 31+i#35
          else:
            maxIndex = 22+i#35
          for lag in range(11):
          #outliers = mad_based_outlier(pow_clean3[12:27], thresh=2.)
          #print("Cuts: ",first_cut*15, last_cut*15)
            outliers = self.mad_based_outlier(dataOut.kabxys_integrated[6][minIndex:maxIndex,lag,0])
            aux = minIndex+numpy.array(outliers.nonzero()).ravel()
            outliers_IDs=numpy.append(outliers_IDs,aux)
          #print(minIndex+numpy.array(outliers.nonzero()).ravel())
        #print(outliers_IDs)
        #exit(1)
        if outliers_IDs != []:
            outliers_IDs=numpy.array(outliers_IDs)
            #outliers_IDs=outliers_IDs.ravel()
            outliers_IDs=outliers_IDs.astype(numpy.dtype('int64'))
            #print(outliers_IDs)
            #exit(1)

            (uniq, freq) = (numpy.unique(outliers_IDs, return_counts=True))
            aux_arr = numpy.column_stack((uniq,freq))
            #print("repetitions: ",aux_arr)
            #exit(1)

        #if aux_arr != []:
            final_index = []
            for i in range(aux_arr.shape[0]):
                if aux_arr[i,1] >= 3*11:
                    final_index.append(aux_arr[i,0])

            if final_index != []:# and len(final_index) > 1:
                #final_index += minIndex
                #print("final_index: ",final_index)
                following_index = final_index[-1]+1 #Remove following index to ensure we remove remaining SpreadF
                previous_index = final_index[0]-1 #Remove previous index to ensure we remove remaning SpreadF
                final_index = numpy.concatenate(([previous_index],final_index,[following_index]))
                final_index = numpy.unique(final_index) #If there was only one outlier
                #print("final_index: ",final_index)
                #exit(1)
                dataOut.kabxys_integrated[4][final_index,:,0] = numpy.nan
                dataOut.kabxys_integrated[6][final_index,:,0] = numpy.nan

                dataOut.flagSpreadF = True

        #Removing echoes greater than 35 dB
        if hasattr(dataOut.pbn, "__len__"):
            maxdB = 10*numpy.log10(dataOut.pbn[0]) + 10 #Lag 0 Noise
        else:
            maxdB = 10*numpy.log10(dataOut.pbn) + 10

        #print(dataOut.kabxys_integrated[6][:,0,0])
        data = numpy.copy(10*numpy.log10(dataOut.kabxys_integrated[6][:,0,0])) #Lag0 ChB
        #print("data: ",data)

        for i in range(12,data.shape[0]):
            #for j in range(data.shape[1]):
            if data[i]>maxdB:
                dataOut.kabxys_integrated[4][i-2:i+3,:,0] = numpy.nan #Debido a que estos ecos son intensos, se
                dataOut.kabxys_integrated[6][i-2:i+3,:,0] = numpy.nan #remueven además dos muestras antes y después
                dataOut.flagSpreadF = True
                #print("Removing Threshold",i)

        #print(final_index+minIndex)
        #print(outliers_IDs)
        #exit(1)
        #print("flagSpreadF",dataOut.flagSpreadF)

        '''
        for lag in range(11):
            #print("Lag: ",lag)
            outliers = self.mad_based_outlier(dataOut.kabxys_integrated[6][minIndex:maxIndex,lag,0], thresh=2.)
            dataOut.kabxys_integrated[6][minIndex:maxIndex,lag,0][outliers == True] = numpy.nan
            '''
        #dataOut.kabxys_integrated[4][minIndex:,:,0] = self.remove_coh(dataOut.kabxys_integrated[4][minIndex:,:,0])
        '''
        import matplotlib.pyplot as plt
        for i in range(11):
            plt.plot(dataOut.kabxys_integrated[6][:,i,0],dataOut.heightList,label='{}'.format(i))
        plt.xlim(0,2000)
        plt.legend()
        plt.grid()
        plt.show()
        '''
        '''
        for nlag in range(11):
            print("BEFORE",dataOut.kabxys_integrated[6][:,nlag,0])
        #exit(1)
        '''
        #dataOut.kabxys_integrated[6][minIndex:,:,0] = self.remove_coh(dataOut.kabxys_integrated[6][minIndex:,:,0])


        '''
        for nlag in range(11):
            print("AFTER",dataOut.kabxys_integrated[6][:,nlag,0])
        exit(1)
        '''

    def run(self,dataOut):
        dataOut.flagSpreadF = False
        #print(gmtime(dataOut.utctime).tm_hour)
        #print(dataOut.ut_Faraday)
        #exit(1)
        if gmtime(dataOut.utctime).tm_hour >= 23. or gmtime(dataOut.utctime).tm_hour < 11.: #18-06 LT
            #print("Inside if we are in SpreadF Time: ",gmtime(dataOut.utctime).tm_hour)
            #if gmtime(dataOut.utctime).tm_hour == 2 and gmtime(dataOut.utctime).tm_min == 10: #Year: 2023, DOY:310
            #if gmtime(dataOut.utctime).tm_hour == 3 and gmtime(dataOut.utctime).tm_min == 10: #Year: 2023, DOY:312
            #if gmtime(dataOut.utctime).tm_hour == 0 and gmtime(dataOut.utctime).tm_min == 0: #Year: 2024, DOY:082
            #if 1: #Year: 2024, DOY:081
            #if gmtime(dataOut.utctime).tm_hour == 0: #Year: 2024, DOY:080
                #pass
            #else:
            self.removeSpreadF(dataOut)
        #exit(1)

        return dataOut

class CleanCohEchoesHP(Operation):
    '''
    Written by R. Flores
    '''
    """Operation to clean coherent echoes.

    Parameters:
    -----------
    None

    Example
    --------

    op = proc_unit.addOperation(name='CleanCohEchoes')

    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

    def remove_coh(self,pow):
        #print("pow inside: ",pow)
        #print(pow.shape)
        q75,q25 = numpy.percentile(pow,[75,25],axis=0)
        #print(q75,q25)
        intr_qr = q75-q25

        max = q75+(1.5*intr_qr)
        min = q25-(1.5*intr_qr)

        pow[pow > max] = numpy.nan

        #print("Max: ",max)
        #print("Min: ",min)

        return pow

    def mad_based_outlier_V0(self, points, thresh=3.5):
        #print("points: ",points)
        if len(points.shape) == 1:
            points = points[:,None]
        median = numpy.nanmedian(points, axis=0)
        diff = numpy.nansum((points - median)**2, axis=-1)
        diff = numpy.sqrt(diff)
        med_abs_deviation = numpy.nanmedian(diff)

        modified_z_score = 0.6745 * diff / med_abs_deviation
        #print(modified_z_score)
        return modified_z_score > thresh

    def mad_based_outlier(self, points, thresh=3.5):

        median = numpy.nanmedian(points)
        diff = (points - median)**2
        diff = numpy.sqrt(diff)
        med_abs_deviation = numpy.nanmedian(diff)

        modified_z_score = 0.6745 * diff / med_abs_deviation

        return modified_z_score > thresh

    def removeSpreadF_V0(self,dataOut):
        for i in range(11):
            print("BEFORE Chb: ",i,dataOut.kabxys_integrated[6][:,i,0])
        #exit(1)

        #Removing echoes greater than 35 dB
        maxdB = 35 #DEBERÍA SER NOISE+ALGO!!!!!!!!!!!!!!!!!!!!!!
        #print(dataOut.kabxys_integrated[6][:,0,0])
        data = numpy.copy(10*numpy.log10(dataOut.kabxys_integrated[6][:,0,0])) #Lag0 ChB
        #print(data)
        for i in range(12,data.shape[0]):
            #for j in range(data.shape[1]):
            if data[i]>maxdB:
                dataOut.kabxys_integrated[4][i-2:i+3,:,0] = numpy.nan #Debido a que estos ecos son intensos, se
                dataOut.kabxys_integrated[6][i-2:i+3,:,0] = numpy.nan #remueve además dos muestras antes y después
                #dataOut.kabxys_integrated[4][i-1,:,0] = numpy.nan
                #dataOut.kabxys_integrated[6][i-1,:,0] = numpy.nan
                #dataOut.kabxys_integrated[4][i+1,:,0] = numpy.nan
                #dataOut.kabxys_integrated[6][i+1,:,0] = numpy.nan
                dataOut.flagSpreadF = True
                print("Removing Threshold",i)
                #print("i: ",i)

        #print("BEFORE Chb: ",dataOut.kabxys_integrated[6][:,0,0])
        #exit(1)

        #Removing outliers from the profile
        nlag = 9
        minHei = 180
        #maxHei = 600
        maxHei = 525
        inda = numpy.where(dataOut.heightList >= minHei)
        indb = numpy.where(dataOut.heightList <= maxHei)
        minIndex = inda[0][0]
        maxIndex = indb[0][-1]
        #l0 = 0
        #print("BEFORE Cha: ",dataOut.kabxys_integrated[4][:,l0,0])
        #print("BEFORE Chb: ",dataOut.kabxys_integrated[6][:,l0,0])
        #exit(1)
        #'''
        l0 = 0
        #print("BEFORE Cha: ",dataOut.kabxys_integrated[4][:,l0,0])
        #print("BEFORE Chb: ",dataOut.kabxys_integrated[6][:,l0,0])

        import matplotlib.pyplot as plt
        for i in range(l0,l0+11):
            plt.plot(dataOut.kabxys_integrated[6][:,i,0],dataOut.heightList,label='{}'.format(i))
        #plt.xlim(1.e5,1.e8)
        plt.legend()
        plt.xlim(0,2000)
        plt.show()
        #'''
        #dataOut.kabxys_integrated[4][minIndex:,:,0] = self.remove_coh(dataOut.kabxys_integrated[4][minIndex:,:,0
        outliers_IDs = []
        '''
        for lag in range(11):
            outliers = self.mad_based_outlier(dataOut.kabxys_integrated[4][minIndex:,lag,0], thresh=3.)
            #print("Outliers: ",outliers)
            #indexes.append(outliers.nonzero())
            #numpy.concatenate((outliers))
            #dataOut.kabxys_integrated[4][minIndex:,lag,0][outliers == True] = numpy.nan
            outliers_IDs=numpy.append(outliers_IDs,outliers.nonzero())
            '''
        for lag in range(11):
            #outliers = self.mad_based_outlier(dataOut.kabxys_integrated[6][minIndex:maxIndex,lag,0], thresh=2.)
            outliers = self.mad_based_outlier(dataOut.kabxys_integrated[6][minIndex:maxIndex,lag,0])
            outliers_IDs=numpy.append(outliers_IDs,outliers.nonzero())
        #print(outliers_IDs)
        #exit(1)
        if outliers_IDs != []:
            outliers_IDs=numpy.array(outliers_IDs)
            outliers_IDs=outliers_IDs.ravel()
            outliers_IDs=outliers_IDs.astype(numpy.dtype('int64'))

            (uniq, freq) = (numpy.unique(outliers_IDs, return_counts=True))
            aux_arr = numpy.column_stack((uniq,freq))
            #print("repetitions: ",aux_arr)

        #if aux_arr != []:
            final_index = []
            for i in range(aux_arr.shape[0]):
                if aux_arr[i,1] >= 10:
                    final_index.append(aux_arr[i,0])

            if final_index != [] and len(final_index) > 1:
                final_index += minIndex
                #print("final_index: ",final_index)
                following_index = final_index[-1]+1 #Remove following index to ensure we remove remaining SpreadF
                previous_index = final_index[0]-1 #Remove previous index to ensure we remove remaning SpreadF
                final_index = numpy.concatenate(([previous_index],final_index,[following_index]))
                final_index = numpy.unique(final_index) #If there was only one outlier
                #print("final_index: ",final_index)
                #exit(1)
                dataOut.kabxys_integrated[4][final_index,:,0] = numpy.nan
                dataOut.kabxys_integrated[6][final_index,:,0] = numpy.nan

                dataOut.flagSpreadF = True

        #print(final_index+minIndex)
        #print(outliers_IDs)
        #exit(1)
        #print("flagSpreadF",dataOut.flagSpreadF)

        '''
        for lag in range(11):
            #print("Lag: ",lag)
            outliers = self.mad_based_outlier(dataOut.kabxys_integrated[6][minIndex:maxIndex,lag,0], thresh=2.)
            dataOut.kabxys_integrated[6][minIndex:maxIndex,lag,0][outliers == True] = numpy.nan
            '''
        #dataOut.kabxys_integrated[4][minIndex:,:,0] = self.remove_coh(dataOut.kabxys_integrated[4][minIndex:,:,0])
        '''
        import matplotlib.pyplot as plt
        for i in range(11):
            plt.plot(dataOut.kabxys_integrated[6][:,i,0],dataOut.heightList,label='{}'.format(i))
        plt.xlim(0,2000)
        plt.legend()
        plt.grid()
        plt.show()
        '''
        '''
        for nlag in range(11):
            print("BEFORE",dataOut.kabxys_integrated[6][:,nlag,0])
        #exit(1)
        '''
        #dataOut.kabxys_integrated[6][minIndex:,:,0] = self.remove_coh(dataOut.kabxys_integrated[6][minIndex:,:,0])


        '''
        for nlag in range(11):
            print("AFTER",dataOut.kabxys_integrated[6][:,nlag,0])
        exit(1)
        '''
        #print("AFTER",dataOut.kabxys_integrated[4][33,:,0])
        #print("AFTER",dataOut.kabxys_integrated[6][33,:,0])
        #exit(1)

    def removeSpreadF(self,dataOut):
        #for i in range(11):
            #print("BEFORE Chb: ",i,dataOut.kabxys_integrated[6][:,i,0])
        #exit(1)

        #for i in range(12,data.shape[0]):
            #for j in range(data.shape[1]):
            #if data[i]>maxdB:
                #dataOut.kabxys_integrated[4][i-2:i+3,:,0] = numpy.nan #Debido a que estos ecos son intensos, se
                #dataOut.kabxys_integrated[6][i-2:i+3,:,0] = numpy.nan #remueven además dos muestras antes y después
                #dataOut.flagSpreadF = True
                #print("Removing Threshold",i)
                #print("i: ",i)

        #print("BEFORE Chb: ",dataOut.kabxys_integrated[6][:,0,0])
        #exit(1)

        #Removing outliers from the profile
        nlag = 9
        minHei = 180
        #maxHei = 600
        maxHei = 525
        inda = numpy.where(dataOut.heightList >= minHei)
        indb = numpy.where(dataOut.heightList <= maxHei)
        minIndex = inda[0][0]
        maxIndex = indb[0][-1]
        #l0 = 0
        #print("BEFORE Cha: ",dataOut.kabxys_integrated[4][:,l0,0])
        #print("BEFORE Chb: ",dataOut.kabxys_integrated[6][:,l0,0])
        #exit(1)
        '''
        l0 = 0
        #print("BEFORE Cha: ",dataOut.kabxys_integrated[4][:,l0,0])
        #print("BEFORE Chb: ",dataOut.kabxys_integrated[6][:,l0,0])

        import matplotlib.pyplot as plt
        for i in range(l0,l0+11):
            plt.plot(dataOut.kabxys_integrated[6][:,i,0],dataOut.heightList,label='{}'.format(i))
        #plt.xlim(1.e5,1.e8)
        plt.legend()
        plt.xlim(0,2000)
        plt.show()
        '''
        #dataOut.kabxys_integrated[4][minIndex:,:,0] = self.remove_coh(dataOut.kabxys_integrated[4][minIndex:,:,0
        outliers_IDs = []
        '''
        for lag in range(11):
            outliers = self.mad_based_outlier(dataOut.kabxys_integrated[4][minIndex:,lag,0], thresh=3.)
            #print("Outliers: ",outliers)
            #indexes.append(outliers.nonzero())
            #numpy.concatenate((outliers))
            #dataOut.kabxys_integrated[4][minIndex:,lag,0][outliers == True] = numpy.nan
            outliers_IDs=numpy.append(outliers_IDs,outliers.nonzero())
            '''
        '''
        for lag in range(11):
            #outliers = self.mad_based_outlier(dataOut.kabxys_integrated[6][minIndex:maxIndex,lag,0], thresh=2.)
            outliers = self.mad_based_outlier(dataOut.kabxys_integrated[6][minIndex:maxIndex,lag,0])
            outliers_IDs=numpy.append(outliers_IDs,outliers.nonzero())
            '''

        for i in range(15):
          minIndex = 12+i#12
          #maxIndex = 22+i#35
          if gmtime(dataOut.utctime).tm_hour >= 23. or gmtime(dataOut.utctime).tm_hour < 3.:
            maxIndex = 31+i#35
          else:
            maxIndex = 22+i#35
          for lag in range(11):
          #outliers = mad_based_outlier(pow_clean3[12:27], thresh=2.)
          #print("Cuts: ",first_cut*15, last_cut*15)
            outliers = self.mad_based_outlier(dataOut.kabxys_integrated[6][minIndex:maxIndex,lag,0])
            aux = minIndex+numpy.array(outliers.nonzero()).ravel()
            outliers_IDs=numpy.append(outliers_IDs,aux)
          #print(minIndex+numpy.array(outliers.nonzero()).ravel())
        #print(outliers_IDs)
        #exit(1)
        if outliers_IDs != []:
            outliers_IDs=numpy.array(outliers_IDs)
            #outliers_IDs=outliers_IDs.ravel()
            outliers_IDs=outliers_IDs.astype(numpy.dtype('int64'))
            #print(outliers_IDs)
            #exit(1)

            (uniq, freq) = (numpy.unique(outliers_IDs, return_counts=True))
            aux_arr = numpy.column_stack((uniq,freq))
            #print("repetitions: ",aux_arr)
            #exit(1)

        #if aux_arr != []:
            final_index = []
            for i in range(aux_arr.shape[0]):
                if aux_arr[i,1] >= 3*11:
                    final_index.append(aux_arr[i,0])

            if final_index != []:# and len(final_index) > 1:
                #final_index += minIndex
                #print("final_index: ",final_index)
                following_index = final_index[-1]+1 #Remove following index to ensure we remove remaining SpreadF
                previous_index = final_index[0]-1 #Remove previous index to ensure we remove remaning SpreadF
                final_index = numpy.concatenate(([previous_index],final_index,[following_index]))
                final_index = numpy.unique(final_index) #If there was only one outlier
                #print("final_index: ",final_index)
                #exit(1)
                dataOut.kabxys_integrated[4][final_index,:,0] = numpy.nan
                dataOut.kabxys_integrated[6][final_index,:,0] = numpy.nan

                dataOut.flagSpreadF = True

        #Removing echoes greater than 35 dB
        maxdB = 10*numpy.log10(dataOut.pbn[0]) + 10 #Lag 0 NOise
        #maxdB = 35 #DEBERÍA SER NOISE+ALGO!!!!!!!!!!!!!!!!!!!!!!
        #print("noise: ",maxdB - 10)
        #print(dataOut.kabxys_integrated[6][:,0,0])
        data = numpy.copy(10*numpy.log10(dataOut.kabxys_integrated[6][:,0,0])) #Lag0 ChB
        #print("data: ",data)

        for i in range(12,data.shape[0]):
            #for j in range(data.shape[1]):
            if data[i]>maxdB:
                dataOut.kabxys_integrated[4][i-2:i+3,:,0] = numpy.nan #Debido a que estos ecos son intensos, se
                dataOut.kabxys_integrated[6][i-2:i+3,:,0] = numpy.nan #remueven además dos muestras antes y después
                dataOut.flagSpreadF = True
                #print("Removing Threshold",i)

        #print(final_index+minIndex)
        #print(outliers_IDs)
        #exit(1)
        #print("flagSpreadF",dataOut.flagSpreadF)

        '''
        for lag in range(11):
            #print("Lag: ",lag)
            outliers = self.mad_based_outlier(dataOut.kabxys_integrated[6][minIndex:maxIndex,lag,0], thresh=2.)
            dataOut.kabxys_integrated[6][minIndex:maxIndex,lag,0][outliers == True] = numpy.nan
            '''
        #dataOut.kabxys_integrated[4][minIndex:,:,0] = self.remove_coh(dataOut.kabxys_integrated[4][minIndex:,:,0])
        '''
        import matplotlib.pyplot as plt
        for i in range(11):
            plt.plot(dataOut.kabxys_integrated[6][:,i,0],dataOut.heightList,label='{}'.format(i))
        plt.xlim(0,2000)
        plt.legend()
        plt.grid()
        plt.show()
        '''
        '''
        for nlag in range(11):
            print("BEFORE",dataOut.kabxys_integrated[6][:,nlag,0])
        #exit(1)
        '''
        #dataOut.kabxys_integrated[6][minIndex:,:,0] = self.remove_coh(dataOut.kabxys_integrated[6][minIndex:,:,0])


        '''
        for nlag in range(11):
            print("AFTER",dataOut.kabxys_integrated[6][:,nlag,0])
        exit(1)
        '''

    def run(self,dataOut):
        dataOut.flagSpreadF = False
        #print(gmtime(dataOut.utctime).tm_hour)
        #print(dataOut.ut_Faraday)
        #exit(1)
        if gmtime(dataOut.utctime).tm_hour >= 23. or gmtime(dataOut.utctime).tm_hour < 11.: #18-06 LT
            #print("Inside if we are in SpreadF Time: ",gmtime(dataOut.utctime).tm_hour)
            self.removeSpreadF(dataOut)
        #exit(1)

        return dataOut


class NoisePower(Operation):
    '''
    Written by R. Flores
    '''
    """Operation to get noise power from the integrated data of the Double Pulse.

    Parameters:
    -----------
    None

    Example
    --------

    op = proc_unit.addOperation(name='NoisePower', optype='other')

    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

    def hildebrand(self,dataOut,data):

        divider=10 # divider was originally 10
        noise=0.0
        n1=0
        n2=int(dataOut.NDP/2)
        sorts= sorted(data)
        nums_min= dataOut.NDP/divider
        if((dataOut.NDP/divider)> 2):
            nums_min= int(dataOut.NDP/divider)

        else:
            nums_min=2
        sump=0.0
        sumq=0.0
        j=0
        cont=1
        while( (cont==1) and (j<n2)):
            sump+=sorts[j+n1]
            sumq+= sorts[j+n1]*sorts[j+n1]
            t3= sump/(j+1)
            j=j+1
            if(j> nums_min):
                rtest= float(j/(j-1)) +1.0/dataOut.NAVG
                t1= (sumq*j)
                t2=(rtest*sump*sump)
                if( (t1/t2) > 0.990):
                    j=j-1
                    sump-= sorts[j+n1]
                    sumq-=sorts[j+n1]*sorts[j+n1]
                    cont= 0

        noise= sump/j
        stdv=numpy.sqrt((sumq- noise*noise)/(j-1))
        return noise

    def run(self,dataOut):

        p=numpy.zeros((dataOut.NR,dataOut.NDP,dataOut.DPL),'float32')
        av=numpy.zeros(dataOut.NDP,'float32')
        dataOut.pnoise=numpy.zeros(dataOut.NR,'float32')

        p[0,:,:]=dataOut.kabxys_integrated[4][:,:,0]+dataOut.kabxys_integrated[5][:,:,0] #total power for channel 0, just  pulse with non-flip
        p[1,:,:]=dataOut.kabxys_integrated[6][:,:,0]+dataOut.kabxys_integrated[7][:,:,0] #total power for channel 1

        for i in range(dataOut.NR):
            dataOut.pnoise[i]=0.0
            for k in range(dataOut.DPL):
                dataOut.pnoise[i]+= self.hildebrand(dataOut,p[i,:,k])

            dataOut.pnoise[i]=dataOut.pnoise[i]/dataOut.DPL


        dataOut.pan=.8*dataOut.pnoise[0] # weights could change
        dataOut.pbn=.8*dataOut.pnoise[1] # weights could change
        '''
        print("pan: ",dataOut.pan)
        print("pbn: ",dataOut.pbn)
        print("pan dB: ",10*numpy.log10(dataOut.pan))
        print("pbn dB: ",10*numpy.log10(dataOut.pbn))
        exit(1)
        '''
        dataOut.power = dataOut.getPower()
        return dataOut


class DoublePulseACFs(Operation):
    '''
    Written by R. Flores
    '''
    """Operation to get the ACFs of the Double Pulse.

    Parameters:
    -----------
    None

    Example
    --------

    op = proc_unit.addOperation(name='DoublePulseACFs', optype='other')

    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)
        self.aux=1

    def run(self,dataOut):

        dataOut.igcej=numpy.zeros((dataOut.NDP,dataOut.DPL),'int32')
        #print("init")
        if self.aux==1:
            dataOut.rhor=numpy.zeros((dataOut.NDP,dataOut.DPL), dtype=float)
            dataOut.rhoi=numpy.zeros((dataOut.NDP,dataOut.DPL), dtype=float)
            dataOut.sdp=numpy.zeros((dataOut.NDP,dataOut.DPL), dtype=float)
            dataOut.sd=numpy.zeros((dataOut.NDP,dataOut.DPL), dtype=float)
            dataOut.p=numpy.zeros((dataOut.NDP,dataOut.DPL), dtype=float)
            dataOut.alag=numpy.zeros(dataOut.NDP,'float32')
            for l in range(dataOut.DPL):
                dataOut.alag[l]=l*dataOut.DH*2.0/150.0
            self.aux=0
        sn4=dataOut.pan*dataOut.pbn
        rhorn=0
        rhoin=0
        panrm=numpy.zeros((dataOut.NDP,dataOut.DPL), dtype=float)

        id = numpy.where(dataOut.heightList>700)[0]
        #print("kabxys: ", numpy.shape(dataOut.kabxys_integrated))
        for i in range(dataOut.NDP):
            for j in range(dataOut.DPL):
                #################  Total power
                pa=numpy.abs(dataOut.kabxys_integrated[4][i,j,0]+dataOut.kabxys_integrated[5][i,j,0])
                #print("pa::",pa)
                pb=numpy.abs(dataOut.kabxys_integrated[6][i,j,0]+dataOut.kabxys_integrated[7][i,j,0])
                st4=pa*pb

                '''
                if i > id[0]:
                    dataOut.p[i,j] pa-dataOut.pan
                else:
                    dataOut.p[i,j]=pa+pb-(dataOut.pan+dataOut.pbn)
                    '''
                #print("init 2.6",pa,dataOut.pan)
                #dataOut.pan = 23600/2
                #dataOut.pbn = 23600/2
                dataOut.p[i,j]=pa+pb-(dataOut.pan+dataOut.pbn)
                #print(i,j,dataOut.p[i,j])
                dataOut.sdp[i,j]=2*dataOut.rnint2[j]*((pa+pb)*(pa+pb))
                ## ACF

                rhorp=dataOut.kabxys_integrated[8][i,j,0]+dataOut.kabxys_integrated[11][i,j,0]
                rhoip=dataOut.kabxys_integrated[10][i,j,0]-dataOut.kabxys_integrated[9][i,j,0]
                '''
                import matplotlib.pyplot as plt
                plt.plot(numpy.abs(dataOut.kabxys_integrated[4][:,j,0]+dataOut.kabxys_integrated[5][:,j,0])+numpy.abs(dataOut.kabxys_integrated[6][:,j,0]+dataOut.kabxys_integrated[7][:,j,0]),dataOut.heightList)
                plt.axvline((dataOut.pan+dataOut.pbn))
                plt.xlim(20000,30000)
                #plt.plot(numpy.abs(dataOut.kabxys_integrated[4][:,j,0]+dataOut.kabxys_integrated[5][:,j,0])+numpy.abs(dataOut.kabxys_integrated[6][:,j,0]+dataOut.kabxys_integrated[7][:,j,0])-(dataOut.pan+dataOut.pbn),dataOut.heightList)
                #plt.xlim(1,10000)
                plt.grid()
                plt.show()
                '''
                if ((pa>dataOut.pan)&(pb>dataOut.pbn)):

                    ss4=numpy.abs((pa-dataOut.pan)*(pb-dataOut.pbn))
                    panrm[i,j]=math.sqrt(ss4)
                    rnorm=1/panrm[i,j]
                    ##  ACF
                    dataOut.rhor[i,j]=rhorp*rnorm
                    dataOut.rhoi[i,j]=rhoip*rnorm
                    #if i==13 and j== 0 or i ==14 and j==0:
                        #print(numpy.sum(rhorp))
                    #exit(1)
                    #############  Compute standard error for ACF
                    stoss4=st4/ss4
                    snoss4=sn4/ss4
                    rp2=((rhorp*rhorp)+(rhoip*rhoip))/st4
                    rn2=((rhorn*rhorn)+(rhoin*rhoin))/sn4
                    rs2=(dataOut.rhor[i,j]*dataOut.rhor[i,j])+(dataOut.rhoi[i,j]*dataOut.rhoi[i,j])
                    st=1.0+rs2*(stoss4-(2*math.sqrt(stoss4*snoss4)))
                    stn=1.0+rs2*(snoss4-(2*math.sqrt(stoss4*snoss4)))
                    dataOut.sd[i,j]=((stoss4*((1.0+rp2)*st+(2.0*rp2*rs2*snoss4)-4.0*math.sqrt(rs2*rp2)))+(0.25*snoss4*((1.0+rn2)*stn+(2.0*rn2*rs2*stoss4)-4.0*math.sqrt(rs2*rn2))))*dataOut.rnint2[j]
                    dataOut.sd[i,j]=numpy.abs(dataOut.sd[i,j])
                    '''
                    if i == 4:
                        print(i,j,dataOut.sd[i,j])
                        '''
                    #print(i,j,dataOut.rhor[i,j])
                    #exit(1)
                else: #default values for bad points
                    rnorm=1/math.sqrt(st4)
                    dataOut.sd[i,j]=1.e30
                    dataOut.ibad[i,j]=4
                    dataOut.rhor[i,j]=rhorp*rnorm
                    dataOut.rhoi[i,j]=rhoip*rnorm
                '''
                if i==47:
                    print("j",j)
                    print("pa",pa/dataOut.pan)
                    print("pb",pb/dataOut.pbn)
                    print((pa/dataOut.pan-1.0))
                    print((pb/dataOut.pbn-1.0))
                    '''
                #'''
                if ((pb/dataOut.pbn-1.0)>2.25*(pa/dataOut.pan-1.0)): #To flag bad points from the pulse and EEJ for lags != 0 for Channel B
                    #print(dataOut.heightList[i],"EJJ")
                    dataOut.igcej[i,j]=1

                elif ((pa/dataOut.pan-1.0)>2.25*(pb/dataOut.pbn-1.0)):
                    #print(dataOut.heightList[i],"EJJ")
                    dataOut.igcej[i,j]=1
                    #'''
                '''
                if ((pa/dataOut.pan-1.0)>2.25*(pb/dataOut.pbn-1.0)):
                    #print("EJJ")
                    dataOut.igcej[i,j]=1
                    '''
            '''
            if i == 4:
                exit(1)
                '''
        #print(numpy.sum(dataOut.kabxys_integrated[8][:,:,0]+dataOut.kabxys_integrated[11][:,:,0]))
        #print(numpy.sum(dataOut.kabxys_integrated[10][:,:,0]-dataOut.kabxys_integrated[9][:,:,0]))
        #print(numpy.sum(dataOut.rhor))
        #print("bad points: ",dataOut.igcej[47,:])
        #exit(1)
        '''
        for l in range(11):
            print("p: ",dataOut.p[:,l])
        exit(1)
        '''
        #print(pa)
        #print("pa: ", numpy.shape(pa))
        #print(numpy.shape(dataOut.heightList))
        '''
        import matplotlib.pyplot as plt
        plt.plot(dataOut.p[:,-1],dataOut.heightList)
        #plt.plot(pa/dataOut.pan-1.,dataOut.heightList)
        #plt.plot(pb/dataOut.pbn-1.,dataOut.heightList)
        plt.grid()
        plt.xlim(0,1e5)
        plt.show()
        #print("p: ",dataOut.p[33,:])
        #exit(1)
        #'''
        '''
        import matplotlib.pyplot as plt
        #plt.plot(numpy.abs(dataOut.kabxys_integrated[4][:,j,0]+dataOut.kabxys_integrated[5][:,j,0])+numpy.abs(dataOut.kabxys_integrated[6][:,j,0]+dataOut.kabxys_integrated[7][:,j,0]),dataOut.heightList)
        #plt.axvline((dataOut.pan+dataOut.pbn))
        print(numpy.shape(dataOut.p))
        plt.plot(dataOut.p[:,0]*dataOut.heightList*dataOut.heightList,dataOut.heightList)

        #plt.xlim(1,100000000)
        plt.xlim(100,100000000)
        plt.grid()
        plt.show()
        '''
        #print(numpy.sum(dataOut.rhor))
        #exit(1)
        return dataOut

class DoublePulseACFs_PerLag(Operation):
    '''
    Written by R. Flores
    '''
    """Operation to get the ACFs of the Double Pulse.

    Parameters:
    -----------
    None

    Example
    --------

    op = proc_unit.addOperation(name='DoublePulseACFs', optype='other')

    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)
        self.aux=1

    def run(self,dataOut):
        
        # flag bad data points
        dataOut.igcej=numpy.zeros((dataOut.NDP,dataOut.DPL),'int32')

        if self.aux==1:
            # Real part of ACF
            dataOut.rhor = numpy.zeros((dataOut.NDP, dataOut.DPL), dtype=float)
            # Imaginary part of ACF
            dataOut.rhoi = numpy.zeros((dataOut.NDP, dataOut.DPL), dtype=float)
            # Standard deviation of power
            dataOut.sdp = numpy.zeros((dataOut.NDP, dataOut.DPL), dtype=float)
            # Standard deviation of ACF
            dataOut.sd = numpy.zeros((dataOut.NDP, dataOut.DPL), dtype=float)
            # Stores the power with noise level removed
            dataOut.p=numpy.zeros((dataOut.NDP,dataOut.DPL), dtype=float)
            # Stores lags for which ACFs are calculated
            dataOut.alag=numpy.zeros(dataOut.NDP,'float32')
            for l in range(dataOut.DPL):
                dataOut.alag[l]=l*dataOut.DH*dataOut.TxLagRate/150.0
            self.aux=0
        # dataOut.pan.- Power noise level of channel A - definned in SpectraDataToFaraday
        # Signal noise
        sn4=dataOut.pan*dataOut.pbn
        rhorn=0
        rhoin=0
        panrm=numpy.zeros((dataOut.NDP,dataOut.DPL), dtype=float)

        id = numpy.where(dataOut.heightList>700)[0]

        PA = numpy.zeros((dataOut.NDP, dataOut.DPL), dtype=float)
        PB = numpy.zeros((dataOut.NDP, dataOut.DPL), dtype=float)

        for i in range(dataOut.NDP): #Heights
            for j in range(dataOut.DPL): # Lags
                #################  Total power
                # Power channel A
                pa=numpy.abs(dataOut.kabxys_integrated[4][i,j,0]+dataOut.kabxys_integrated[5][i,j,0])
                # Power channel B
                pb=numpy.abs(dataOut.kabxys_integrated[6][i,j,0]+dataOut.kabxys_integrated[7][i,j,0])
                st4=pa*pb
                '''
                if i > id[0]:
                    dataOut.p[i,j] pa-dataOut.pan
                else:
                    dataOut.p[i,j]=pa+pb-(dataOut.pan+dataOut.pbn)
                    '''
                dataOut.p[i,j]=pa+pb-(dataOut.pan[j]+dataOut.pbn[j])
                dataOut.sdp[i,j]=2*dataOut.rnint2[j]*((pa+pb)*(pa+pb))
                ## ACF
                rhorp=dataOut.kabxys_integrated[8][i,j,0]+dataOut.kabxys_integrated[11][i,j,0]
                rhoip=dataOut.kabxys_integrated[10][i,j,0]-dataOut.kabxys_integrated[9][i,j,0]
                #PA[i,j] = pa; PB[i,j] = pb
                if ((pa>dataOut.pan[j])&(pb>dataOut.pbn[j])):
                    # panrm is RMS of power, used to normalize ACFs
                    ss4=numpy.abs((pa-dataOut.pan[j])*(pb-dataOut.pbn[j]))
                    panrm[i,j]=math.sqrt(ss4)
                    rnorm=1/panrm[i,j]
                    ##  ACF
                    dataOut.rhor[i,j]=rhorp*rnorm
                    dataOut.rhoi[i,j]=rhoip*rnorm
                    #if i==13 and j== 0 or i ==14 and j==0:
                        #print(numpy.sum(rhorp))
                    #exit(1)
                    #############  Compute standard error for ACF
                    stoss4=st4/ss4
                    snoss4=sn4[j]/ss4
                    rp2=((rhorp*rhorp)+(rhoip*rhoip))/st4
                    rn2=((rhorn*rhorn)+(rhoin*rhoin))/sn4[j]
                    rs2=(dataOut.rhor[i,j]*dataOut.rhor[i,j])+(dataOut.rhoi[i,j]*dataOut.rhoi[i,j])
                    st=1.0+rs2*(stoss4-(2*math.sqrt(stoss4*snoss4)))
                    stn=1.0+rs2*(snoss4-(2*math.sqrt(stoss4*snoss4)))
                    dataOut.sd[i,j]=((stoss4*((1.0+rp2)*st+(2.0*rp2*rs2*snoss4)-4.0*math.sqrt(rs2*rp2)))+(0.25*snoss4*((1.0+rn2)*stn+(2.0*rn2*rs2*stoss4)-4.0*math.sqrt(rs2*rn2))))*dataOut.rnint2[j]
                    dataOut.sd[i,j]=numpy.abs(dataOut.sd[i,j])
                    '''
                    if i == 4:
                        print(i,j,dataOut.sd[i,j])
                        '''
                    #print(i,j,dataOut.rhor[i,j])
                    #exit(1)
                else: #default values for bad points
                    rnorm=1/math.sqrt(st4)
                    dataOut.sd[i,j]=1.e30
                    dataOut.ibad[i,j]=4
                    dataOut.rhor[i,j]=rhorp*rnorm
                    dataOut.rhoi[i,j]=rhoip*rnorm
                '''
                if i==47:
                    print("j",j)
                    print("pa",pa/dataOut.pan)
                    print("pb",pb/dataOut.pbn)
                    print((pa/dataOut.pan-1.0))
                    print((pb/dataOut.pbn-1.0))
                    '''
                #'''
                if ((pb/dataOut.pbn[j]-1.0)>2.25*(pa/dataOut.pan[j]-1.0)): #To flag bad points from the pulse and EEJ for lags != 0 for Channel B
                    #print(dataOut.heightList[i],j,"EJJ")
                    dataOut.igcej[i,j]=1

                elif ((pa/dataOut.pan[j]-1.0)>2.25*(pb/dataOut.pbn[j]-1.0)):
                    #print(dataOut.heightList[i],j,"EJJ")
                    dataOut.igcej[i,j]=1
                    #'''
                '''
                if ((pa/dataOut.pan-1.0)>2.25*(pb/dataOut.pbn-1.0)):
                    #print("EJJ")
                    dataOut.igcej[i,j]=1
                    #'''
        '''import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, dataOut.DPL, figsize=(18, 6), sharex=True, sharey=True)

        for i in range(dataOut.DPL):
            axes[0,i].plot(PA[:, i], dataOut.heightList, label=f'PA {i+1}')
            axes[0, i].axvline(dataOut.pan[i], color='gray', linestyle='--', linewidth=1) 
            axes[0, i].set_title(f'Lag {i+1}')
            axes[0, i].set_xscale("log")  # Log scale for y-axis
            axes[0, i].set_xlim([0,1e+7])
            axes[1,i].plot(PB[:, i], dataOut.heightList, label=f'PB {i+1}')
            axes[1, i].axvline(dataOut.pbn[i], color='gray', linestyle='--', linewidth=1) 
            axes[1, i].set_xscale("log")  # Log scale for y-axis
            axes[1, i].set_xlim([0,1e+7])
            
        
        plt.tight_layout()
        plt.show()'''

        
        #print("dataOut.p",datetime.datetime.utcfromtimestamp(dataOut.utctime), dataOut.p)

        #print(numpy.sum(dataOut.kabxys_integrated[8][:,:,0]+dataOut.kabxys_integrated[11][:,:,0]))
        #print(numpy.sum(dataOut.kabxys_integrated[10][:,:,0]-dataOut.kabxys_integrated[9][:,:,0]))
        #print(numpy.sum(dataOut.rhor))
        #print("bad points: ",dataOut.igcej[47,:])
        #exit(1)
        '''
        for l in range(11):
            print("p: ",dataOut.p[:,l])
        exit(1)
        '''
        #print(pa)
        '''
        import matplotlib.pyplot as plt
        #plt.plot(dataOut.p[:,-1],dataOut.heightList)
        plt.plot(pa/dataOut.pan-1.,dataOut.heightList)
        plt.plot(pb/dataOut.pbn-1.,dataOut.heightList)
        plt.grid()
        #plt.xlim(0,1e5)
        plt.show()
        #print("p: ",dataOut.p[33,:])
        #exit(1)
        '''
        return dataOut

class FaradayAngleAndDPPower(Operation):
    '''
    Written by R. Flores
    '''
    """Operation to calculate Faraday angle and Double Pulse power.

    Parameters:
    -----------
    None

    Example
    --------

    op = proc_unit.addOperation(name='FaradayAngleAndDPPower', optype='other')

    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)
        self.aux=1

    def run(self,dataOut):

        if self.aux==1:
            dataOut.h2=numpy.zeros(dataOut.MAXNRANGENDT,'float32')
            dataOut.range1=numpy.zeros(dataOut.MAXNRANGENDT,order='F',dtype='float32')
            dataOut.sdn2=numpy.zeros(dataOut.NDP,'float32')
            dataOut.ph2=numpy.zeros(dataOut.NDP,'float32')
            dataOut.sdp2=numpy.zeros(dataOut.NDP,'float32')
            dataOut.ibd=numpy.zeros(dataOut.NDP,'float32')
            dataOut.phi=numpy.zeros(dataOut.NDP,'float32')

            self.aux=0

        for i in range(dataOut.MAXNRANGENDT):
            dataOut.range1[i]=dataOut.H0 + i*dataOut.DH
            dataOut.h2[i]=dataOut.range1[i]**2
        #print("shape ph2",numpy.shape(dataOut.ph2))
        for j in range(dataOut.NDP):
            dataOut.ph2[j]=0.
            dataOut.sdp2[j]=0.
            ri=dataOut.rhoi[j][0]/dataOut.sd[j][0]
            rr=dataOut.rhor[j][0]/dataOut.sd[j][0]
            dataOut.sdn2[j]=1./dataOut.sd[j][0]

            pt=0.# // total power
            st=0.# // total signal
            ibt=0# // bad lags
            ns=0#  // no. good lags
            #print(dataOut.heightList[j])
            for l in range(dataOut.DPL):
                 #add in other lags if outside of e-jet contamination
                if( (dataOut.igcej[j][l] == 0) and (dataOut.ibad[j][l] ==  0) ):
                    #print("l", l, dataOut.range1[l])
                    dataOut.ph2[j]+=dataOut.p[j][l]/dataOut.sdp[j][l]
                    dataOut.sdp2[j]=dataOut.sdp2[j]+1./dataOut.sdp[j][l]
                    ns+=1
                #if dataOut.igcej[j][l] != 0:
                    #print(l)
                pt+=dataOut.p[j][l]/dataOut.sdp[j][l]
                st+=1./dataOut.sdp[j][l]
                ibt|=dataOut.ibad[j][l];
            #print(dataOut.sdp2[j],st)
            if(ns!= 0):
                #print("Good lags: ",j,ns)
                dataOut.ibd[j]=0
                dataOut.ph2[j]=dataOut.ph2[j]/dataOut.sdp2[j]
                dataOut.sdp2[j]=1./dataOut.sdp2[j]
                #print("j", dataOut.range1[j])
                #print(dataOut.sdp2[j])
            else:
                #print("Bad lags: ",j)
                #print("Bad j", dataOut.range1[j])
                dataOut.ibd[j]=ibt
                dataOut.ph2[j]=pt/st
                dataOut.sdp2[j]=1./st
                #dataOut.sdp2[j]**=2 #Added on May 22, 2024 by R.Flores #To Validate
                #print(dataOut.sdp2[j])

            dataOut.ph2[j]=dataOut.ph2[j]*dataOut.h2[j]
            dataOut.sdp2[j]=numpy.sqrt(dataOut.sdp2[j])*dataOut.h2[j]
            rr=rr/dataOut.sdn2[j]
            ri=ri/dataOut.sdn2[j]
            #rm[j]=np.sqrt(rr*rr + ri*ri) it is not used in c program
            dataOut.sdn2[j]=1./(dataOut.sdn2[j]*(rr*rr + ri*ri))
            if( (ri == 0.) and (rr == 0.) ):
                dataOut.phi[j]=0.
            else:
                dataOut.phi[j]=math.atan2( ri , rr )

        dataOut.flagTeTiCorrection = False
        #print("ph2: ", numpy.sum(dataOut.ph2[:16]))
        #print("ph2: ", numpy.sum(dataOut.ph2[16:32]))
        
        '''import matplotlib.pyplot as plt
        #plt.plot(numpy.abs(dataOut.kabxys_integrated[4][:,j,0]+dataOut.kabxys_integrated[5][:,j,0])+numpy.abs(dataOut.kabxys_integrated[6][:,j,0]+dataOut.kabxys_integrated[7][:,j,0]),dataOut.heightList)
        #plt.axvline((dataOut.pan+dataOut.pbn))
        #print(numpy.shape(dataOut.p))
        plt.plot(dataOut.ph2,dataOut.heightList, 'r', label='$P_{H_2}$')
        plt.plot(dataOut.phi,dataOut.heightList, 'k--', label='$\Phi$')
        plt.legend()
        plt.xlim(1000,1000000000)
        #plt.ylim(50,400)
        plt.grid()
        plt.show()
        #exit(1)'''
        
        return dataOut

class ElectronDensityFaraday(Operation):
    '''
    Written by R. Flores
    '''
    """Operation to calculate electron density from Faraday angle.

    Parameters:
    -----------
    NSHTS : int
        .*
    RATE : float
        .*

    Example
    --------

    op = proc_unit.addOperation(name='ElectronDensityFaraday', optype='other')
    op.addParameter(name='NSHTS', value='50', format='int')
    op.addParameter(name='RATE', value='1.8978873e-6', format='float')

    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)
        self.aux=1

    def run(self,dataOut,NSHTS=50,RATE=1.8978873e-6):

        dataOut.NSHTS=NSHTS
        dataOut.RATE=RATE

        if self.aux==1:
            dataOut.dphi=numpy.zeros(dataOut.NDP,'float32')
            #dataOut.dphi_uc=numpy.zeros(dataOut.NDP,'float32')
            dataOut.sdn1=numpy.zeros(dataOut.NDP,'float32')
            self.aux=0
        theta=numpy.zeros(dataOut.NDP,dtype=numpy.complex_)
        thetai=numpy.zeros(dataOut.NDP,dtype=numpy.complex_)
        # use complex numbers for phase
        '''
        for i in range(dataOut.NSHTS):
            theta[i]=math.cos(dataOut.phi[i])+math.sin(dataOut.phi[i])*1j
            thetai[i]=-math.sin(dataOut.phi[i])+math.cos(dataOut.phi[i])*1j
            ''' #Old Method
        
        # differentiate and convert to number density
        ndphi=dataOut.NSHTS-4
        #print(dataOut.phi)
        #exit(1)
        #'''
        if hasattr(dataOut, 'flagSpreadF') and dataOut.flagSpreadF:
        #if dataOut.flagSpreadF:
            nanindex = numpy.argwhere(numpy.isnan(dataOut.phi))
            i1 = nanindex[-1][0]
            #Analizar cuando SpreadF es Pluma

        #print(i1)
            dataOut.phi[i1+1:]=numpy.unwrap(dataOut.phi[i1+1:]) #Better results
        else:
            #dataOut.phi_uwrp = dataOut.phi.copy()
            dataOut.phi[:]=numpy.unwrap(dataOut.phi[:]) #Better results
            #'''
        #print(dataOut.phi)
        #print(dataOut.ph2)
        #exit(1)
        
        '''
        #if dataOut.flagDecodeData:
        if 1:
            import matplotlib.pyplot as plt
            plt.plot(dataOut.phi,dataOut.heightList,'*-')
            #plt.ylim(60,95)
            plt.grid()
            plt.show()
            '''
        #print(dataOut.bki)
        #print(dataOut.NDP,dataOut.NSHTS)
        #print("phi: ", dataOut.phi)
        for i in range(2,dataOut.NSHTS-2):
            fact=(-0.5/(dataOut.RATE*dataOut.DH))*dataOut.bki[i]
            #print("fact: ", fact,dataOut.RATE,dataOut.DH,dataOut.bki[i])
            #four-point derivative, no phase unwrapping necessary
            #####dataOut.dphi[i]=((((theta[i+1]-theta[i-1])+(2.0*(theta[i+2]-theta[i-2])))/thetai[i])).real/10.0 #Original from C program

            ##dataOut.dphi[i]=((((theta[i-2]-theta[i+2])+(8.0*(theta[i+1]-theta[i-1])))/thetai[i])).real/12.0
            dataOut.dphi[i]=((dataOut.phi[i+1]-dataOut.phi[i-1])+(2.0*(dataOut.phi[i+2]-dataOut.phi[i-2])))/10.0 #Better results

            #dataOut.dphi_uc[i] = abs(dataOut.phi[i]*dataOut.bki[i]*(-0.5)/dataOut.DH)
            #dataOut.dphi[i]=abs(dataOut.dphi[i]*fact)
            dataOut.dphi[i]=dataOut.dphi[i]*abs(fact)
            dataOut.sdn1[i]=(4.*(dataOut.sdn2[i-2]+dataOut.sdn2[i+2])+dataOut.sdn2[i-1]+dataOut.sdn2[i+1])
            dataOut.sdn1[i]=numpy.sqrt(dataOut.sdn1[i])*fact

        #print("dphi: ", dataOut.dphi)
        
        '''if True: #dataOut.flagDecodeData:
            #exit(1)
            import matplotlib.pyplot as plt
            plt.plot(abs(dataOut.dphi),dataOut.heightList)
            plt.grid()
        #plt.xlim(0,1e7)
            plt.show()'''

            
        #print("dH: ", dataOut.heightList[1]-dataOut.heightList[0])
        
    
        

        return dataOut


class NormalizeDPPowerRoberto_V2(Operation):
    '''
    Written by R. Flores
    '''
    """Operation to normalize relative electron density from power with total electron density from Farday angle.

    Parameters:
    -----------
    None

    Example
    --------

    op = proc_unit.addOperation(name='NormalizeDPPower', optype='other')

    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)
        self.aux=1

    def normal(self,a,b,n,m):
        chmin=1.0e30
        chisq=numpy.zeros(150,'float32')
        temp=numpy.zeros(150,'float32')

        for i in range(2*m-1):
            an=al=be=chisq[i]=0.0
            for j in range(int(n/m)):
                k=int(j+i*n/(2*m))
                if(a[k]>0.0 and b[k]>0.0):
                    al+=a[k]*b[k]
                    be+=b[k]*b[k]

            if(be>0.0):
                temp[i]=al/be
            else:
                temp[i]=1.0

            for j in range(int(n/m)):
                k=int(j+i*n/(2*m))
                #print("a,b",a[k],b[k])
                if(a[k]>0.0 and b[k]>0.0):
                    chisq[i]+=(numpy.log10(b[k]*temp[i]/a[k]))**2
                    an=an+1

            if(chisq[i]>0.0):
                chisq[i]/=an

        for i in range(int(2*m-1)):
            #print("xi",chisq[i])
            if(chisq[i]<chmin and chisq[i]>1.0e-6):
                chmin=chisq[i]
                cf=temp[i]
        return cf

 

    def normalize(self,dataOut):
        # cf .- constant factor of normalization
        
        if self.aux==1:
            dataOut.cf=numpy.zeros(1,'float32')
            dataOut.cflast=numpy.zeros(1,'float32')
            self.aux=0

        print(dataOut.ut_Faraday)

        if (dataOut.ut_Faraday>=11.5 and dataOut.ut_Faraday<23): # 6 30am to 6pm
            i2=(500.-dataOut.range1[0])/dataOut.DH
            i1=(200.-dataOut.range1[0])/dataOut.DH

        elif(dataOut.ut_Faraday>=5 and dataOut.ut_Faraday<8): # 0 am to 3am
            inda = numpy.where(dataOut.heightList >= 330) # 260 #200 km 
            minIndex = inda[0][0]
            indb = numpy.where(dataOut.heightList < 470) # 350 # 700 km
            maxIndex = indb[0][-1]
            print(minIndex)
            print(dataOut.heightList)

            ph2max_idx = numpy.nanargmax(dataOut.ph2[minIndex:maxIndex])
            #print("dataOut.ph2[minIndex:maxIndex]: ", dataOut.ph2[minIndex:maxIndex])
            #print("dataOut.ph2: ", dataOut.ph2)
            #print("dataOut.phi: ", dataOut.phi)
            print("minIndex", minIndex, "maxIndex", maxIndex)
            print(ph2max_idx)

            ph2max_idx += minIndex

            i2 = maxIndex #ph2max_idx + 6
            i1 = minIndex #ph2max_idx - 6

            print("ELSE^^^^^^^^^^^^^^^^^^^^^")
            print(dataOut.heightList[i1])
            print(dataOut.heightList[i2])
        elif(dataOut.ut_Faraday>=8 and dataOut.ut_Faraday<11.5): # 3 am to 6 30am ADDED
            inda = numpy.where(dataOut.heightList >= 260) # 260 #200 km 
            minIndex = inda[0][0]
            indb = numpy.where(dataOut.heightList < 350) # 350 # 700 km
            maxIndex = indb[0][-1]

            ph2max_idx = numpy.nanargmax(dataOut.ph2[minIndex:maxIndex])
            #print("dataOut.ph2[minIndex:maxIndex]: ", dataOut.ph2[minIndex:maxIndex])
            #print("dataOut.ph2: ", dataOut.ph2)
            print("minIndex", minIndex, "maxIndex", maxIndex)
            print(ph2max_idx)

            ph2max_idx += minIndex

            i2 = maxIndex #ph2max_idx + 6
            i1 = minIndex #ph2max_idx - 6

            print("ELSE^^^^^^^^^^^^^^^^^^^^^")
            print(dataOut.heightList[i1])
            print(dataOut.heightList[i2])
        else: #6pm to 12am
            inda = numpy.where(dataOut.heightList >= 200) #200 km 
            minIndex = inda[0][0]
            indb = numpy.where(dataOut.heightList < 700) # 700 km
            maxIndex = indb[0][-1]
            print(minIndex)
            print(dataOut.heightList)

            ph2max_idx = numpy.nanargmax(dataOut.ph2[minIndex:maxIndex])
            #print("dataOut.ph2[minIndex:maxIndex]: ", dataOut.ph2[minIndex:maxIndex])
            #print("dataOut.ph2: ", dataOut.ph2)
            '''
            #if dataOut.flagTeTiCorrection:
            if 1:
                import matplotlib.pyplot as plt
                plt.figure()
                plt.plot(dataOut.ph2[minIndex:maxIndex],dataOut.heightList[minIndex:maxIndex],'*-')
                plt.show()
                '''
            ph2max_idx += minIndex

            i2 = ph2max_idx + 6
            i1 = ph2max_idx - 6

            print(dataOut.heightList[i1])
            print(dataOut.heightList[i2])
        '''
        elif(dataOut.ut_Faraday>=1.66 and dataOut.ut_Faraday<3.16): # 20 40 - 22 10 (1 40 - 3 10)
            inda = numpy.where(dataOut.heightList >= 435) #200 km 
            minIndex = inda[0][0]
            indb = numpy.where(dataOut.heightList < 480) # 700 km
            maxIndex = indb[0][-1]
            print(minIndex)
            print(dataOut.heightList)

            ph2max_idx = numpy.nanargmax(dataOut.ph2[minIndex:maxIndex])
            #print("dataOut.ph2[minIndex:maxIndex]: ", dataOut.ph2[minIndex:maxIndex])
            #print("dataOut.ph2: ", dataOut.ph2)
            print("minIndex", minIndex, "maxIndex", maxIndex)
            print(ph2max_idx)
            
            ph2max_idx += minIndex

            i2 = maxIndex #ph2max_idx + 6
            i1 = minIndex #ph2max_idx - 6

            print("ELSE^^^^^^^^^^^^^^^^^^^^^")
            print(dataOut.heightList[i1])
            print(dataOut.heightList[i2])
        '''
        

        try:
            dataOut.heightList[i2]
        except:
            i2 -= 1

        '''
        if not dataOut.flagSpreadF:
            i2=(420-dataOut.range1[0])/dataOut.DH
        else:
            i2=(620-dataOut.range1[0])/dataOut.DH
            '''
        #i1=(200 -dataOut.range1[0])/dataOut.DH
        ##print(i1*dataOut.DH)
        ##print(i2*dataOut.DH)

        i1=int(i1)
        i2=int(i2)
        print("Bounds 1: ", dataOut.heightList[i1],dataOut.heightList[i2])
        '''
        print(dataOut.ph2)
        import matplotlib.pyplot as plt
        plt.plot(dataOut.ph2,dataOut.heightList)
        plt.xlim(1.e5,1.e8)
        plt.show()
        '''
        #print("Flag: ",dataOut.flagTeTiCorrection)
        #print(dataOut.dphi[i1::])
        #print(dataOut.ph2[:])

        if dataOut.flagTeTiCorrection:
            for i in range(dataOut.NSHTS):
                dataOut.ph2[i]/=dataOut.cf
                dataOut.sdp2[i]/=dataOut.cf

        #'''
        #if dataOut.flagSpreadF:
        if hasattr(dataOut, 'flagSpreadF') and dataOut.flagSpreadF:
            print("flagSpreadF activated!!")
            i2=int((700-dataOut.range1[0])/dataOut.DH)
            #print(dataOut.ph2)
            #print(dataOut.heightList)
            nanindex = numpy.argwhere(numpy.isnan(dataOut.ph2))
            #print("nanindex",nanindex)
            i1 = nanindex[-1][0] #VER CUANDO i1>i2
            if i1 != numpy.shape(dataOut.heightList)[0]:
                i1 += 1+2 #Se suma uno para no tomar el nan, se suma 2 para no tomar datos nan de "phi" debido al calculo de la derivada
            if i1 >= i2:
                i1 = i2-4
        #print("i1, i2",i1,i2)
        #print(dataOut.heightList)
        #print("Bounds: ", dataOut.heightList[i1],dataOut.heightList[i2])
        #print(dataOut.dphi[33])
        #print(dataOut.ph2[33])
        #print(dataOut.dphi[i1::])
        #print(dataOut.ph2[i1::])
        #'''
        print("Bounds 2: ", dataOut.heightList[i1],dataOut.heightList[i2])

        print("\n========== ANTES DE NORMALIZAR ==========")
        print(f"Tiempo UTC: {datetime.datetime.utcfromtimestamp(dataOut.utctime)}")
        print("Rango de alturas impreso: 300.0 - 500.0 km")

        print(f"{'Altura(km)':>12} {'Faraday':>15} {'Power':>15}")

        # Índices correspondientes a 200 y 500 km
        idx1 = numpy.where(dataOut.heightList >= 200)[0][0]
        idx2 = numpy.where(dataOut.heightList <= 500)[0][-1]

        for i in range(idx1, idx2 + 1):
            print(f"{dataOut.heightList[i]:12.1f} "
                f"{dataOut.dphi[i]:15.6e} "
                f"{dataOut.ph2[i]:15.6e}")

        try:
            dataOut.cf=self.normal(dataOut.dphi[i1::], dataOut.ph2[i1::], i2-i1, 1)

        except:
            print("except: chi factor not achieved in normalization")
            dataOut.cf = numpy.nan

        #print("cf: ",dataOut.cf)
        #print(dataOut.ph2)
        #input()
        #  in case of spread F, normalize much higher
        #print("dens: ", dataOut.dphi,dataOut.ph2)

        night_first1= 300.0#350.0
        night_end= 450.0
        night_first1= 220.0#350.0
        night_end= 400.0

        if(dataOut.cf<dataOut.cflast[0]/10.0):
            i1=(night_first1-dataOut.range1[0])/dataOut.DH
            i2=(night_end-dataOut.range1[0])/dataOut.DH
            i1=int(i1)
            i2=int(i2)
            try:
                dataOut.cf=self.normal(dataOut.dphi[int(i1)::], dataOut.ph2[int(i1)::], int(i2-i1), 1)
            except:
                pass

        #print(dataOut.cf,dataOut.cflast[0])
        

        ### Manual cf correction ###
        # To Do: Pasarlo a funcion aparte - usar una sola rutina para DP e Hybrid - local path

        flagcfcorrection = True

        time_text = datetime.datetime.utcfromtimestamp(dataOut.utctime) # time in UTC
        DOY = time_text.timetuple().tm_yday
        print("Bounds 3: ", dataOut.heightList[i1], dataOut.heightList[i2])
        print('time text', time_text, DOY)
        if flagcfcorrection:
            print("***Cleaning*** cflast: ", dataOut.cflast[0])
            print("***Cleaning*** cf: ", dataOut.cf)
            path = os.path.join(os.path.dirname(__file__), 'normalize_factor2.json')
            with open(path) as f:
                jsondata= json.load(f)

            corrections = {}
            for condition in jsondata['conditions']:
                year = condition['year']
                doy = condition['doy']
                cf = condition['cf']
                
                for time_condition in condition['time']:
                    hour = time_condition[0]
                    minute = time_condition[1]
                    corrections[(year, doy, hour, minute)] = cf
            key = (time_text.year, DOY, time_text.hour, time_text.minute)

            if key in corrections:
                cf_value = corrections[key]
                
                if isinstance(cf_value, str) and "cflast" in cf_value:
                    dataOut.cf = eval(cf_value)  #ast.literal_eval(cf_value) 
                else:
                    dataOut.cf = float(cf_value)
                
                print(f"Correction applied: {dataOut.cf}")
            print("***Cleaning*** cf After: ", dataOut.cf)
        
        ###
        import pandas as pd
        localtime = dataOut.utctime - 5*3600
        DOY = (datetime.datetime.utcfromtimestamp(localtime)).timetuple().tm_yday # DOY in UTC-5
        YEAR = (datetime.datetime.utcfromtimestamp(localtime)).year # YEAR in UTC-5
        path = os.path.join(os.path.dirname(__file__), f'cf{str(YEAR)}{str(DOY)}.csv')

        if ~os.path.exists(path): flagcfcorrection = False; print("No cf path")

        if flagcfcorrection:
            df = pd.read_csv(path)
            cf_time = df['timestamp'].to_numpy() 
            cf_cf = df['cf'].to_numpy() 

            dt_num = time_text   # Data time
            dt_array = pd.to_datetime(cf_time, unit="s", utc=True) # File time

            mask = (
                (dt_array.year.to_numpy() == dt_num.year) &
                (dt_array.dayofyear.to_numpy() == dt_num.timetuple().tm_yday) &
                (dt_array.hour.to_numpy() == dt_num.hour) &
                (dt_array.minute.to_numpy() == dt_num.minute)
            )

            indices = numpy.where(mask)[0]

            print("CHECK" ,dataOut.utctime, time_text, indices)
            print("CHECK2", cf_time)
        # Change cf
        try:
            cf_index = indices[0]
            print("cf changed")
            dataOut.cf = cf_cf[cf_index]
            print("***Cleaning*** cf After: ", dataOut.cf)
        except:
            print("not changed")
        ###

        # Update cf last
        dataOut.cflast[0]=dataOut.cf

        #print(dataOut.ph2)
        #print(dataOut.sdp2)

        ## normalize double pulse power and error bars to Faraday
        for i in range(dataOut.NSHTS):
            dataOut.ph2[i]*=dataOut.cf
            dataOut.sdp2[i]*=dataOut.cf
        #print(dataOut.ph2)
        #print(dataOut.sdp2)
        #exit(1)

        print("\n========== DESPUÉS DE NORMALIZAR ==========")
        print(f"cf = {dataOut.cf:.6e}")
        print(f"{'Altura(km)':>12} {'Faraday':>15} {'Power norm':>15} {'Relación':>12}")

        # Índices para 200-500 km
        idx1 = numpy.where(dataOut.heightList >= 100)[0][0]
        idx2 = numpy.where(dataOut.heightList <= 800)[0][-1]

        for i in range(idx1, idx2 + 1):
            ratio = dataOut.ph2[i] / dataOut.dphi[i] if dataOut.dphi[i] != 0 else numpy.nan

            print(f"{dataOut.heightList[i]:12.1f} "
                f"{dataOut.dphi[i]:15.6e} "
                f"{dataOut.ph2[i]:15.6e} "
                f"{ratio:12.4f} "
                f"{i}")

        for i in range(dataOut.NSHTS):
            dataOut.ph2[i]=(max(1.0, dataOut.ph2[i]))
            dataOut.dphi[i]=(max(1.0, dataOut.dphi[i]))

    def run(self,dataOut):

        self.normalize(dataOut)
        #print(dataOut.ph2)
        #print(dataOut.sdp2)
        #input()
        #print("shape before" ,numpy.shape(dataOut.ph2))

        return dataOut


class suppress_stdout_stderr(object):
    '''
    A context manager for doing a "deep suppression" of stdout and stderr in
    Python, i.e. will suppress all print, even if the print originates in a
    compiled C/Fortran sub-function.
       This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).

    '''
    def __init__(self):
        # Open a pair of null files
        self.null_fds =  [os.open(os.devnull,os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = [os.dup(1), os.dup(2)]

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0],1)
        os.dup2(self.null_fds[1],2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0],1)
        os.dup2(self.save_fds[1],2)
        # Close all file descriptors
        for fd in self.null_fds + self.save_fds:
            os.close(fd)


class DPTemperaturesEstimation(Operation):
    '''
    Written by R. Flores
    '''
    """Operation to estimate temperatures for Double Pulse data.

    Parameters:
    -----------
    IBITS : int
        .*

    Example
    --------

    op = proc_unit.addOperation(name='DPTemperaturesEstimation', optype='other')
    op.addParameter(name='IBITS', value='16', format='int')

    """
    '''
    NSTHS Number of sample heights (input), for temperature processing
    NDP Number of Data Points (nHeights dependent)
    NSTHS < NDP

    Input/Output Data:

        te2, ti2: Estimated electron and ion temperatures.

        ete2, eti2: Errors in the estimated temperatures.

        phy2, ephy2: Physical parameter and its error.

    Fitting Process:

        ifit: Flags for which parameters are being fitted.

        params: Initial guesses and fitted parameters.

        cov, covinv: Covariance matrix and its inverse for uncertainty estimation.

    Metadata/Status:

        m: Status or counter for the fitting process.

        info2: Additional information about the fitting results.

    '''


    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

        self.aux=1

    def Estimation(self,dataOut):
    #with suppress_stdout_stderr():

        if self.aux==1:
            dataOut.ifit=numpy.zeros(5,order='F',dtype='int32')
            dataOut.m=numpy.zeros(1,order='F',dtype='int32')
            dataOut.te2=numpy.zeros(dataOut.NSHTS,order='F',dtype='float32')
            dataOut.ti2=numpy.zeros(dataOut.NSHTS,order='F',dtype='float32')
            dataOut.ete2=numpy.zeros(dataOut.NSHTS,order='F',dtype='float32')
            dataOut.eti2=numpy.zeros(dataOut.NSHTS,order='F',dtype='float32')

            self.aux=0

        dataOut.phy2=numpy.zeros(dataOut.NSHTS,order='F',dtype='float32')
        dataOut.ephy2=numpy.zeros(dataOut.NSHTS,order='F',dtype='float32')
        dataOut.info2=numpy.zeros(dataOut.NDP,order='F',dtype='float32')
        dataOut.params=numpy.zeros(10,order='F',dtype='float32')
        dataOut.cov=numpy.zeros(dataOut.IBITS*dataOut.IBITS,order='F',dtype='float32')
        dataOut.covinv=numpy.zeros(dataOut.IBITS*dataOut.IBITS,order='F',dtype='float32')

        #null_fd = os.open(os.devnull, os.O_RDWR)
        #os.dup2(null_fd, 1)
        ymin_index =  numpy.abs(dataOut.heightList - 150).argmin()  #no point below 150 km
        for i in range(ymin_index,dataOut.NSHTS): 

            #some definitions
            iflag=0 # inicializado a cero?
            wl = 3.0
            x=numpy.zeros(dataOut.DPL+dataOut.IBITS,order='F',dtype='float32')
            y=numpy.zeros(dataOut.DPL+dataOut.IBITS,order='F',dtype='float32')
            e=numpy.zeros(dataOut.DPL+dataOut.IBITS,order='F',dtype='float32')
            eb=numpy.zeros(5,order='F',dtype='float32')
            zero=numpy.zeros(1,order='F',dtype='float32')
            depth=numpy.zeros(1,order='F',dtype='float32')
            t1=numpy.zeros(1,order='F',dtype='float32')
            t2=numpy.zeros(1,order='F',dtype='float32')
            
            '''
            x lag time      y correlation       e their uncertanities
            t1 t2 initial guesses               eb errors
            '''
            if i>ymin_index and l1>=0:
                if l1==0:
                    l1=1

                dataOut.cov=numpy.reshape(dataOut.cov,l1*l1)
                dataOut.cov=numpy.resize(dataOut.cov,dataOut.DPL*dataOut.DPL)
                dataOut.covinv=numpy.reshape(dataOut.covinv,l1*l1)
                dataOut.covinv=numpy.resize(dataOut.covinv,dataOut.DPL*dataOut.DPL)

            for l in range(dataOut.DPL*dataOut.DPL):
                dataOut.cov[l]=0.0
            acfm= (dataOut.rhor[i][0])**2 + (dataOut.rhoi[i][0])**2
            if acfm> 0.0:
                cc=dataOut.rhor[i][0]/acfm
                ss=dataOut.rhoi[i][0]/acfm
            else:
                cc=1.
                ss=0.
            # keep only uncontaminated data, don't pass zero lag to fitter
            l1=0
            '''
            if i==13 or i==14:
                print(numpy.sum(dataOut.rhor))
                print(numpy.sum(dataOut.rhoi))
                print(acfm)
                #exit(1)
                '''
            for l in range(0+1,dataOut.DPL):
                if dataOut.igcej[i][l]==0 and dataOut.ibad[i][l]==0:
                    y[l1]=dataOut.rhor[i][l]*cc + dataOut.rhoi[i][l]*ss
                    x[l1]=dataOut.alag[l]*1.0e-3 # *1.0e-3
                    dataOut.sd[i][l]=dataOut.sd[i][l]/((acfm)**2)# important
                    e[l1]=dataOut.sd[i][l] #this is the variance, not the st. dev.
                    l1=l1+1

            for l in range(l1*(l1+1)):
                dataOut.cov[l]=0.0
            for l in range(l1):
                dataOut.cov[l*(1+l1)]=e[l]
            angle=dataOut.thb[i]*0.01745
            bm=dataOut.bfm[i]
            dataOut.params[0]=1.0     #norm
            dataOut.params[1]=1000.0  #te
            dataOut.params[2]=800.0   #ti
            dataOut.params[3]=0.00    #ph
            dataOut.params[4]=0.00    #phe

            if l1!=0:
                x=numpy.resize(x,l1)
                y=numpy.resize(y,l1)
            else:
                x=numpy.resize(x,1)
                y=numpy.resize(y,1)

            if True: #len(y)!=0:
                with suppress_stdout_stderr():
                    fitacf_guess.guess(y,x,zero,depth,t1,t2,len(y)) #t1 = te , t2 = tr = te/ti
                t2=t1/t2  # ti

                if (t1<5000.0 and t1> 600.0):
                    dataOut.params[1]=t1
                    dataOut.params[2]=min(t2,t1)
                dataOut.ifit[1]=dataOut.ifit[2]=1
                dataOut.ifit[0]=dataOut.ifit[3]=dataOut.ifit[4]=0

                if dataOut.ut_Faraday<10.0 and dataOut.ut_Faraday>=0.5: # 6 30 pm to 5 am LT
                    dataOut.ifit[2]=0

                den=dataOut.ph2[i]

                if l1!=0:
                    dataOut.covinv=dataOut.covinv[0:l1*l1].reshape((l1,l1))
                    dataOut.cov=dataOut.cov[0:l1*l1].reshape((l1,l1))
                    e=numpy.resize(e,l1)
                else:
                    dataOut.covinv=numpy.resize(dataOut.covinv,1)
                    dataOut.cov=numpy.resize(dataOut.cov,1)
                    e=numpy.resize(e,1)

                eb=numpy.resize(eb,10)
                dataOut.ifit=numpy.resize(dataOut.ifit,10)
                #print("*********************FITACF_FIT*********************",dataOut.covinv,e,dataOut.params,eb,dataOut.m)
                #exit(1)
                with suppress_stdout_stderr():
                    dataOut.covinv,e,dataOut.params,eb,dataOut.m=fitacf_fit_short.fit(wl,x,y,dataOut.cov,dataOut.covinv,e,dataOut.params,bm,angle,den,dataOut.range1[i],dataOut.year,dataOut.ifit,dataOut.m,l1) #
                #exit(1)
                if dataOut.params[2]>dataOut.params[1]*1.05:
                    dataOut.ifit[2]=0
                    dataOut.params[1]=dataOut.params[2]=t1
                    with suppress_stdout_stderr():
                        dataOut.covinv,e,dataOut.params,eb,dataOut.m=fitacf_fit_short.fit(wl,x,y,dataOut.cov,dataOut.covinv,e,dataOut.params,bm,angle,den,dataOut.range1[i],dataOut.year,dataOut.ifit,dataOut.m,l1) #
                #print("*********************FIT SUCCESS*********************",dataOut.covinv,e,dataOut.params,eb,dataOut.m)
                #exit(1)
                if (dataOut.ifit[2]==0):
                     dataOut.params[2]=dataOut.params[1]
                if (dataOut.ifit[3]==0 and iflag==0):
                    dataOut.params[3]=0.0
                if (dataOut.ifit[4]==0):
                    dataOut.params[4]=0.0
                dataOut.te2[i]=dataOut.params[1]
                #if i == 13 or i ==14 or i ==15:
                    #print(dataOut.te2[i])
                dataOut.ti2[i]=dataOut.params[2]
                dataOut.ete2[i]=eb[1]
                dataOut.eti2[i]=eb[2]

                if dataOut.eti2[i]==0:
                    dataOut.eti2[i]=dataOut.ete2[i]

                dataOut.phy2[i]=dataOut.params[3]
                dataOut.ephy2[i]=eb[3]
                if(iflag==1):
                    dataOut.ephy2[i]=0.0

                if (dataOut.m<=3 and dataOut.m!= 0 and dataOut.te2[i]>400.0):
                    dataOut.info2[i]=1
                else:
                    dataOut.info2[i]=0

    def gaussian(self, x, a, b, c, d):
        val = a * numpy.exp(-(x - b)**2 / (2*c**2)) + d
        return val

    def run(self,dataOut,IBITS=16):

        dataOut.IBITS = IBITS
        self.Estimation(dataOut)
        '''
        from scipy.optimize import least_squares

        ratio = dataOut.te2/dataOut.ti2

        ratio[:11] = ratio[:11] = 1.
        ratio[20:] = ratio[20:] = 1.

        #print(ratio)
        #print(ratio-self.gaussian(dataOut.heightList[:dataOut.NSHTS],1,250,20,1))
        def lsq_func(params):
            return (ratio-self.gaussian(dataOut.heightList[:dataOut.NSHTS],params[0],params[1],params[2],params[3]))


        x0_value = numpy.array([1,250,20,1])

        popt = least_squares(lsq_func,x0=x0_value,verbose=0)

        A = popt.x[0]; B = popt.x[1]; C = popt.x[2]; D = popt.x[3]

        aux = self.gaussian(dataOut.heightList[:dataOut.NSHTS], A, B, C, D)

        import matplotlib.pyplot as plt
        #plt.plot(dataOut.te2,dataOut.heightList[:dataOut.NSHTS])
        #plt.plot(dataOut.ti2,dataOut.heightList[:dataOut.NSHTS])
        #plt.xlim(0,5000)
        plt.plot(ratio,dataOut.heightList[:dataOut.NSHTS])
        plt.plot(aux,dataOut.heightList[:dataOut.NSHTS])
        #plt.plot(self.gaussian(dataOut.heightList[:dataOut.NSHTS],1,250,20,1),dataOut.heightList[:dataOut.NSHTS])

        plt.show()
        '''

        return dataOut



class DenCorrection(NormalizeDPPowerRoberto_V2):
    '''
    Written by R. Flores
    '''
    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)
        self.aux = 0
        self.csv_flag = 1

    def gaussian(self, x, a, b, c):
        val = a * numpy.exp(-(x - b)**2 / (2*c**2))
        return val

    def TeTiEstimation(self,dataOut):

        #dataOut.DPL = 2 #for MST
        y=numpy.zeros(dataOut.DPL,order='F',dtype='float32')

        #y_aux = numpy.zeros(1,,dtype='float32')
        for i in range(dataOut.NSHTS):
            y[0]=y[1]=dataOut.range1[i]

        y = y.astype(dtype='float64',order='F')
        three=int(3)
        wl = 3.0
        tion=numpy.zeros(three,order='F',dtype='float32')
        fion=numpy.zeros(three,order='F',dtype='float32')
        nui=numpy.zeros(three,order='F',dtype='float32')
        wion=numpy.zeros(three,order='F',dtype='int32')
        bline=0.0
        #bline=numpy.zeros(1,order='F',dtype='float32')
        my_aux = numpy.ones(dataOut.NSHTS,order='F',dtype='float32')
        acf_Temps = numpy.ones(dataOut.NSHTS,order='F',dtype='float32')*numpy.nan
        acf_no_Temps = numpy.ones(dataOut.NSHTS,order='F',dtype='float32')*numpy.nan

        from scipy import signal

        #def func(params):
        #    return (ratio2-self.gaussian(dataOut.heightList[:dataOut.NSHTS],params[0],params[1],params[2]))

        #print("Before loop")
        dataOut.info2[0] = 1
        for i in range(dataOut.NSHTS):
            if dataOut.info2[i]==1:
                angle=dataOut.thb[i]*0.01745
                nue=nui[0]=nui[1]=nui[2]=0.0#nui[3]=0.0
                wion[0]=16 #O
                wion[1]=1 #H
                wion[2]=4 #He
                tion[0]=tion[1]=tion[2]=dataOut.ti2[i]
                #tion[0]=tion[1]=tion[2]=ti2_smooth[i]
                fion[0]=1.0-dataOut.phy2[i] #1
                fion[1]=dataOut.phy2[i] #0
                fion[2]=0.0 #0
                for j in range(dataOut.DPL):
                    tau=dataOut.alag[j]*1.0e-3
                    #print("***********ACF2***********")
                    with suppress_stdout_stderr():#The smoothness in range of "y" depends on the smoothness of the input parameters
                        y[j]=fitacf_acf2.acf2(wl,tau,dataOut.te2[i],tion,fion,nue,nui,wion,angle,dataOut.ph2[i],dataOut.bfm[i],y[j],three)
                        #y[j]=fitacf_acf2.acf2(wl,tau,te2_smooth[i],tion,fion,nue,nui,wion,angle,dataOut.ph2[i],dataOut.bfm[i],y[j],three)
                    #print("i: ", i, "j: ", j)
                    #if i == 0 and j == 1:
                        #print("Params: ",wl,tau,dataOut.te2[i],tion,fion,nue,nui,wion,angle,dataOut.ph2[i],dataOut.bfm[i],y[j],three)
                    #y[j]=fitacf_acf2.acf2(wl,tau,my_te2[i],tion,fion,nue,nui,wion,angle,dataOut.ph2[i],dataOut.bfm[i],y[j],three)
                        #exit(1)
                #if dataOut.ut_Faraday>11.0 and dataOut.range1[i]>150.0 and dataOut.range1[i]<400.0:

                if dataOut.ut_Faraday>11.0 and dataOut.range1[i]>150.0 and dataOut.range1[i]<300.0:
                #if dataOut.ut_Faraday>11.0 and dataOut.range1[i]>150.0 and dataOut.range1[i]<400.0:
                    tau=0.0
                    with suppress_stdout_stderr():
                        bline=fitacf_acf2.acf2(wl,tau,tion,tion,fion,nue,nui,wion,angle,dataOut.ph2[i],dataOut.bfm[i],bline,three)

                    #if i == 0 and j == 1:
                        #print("bline: ",bline)
                    #y[j]=fitacf_acf2.acf2(wl,tau,my_te2[i],tion,fion,nue,nui,wion,angle,dataOut.ph2[i],dataOut.bfm[i],y[j],three)
                        #exit(1)
                    cf=min(1.2,max(1.0,bline/y[0])) #FACTOR DE EFICIENCIA
                    #cf = bline/y[0]
                    #cf=min(2.,max(1.0,bline/y[0]))
                    my_aux[i] = cf
                    acf_Temps[i] = y[0]
                    acf_no_Temps[i] = bline
                    #dataOut.ph2[i]=cf*dataOut.ph2[i] #Instead we adjust the curve "cf" into a Gaussian,
                    #dataOut.sdp2[i]=cf*dataOut.sdp2[i] #in order to get smoother values of density
                for j in range(1,dataOut.DPL):
                    #y[j]=(y[j]/y[0])*dataOut.DH+dataOut.range1[i]
                    y[j]=min(max((y[j]/y[0]),-1.0),1.0)*dataOut.DH+dataOut.range1[i]
                y[0]=dataOut.range1[i]+dataOut.DH


        ratio = my_aux-1
        #ratio = dataOut.te2[:dataOut.NSHTS]/dataOut.ti2[:dataOut.NSHTS]
        def lsq_func(params):
            return (ratio-self.gaussian(dataOut.heightList[:dataOut.NSHTS],params[0],params[1],params[2]))

        x0_value = numpy.array([max(ratio),250,20])

        popt = least_squares(lsq_func,x0=x0_value,verbose=0)

        A = popt.x[0]; B = popt.x[1]; C = popt.x[2]

        aux = self.gaussian(dataOut.heightList[:dataOut.NSHTS], A, B, C) + 1 #ratio + 1

        dataOut.ph2[:dataOut.NSHTS]*=aux
        dataOut.sdp2[:dataOut.NSHTS]*=aux
        #dataOut.ph2[:26]*=aux[:26]
        #dataOut.sdp2[:26]*=aux[:26]
        #print(aux)
        #print("inside correction",dataOut.ph2)

    def run(self,dataOut,savecf=0):
        #print("hour",gmtime(dataOut.utctime).tm_hour)
        if gmtime(dataOut.utctime).tm_hour < 24. and gmtime(dataOut.utctime).tm_hour >= 11.:
            if hasattr(dataOut, 'flagSpreadF') and dataOut.flagSpreadF:
                pass
            else:
                #print("inside")
                self.TeTiEstimation(dataOut)
                dataOut.flagTeTiCorrection = True
                self.normalize(dataOut)
        #'''
        #Here save dataOut.cf
        if savecf:
            try:
                import pandas as pd
                if self.csv_flag:
                    if not os.path.exists("./cf"):
                        os.makedirs("./cf")
                    self.doy_csv = datetime.datetime.fromtimestamp(dataOut.utctime).strftime('%j')
                    self.year_csv = datetime.datetime.fromtimestamp(dataOut.utctime).strftime('%Y')
                file = open("./cf/cf{0}{1}.csv".format(self.year_csv,self.doy_csv), "x")
                f = csv.writer(file)
                f.writerow(numpy.array(["timestamp", "time_string", 'cf']))  # Added time_string header
                self.csv_flag = 0
                print("Creating cf File")
                print("Writing cf File")
            except:
                file = open("./cf/cf{0}{1}.csv".format(self.year_csv,self.doy_csv), "a")
                f = csv.writer(file)
                print("Writing cf File")
            
            # Create time string from timestamp
            time_str = datetime.datetime.fromtimestamp(dataOut.utctime).strftime('%Y-%m-%d %H:%M:%S')
            cf = numpy.array([dataOut.utctime, time_str, dataOut.cf])  # Added time string
            f.writerow(cf)
        #'''

        return dataOut

class DataPlotCleaner(Operation):
    '''
    Written by R. Flores
    '''
    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

    def run(self,dataOut):

        THRESH_MIN_POW=10000
        THRESH_MAX_POW=10000000
        THRESH_MIN_TEMP=500
        THRESH_MAX_TEMP=4000
        dataOut.DensityClean=numpy.zeros((1,dataOut.NDP))
        dataOut.EDensityClean=numpy.zeros((1,dataOut.NDP))
        dataOut.ElecTempClean=numpy.zeros((1,dataOut.NDP))
        dataOut.EElecTempClean=numpy.zeros((1,dataOut.NDP))
        dataOut.IonTempClean=numpy.zeros((1,dataOut.NDP))
        dataOut.EIonTempClean=numpy.zeros((1,dataOut.NDP))

        dataOut.DensityClean[0]=numpy.copy(dataOut.ph2)
        dataOut.EDensityClean[0]=numpy.copy(dataOut.sdp2)
        dataOut.ElecTempClean[0,:dataOut.NSHTS]=numpy.copy(dataOut.te2)
        dataOut.EElecTempClean[0,:dataOut.NSHTS]=numpy.copy(dataOut.ete2)
        dataOut.IonTempClean[0,:dataOut.NSHTS]=numpy.copy(dataOut.ti2)
        dataOut.EIonTempClean[0,:dataOut.NSHTS]=numpy.copy(dataOut.eti2)

        for i in range(dataOut.NDP):
            if dataOut.DensityClean[0,i]<THRESH_MIN_POW:
                dataOut.DensityClean[0,i]=THRESH_MIN_POW

        for i in range(dataOut.NDP):
            if dataOut.DensityClean[0,i]>THRESH_MAX_POW:
                dataOut.DensityClean[0,i]=THRESH_MAX_POW

        for i in range(dataOut.NSHTS):
            dataOut.ElecTempClean[0,i]=(max(1.0, dataOut.ElecTempClean[0,i]))
            dataOut.IonTempClean[0,i]=(max(1.0, dataOut.IonTempClean[0,i]))
        for i in range(dataOut.NSHTS):
            if dataOut.ElecTempClean[0,i]<THRESH_MIN_TEMP:
                dataOut.ElecTempClean[0,i]=THRESH_MIN_TEMP
            if dataOut.IonTempClean[0,i]<THRESH_MIN_TEMP:
                dataOut.IonTempClean[0,i]=THRESH_MIN_TEMP
        for i in range(dataOut.NSHTS):
            if dataOut.ElecTempClean[0,i]>THRESH_MAX_TEMP:
                dataOut.ElecTempClean[0,i]=THRESH_MAX_TEMP
            if dataOut.IonTempClean[0,i]>THRESH_MAX_TEMP:
                dataOut.IonTempClean[0,i]=THRESH_MAX_TEMP
        for i in range(dataOut.NSHTS):
            if dataOut.EElecTempClean[0,i]>500:#
                dataOut.ElecTempClean[0,i]=500
            if dataOut.EIonTempClean[0,i]>500:#
                dataOut.IonTempClean[0,i]=500

        missing=numpy.nan

        for i in range(dataOut.NSHTS,dataOut.NDP):

            dataOut.ElecTempClean[0,i]=missing
            dataOut.EElecTempClean[0,i]=missing
            dataOut.IonTempClean[0,i]=missing
            dataOut.EIonTempClean[0,i]=missing

        return dataOut


class DataSaveCleaner(Operation):
    '''
    Written by R. Flores
    '''
    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)
        self.csv_flag = 1

    def run(self,dataOut,bypass=False,savecfclean=0):
        #print(dataOut.heightList)
        #exit(1)
        dataOut.DensityFinal=numpy.zeros((1,dataOut.NDP))
        dataOut.dphiFinal=numpy.zeros((1,dataOut.NDP))
        dataOut.EDensityFinal=numpy.zeros((1,dataOut.NDP))
        dataOut.ElecTempFinal=numpy.zeros((1,dataOut.NDP))
        dataOut.EElecTempFinal=numpy.zeros((1,dataOut.NDP))
        dataOut.IonTempFinal=numpy.zeros((1,dataOut.NDP))
        dataOut.EIonTempFinal=numpy.zeros((1,dataOut.NDP))
        dataOut.PhyFinal=numpy.zeros((1,dataOut.NDP))
        dataOut.EPhyFinal=numpy.zeros((1,dataOut.NDP))

        dataOut.DensityFinal[0]=numpy.copy(dataOut.ph2)
        dataOut.dphiFinal[0]=numpy.copy(dataOut.dphi)
        dataOut.EDensityFinal[0]=numpy.copy(dataOut.sdp2)
        dataOut.ElecTempFinal[0,:dataOut.NSHTS]=numpy.copy(dataOut.te2)
        dataOut.EElecTempFinal[0,:dataOut.NSHTS]=numpy.copy(dataOut.ete2)
        dataOut.IonTempFinal[0,:dataOut.NSHTS]=numpy.copy(dataOut.ti2)
        dataOut.EIonTempFinal[0,:dataOut.NSHTS]=numpy.copy(dataOut.eti2)
        dataOut.PhyFinal[0,:dataOut.NSHTS]=numpy.copy(dataOut.phy2)
        dataOut.EPhyFinal[0,:dataOut.NSHTS]=numpy.copy(dataOut.ephy2)

        # Bypass
        if bypass == True:
            # for MAD
            dataOut.flagNoData = numpy.all(numpy.isnan(dataOut.DensityFinal)) #Si todos los valores son NaN no se prosigue
            '''Save of clean data information in temp csv file for den correction'''
            if not dataOut.flagNoData:
                if savecfclean:
                    try:
                        import pandas as pd
                        if self.csv_flag:
                            if not os.path.exists("./cfclean"):
                                os.makedirs("./cfclean")
                            self.doy_csv = datetime.datetime.fromtimestamp(dataOut.utctime).strftime('%j')
                            self.year_csv = datetime.datetime.fromtimestamp(dataOut.utctime).strftime('%Y')
                        file = open("./cfclean/cfclean{0}{1}.csv".format(self.year_csv,self.doy_csv), "x")
                        f = csv.writer(file)
                        f.writerow(numpy.array(["timestamp",'cf']))
                        self.csv_flag = 0
                        print("Creating cf clean File")
                        print("Writing cf clean File")
                    except:
                        file = open("./cfclean/cfclean{0}{1}.csv".format(self.year_csv,self.doy_csv), "a")
                        f = csv.writer(file)
                        print("Writing cf clean File")
                    cf = numpy.array([dataOut.utctime,dataOut.cf])
                    f.writerow(cf)
                    file.close()
            # for plot
            #dataOut.flagNoData = False #Descomentar solo para ploteo #Comentar para MADWriter

            dataOut.DensityFinal *= 1.e6 #Convert units to m^⁻3
            dataOut.EDensityFinal *= 1.e6 #Convert units to m^⁻3
            print("Cleaning bypassed")
            return dataOut

        missing=numpy.nan
        #print("den1: ",dataOut.DensityFinal)
        temp_min=100.0
        temp_max=3000.0#6000.0e
        #print("Density: ",dataOut.DensityFinal[0])
        #print("Error: ",dataOut.EDensityFinal[0])
        #print(100*dataOut.EDensityFinal[0]/dataOut.DensityFinal[0])
        den_err_percent = 100*dataOut.EDensityFinal[0]/dataOut.DensityFinal[0]
        max_den_err_per = 35#30 #Densidades con error mayor al 35% se setean en NaN
        for i in range(dataOut.NSHTS):

            if den_err_percent[i] >= max_den_err_per:
                dataOut.DensityFinal[0,i]=dataOut.EDensityFinal[0,i]=missing
                if i > 40: #Alturas mayores que 600
                    dataOut.DensityFinal[0,i:]=dataOut.EDensityFinal[0,i:]=missing
            # Additional info - number of zeros, temp limit
            if dataOut.info2[i]!=1:
                dataOut.ElecTempFinal[0,i]=dataOut.EElecTempFinal[0,i]=dataOut.IonTempFinal[0,i]=dataOut.EIonTempFinal[0,i]=missing

            #Limit electron temperature range
            if dataOut.ElecTempFinal[0,i]<=temp_min or dataOut.ElecTempFinal[0,i]>temp_max or dataOut.EElecTempFinal[0,i]>temp_max:
                dataOut.ElecTempFinal[0,i]=dataOut.EElecTempFinal[0,i]=missing

            #Limit ion temperature range
            if dataOut.IonTempFinal[0,i]<=temp_min or dataOut.IonTempFinal[0,i]>temp_max or dataOut.EIonTempFinal[0,i]>temp_max:
                dataOut.IonTempFinal[0,i]=dataOut.EIonTempFinal[0,i]=missing

            #If there are 5 or less acf points to plot
            if dataOut.lags_to_plot[i,:][~numpy.isnan(dataOut.lags_to_plot[i,:])].shape[0]<6:
                dataOut.ElecTempFinal[0,i]=dataOut.EElecTempFinal[0,i]=dataOut.IonTempFinal[0,i]=dataOut.EIonTempFinal[0,i]=missing
            
            #If exists at least one big acf error in the morning
            if dataOut.ut_Faraday>4 and dataOut.ut_Faraday<11:
                if numpy.nanmax(dataOut.acfs_error_to_plot[i,:])>=10: #10 (max is dH -> hardcoded)
                    dataOut.ElecTempFinal[0,i]=dataOut.EElecTempFinal[0,i]=dataOut.IonTempFinal[0,i]=dataOut.EIonTempFinal[0,i]=missing

            if dataOut.EPhyFinal[0,i]<0.0 or dataOut.EPhyFinal[0,i]>1.0:
                dataOut.PhyFinal[0,i]=dataOut.EPhyFinal[0,i]=missing
            #'''
            if dataOut.EDensityFinal[0,i]>0.0 and dataOut.DensityFinal[0,i]>0.0 and dataOut.DensityFinal[0,i]<9.9e6:
                dataOut.EDensityFinal[0,i]=max(dataOut.EDensityFinal[0,i],1000.0)
            else:
                dataOut.DensityFinal[0,i]=dataOut.EDensityFinal[0,i]=missing
                #'''

            if dataOut.PhyFinal[0,i]==0 or dataOut.PhyFinal[0,i]>0.4:
                dataOut.PhyFinal[0,i]=dataOut.EPhyFinal[0,i]=missing
            if dataOut.ElecTempFinal[0,i]==dataOut.IonTempFinal[0,i]:
                dataOut.EElecTempFinal[0,i]=dataOut.EIonTempFinal[0,i]
            if numpy.isnan(dataOut.ElecTempFinal[0,i]):
                dataOut.EElecTempFinal[0,i]=missing
            if numpy.isnan(dataOut.IonTempFinal[0,i]):
                dataOut.EIonTempFinal[0,i]=missing
            if numpy.isnan(dataOut.ElecTempFinal[0,i]) or numpy.isnan(dataOut.EElecTempFinal[0,i]):
                dataOut.ElecTempFinal[0,i]=dataOut.EElecTempFinal[0,i]=dataOut.IonTempFinal[0,i]=dataOut.EIonTempFinal[0,i]=missing

        for i in range(12,dataOut.NSHTS-1):

            if numpy.isnan(dataOut.ElecTempFinal[0,i-1]) and numpy.isnan(dataOut.ElecTempFinal[0,i+1]):
                dataOut.ElecTempFinal[0,i]=dataOut.EElecTempFinal[0,i]=missing

            if numpy.isnan(dataOut.IonTempFinal[0,i-1]) and numpy.isnan(dataOut.IonTempFinal[0,i+1]):
                dataOut.IonTempFinal[0,i]=dataOut.EIonTempFinal[0,i]=missing

            if dataOut.ut_Faraday>4 and dataOut.ut_Faraday<11:

                if numpy.isnan(dataOut.ElecTempFinal[0,i-1]) and numpy.isnan(dataOut.ElecTempFinal[0,i-2]) and numpy.isnan(dataOut.ElecTempFinal[0,i+2]) and numpy.isnan(dataOut.ElecTempFinal[0,i+3]): #and numpy.isnan(dataOut.ElecTempFinal[0,i-5]):

                    dataOut.ElecTempFinal[0,i]=dataOut.EElecTempFinal[0,i]=missing
                if numpy.isnan(dataOut.IonTempFinal[0,i-1]) and numpy.isnan(dataOut.IonTempFinal[0,i-2]) and numpy.isnan(dataOut.IonTempFinal[0,i+2]) and numpy.isnan(dataOut.IonTempFinal[0,i+3]): #and numpy.isnan(dataOut.IonTempFinal[0,i-5]):

                    dataOut.IonTempFinal[0,i]=dataOut.EIonTempFinal[0,i]=missing



            if i>25:
                if numpy.isnan(dataOut.ElecTempFinal[0,i-1]) and numpy.isnan(dataOut.ElecTempFinal[0,i-2]) and numpy.isnan(dataOut.ElecTempFinal[0,i-3]) and numpy.isnan(dataOut.ElecTempFinal[0,i-4]): #and numpy.isnan(dataOut.ElecTempFinal[0,i-5]):
                    dataOut.ElecTempFinal[0,i]=dataOut.EElecTempFinal[0,i]=missing
                if numpy.isnan(dataOut.IonTempFinal[0,i-1]) and numpy.isnan(dataOut.IonTempFinal[0,i-2]) and numpy.isnan(dataOut.IonTempFinal[0,i-3]) and numpy.isnan(dataOut.IonTempFinal[0,i-4]): #and numpy.isnan(dataOut.IonTempFinal[0,i-5]):

                    dataOut.IonTempFinal[0,i]=dataOut.EIonTempFinal[0,i]=missing

            if numpy.isnan(dataOut.ElecTempFinal[0,i]) or numpy.isnan(dataOut.EElecTempFinal[0,i]):

                dataOut.ElecTempFinal[0,i]=dataOut.EElecTempFinal[0,i]=dataOut.IonTempFinal[0,i]=dataOut.EIonTempFinal[0,i]=missing

        for i in range(12,dataOut.NSHTS-1):

            if numpy.isnan(dataOut.ElecTempFinal[0,i-1]) and numpy.isnan(dataOut.ElecTempFinal[0,i+1]):
                dataOut.ElecTempFinal[0,i]=dataOut.EElecTempFinal[0,i]=missing

            if numpy.isnan(dataOut.IonTempFinal[0,i-1]) and numpy.isnan(dataOut.IonTempFinal[0,i+1]):
                dataOut.IonTempFinal[0,i]=dataOut.EIonTempFinal[0,i]=missing

            if numpy.isnan(dataOut.DensityFinal[0,i-1]) and numpy.isnan(dataOut.DensityFinal[0,i+1]): ##NEW
                dataOut.DensityFinal[0,i]=dataOut.EDensityFinal[0,i]=missing ##NEW

            if numpy.isnan(dataOut.ElecTempFinal[0,i]) or numpy.isnan(dataOut.EElecTempFinal[0,i]):

                dataOut.ElecTempFinal[0,i]=dataOut.EElecTempFinal[0,i]=dataOut.IonTempFinal[0,i]=dataOut.EIonTempFinal[0,i]=missing

        if numpy.count_nonzero(~numpy.isnan(dataOut.ElecTempFinal[0,12:50]))<5:
            dataOut.ElecTempFinal[0,:]=dataOut.EElecTempFinal[0,:]=missing
        if numpy.count_nonzero(~numpy.isnan(dataOut.IonTempFinal[0,12:50]))<5:
            dataOut.IonTempFinal[0,:]=dataOut.EIonTempFinal[0,:]=missing


        '''
        nanindex = numpy.argwhere(numpy.isnan(dataOut.DensityFinal[0,:]))
        #print(nanindex)
        i1 = nanindex[-1][0]
        #print("i1",i1)
        dataOut.EDensityFinal[0,i1:] = dataOut.DensityFinal[0,i1:] = missing
        '''

        '''
        #print("den2: ",dataOut.DensityFinal)
        if gmtime(dataOut.utctime).tm_hour >= 23. or gmtime(dataOut.utctime).tm_hour < 5.: #18-00 LT
            nanindex = numpy.argwhere(numpy.isnan(dataOut.DensityFinal[0,:33]))
            #print(nanindex)
            i1 = nanindex[-1][0]
            #print("i1",i1)
            dataOut.EDensityFinal[0,:i1] = dataOut.DensityFinal[0,:i1] = missing
            #print("den3: ",dataOut.DensityFinal)
        elif gmtime(dataOut.utctime).tm_hour >= 6. or gmtime(dataOut.utctime).tm_hour < 11.: #18-00 LT
            nanindex = numpy.argwhere(numpy.isnan(dataOut.DensityFinal[0,:20]))
            #print(nanindex)
            i1 = nanindex[-1][0]
            #print("i1",i1)
            dataOut.EDensityFinal[0,:i1] = dataOut.DensityFinal[0,:i1] = missing
            '''
        #print("den_nans: ",dataOut.DensityFinal[0,12:50])
        if numpy.count_nonzero(~numpy.isnan(dataOut.DensityFinal[0,12:50]))<=5:
            dataOut.DensityFinal[0,:]=dataOut.EDensityFinal[0,:]=missing
        #for i in range(dataOut.NSHTS,dataOut.NDP):
        #for i in range(40,dataOut.NDP):
        #print("den2: ",dataOut.DensityFinal)
        dataOut.DensityFinal[0,dataOut.NSHTS:]=missing
        dataOut.EDensityFinal[0,dataOut.NSHTS:]=missing
        dataOut.ElecTempFinal[0,dataOut.NSHTS:]=missing
        dataOut.EElecTempFinal[0,dataOut.NSHTS:]=missing
        dataOut.IonTempFinal[0,dataOut.NSHTS:]=missing
        dataOut.EIonTempFinal[0,dataOut.NSHTS:]=missing
        dataOut.PhyFinal[0,dataOut.NSHTS:]=missing
        dataOut.EPhyFinal[0,dataOut.NSHTS:]=missing

        #To be tested
        '''
        nanindex = numpy.argwhere(numpy.isnan(dataOut.DensityFinal))
        i1 = nanindex[-1][0]
        print("i1",i1)
        dataOut.DensityFinal[0,i1+1:]=missing
        dataOut.EDensityFinal[0,i1+1:]=missing
        dataOut.ElecTempFinal[0,i1+1:]=missing
        dataOut.EElecTempFinal[0,i1+1:]=missing
        dataOut.IonTempFinal[0,i1+1:]=missing
        dataOut.EIonTempFinal[0,i1+1:]=missing
        dataOut.PhyFinal[0,i1+1:]=missing
        dataOut.EPhyFinal[0,i1+1:]=missing
        '''
        #'''
        if gmtime(dataOut.utctime).tm_hour >= 13. and gmtime(dataOut.utctime).tm_hour < 21.: #07-16 LT
            dataOut.DensityFinal[0,:13]=missing
            dataOut.EDensityFinal[0,:13]=missing
            dataOut.ElecTempFinal[0,:13]=missing
            dataOut.EElecTempFinal[0,:13]=missing
            dataOut.IonTempFinal[0,:13]=missing
            dataOut.EIonTempFinal[0,:13]=missing
            dataOut.PhyFinal[0,:13]=missing
            dataOut.EPhyFinal[0,:13]=missing
            #'''
        else:
            if gmtime(dataOut.utctime).tm_hour == 9 and gmtime(dataOut.utctime).tm_min == 20:
                pass
            else:
                dataOut.DensityFinal[0,:dataOut.min_id_eej+1]=missing
                dataOut.EDensityFinal[0,:dataOut.min_id_eej+1]=missing
                dataOut.ElecTempFinal[0,:dataOut.min_id_eej+1]=missing
                dataOut.EElecTempFinal[0,:dataOut.min_id_eej+1]=missing
                dataOut.IonTempFinal[0,:dataOut.min_id_eej+1]=missing
                dataOut.EIonTempFinal[0,:dataOut.min_id_eej+1]=missing
                dataOut.PhyFinal[0,:dataOut.min_id_eej+1]=missing
                dataOut.EPhyFinal[0,:dataOut.min_id_eej+1]=missing
        '''
        if gmtime(dataOut.utctime).tm_hour >= 11. or gmtime(dataOut.utctime).tm_hour < 23.: #06-18 LT
            dataOut.DensityFinal[0,:13]=missing
            dataOut.EDensityFinal[0,:13]=missing
            dataOut.ElecTempFinal[0,:13]=missing
            dataOut.EElecTempFinal[0,:13]=missing
            dataOut.IonTempFinal[0,:13]=missing
            dataOut.EIonTempFinal[0,:13]=missing
            dataOut.PhyFinal[0,:13]=missing
            dataOut.EPhyFinal[0,:13]=missing
        else:
        #for i in range(12):
            dataOut.DensityFinal[0,:12]=missing
            dataOut.EDensityFinal[0,:12]=missing
            dataOut.ElecTempFinal[0,:12]=missing
            dataOut.EElecTempFinal[0,:12]=missing
            dataOut.IonTempFinal[0,:12]=missing
            dataOut.EIonTempFinal[0,:12]=missing
            dataOut.PhyFinal[0,:12]=missing
            dataOut.EPhyFinal[0,:12]=missing
            '''
        #print(dataOut.EDensityFinal)
        #exit(1)
        '''
        print(dataOut.ElecTempFinal)
        print(dataOut.heightList)
        exit(1)
        '''
         ### Manual Data Cleaning
        
        time_text = datetime.datetime.utcfromtimestamp(dataOut.utctime)
        DOY = time_text.timetuple().tm_yday
        '''
        # 2 Sep 24
        if (time_text.hour == 1 and (time_text.minute == 40)): #06-18 LT
            #dataOut.DensityFinal[0,27:]=missing
            #dataOut.EDensityFinal[0,27:]=missing
            dataOut.ElecTempFinal[0,:]=missing
            dataOut.EElecTempFinal[0,:]=missing
            dataOut.IonTempFinal[0,:]=missing
            dataOut.EIonTempFinal[0,:]=missing
            dataOut.PhyFinal[0,:]=missing
            dataOut.EPhyFinal[0,:] = missing
        if (time_text.hour == 5 and (time_text.minute >= 20)) or (time_text.hour == 6 and (time_text.minute <= 5)): #06-18 LT
            #dataOut.DensityFinal[0,27:]=missing
            #dataOut.EDensityFinal[0,27:]=missing
            dataOut.ElecTempFinal[0,27:]=missing
            dataOut.EElecTempFinal[0,27:]=missing
            dataOut.IonTempFinal[0,27:]=missing
            dataOut.EIonTempFinal[0,27:]=missing
            dataOut.PhyFinal[0,27:]=missing
            dataOut.EPhyFinal[0, 27:] = missing
        '''
        '''
        # 5 Sep 24
        
        if (time_text.hour == 3 and (time_text.minute == 35)) or (time_text.hour == 3 and (time_text.minute == 20)): #06-18 LT
            #dataOut.DensityFinal[0,27:]=missing
            #dataOut.EDensityFinal[0,27:]=missing
            dataOut.ElecTempFinal[0,:]=missing
            dataOut.EElecTempFinal[0,:]=missing
            dataOut.IonTempFinal[0,:]=missing
            dataOut.EIonTempFinal[0,:]=missing
            dataOut.PhyFinal[0,:]=missing
            dataOut.EPhyFinal[0,:] = missing
        
        if (time_text.hour == 5 and (time_text.minute >= 20)) or (time_text.hour == 6 and (time_text.minute <= 5)): #06-18 LT
            #dataOut.DensityFinal[0,27:]=missing
            #dataOut.EDensityFinal[0,27:]=missing
            dataOut.ElecTempFinal[0,27:]=missing
            dataOut.EElecTempFinal[0,27:]=missing
            dataOut.IonTempFinal[0,27:]=missing
            dataOut.EIonTempFinal[0,27:]=missing
            dataOut.PhyFinal[0,27:]=missing
            dataOut.EPhyFinal[0,27:] = missing
        '''
        '''
        # 6 Sep 24
        if (time_text.hour == 5 and (time_text.minute <= 25)) or (time_text.hour == 6 and (time_text.minute >= 5)): #06-18 LT
            #dataOut.DensityFinal[0,27:]=missing
            #dataOut.EDensityFinal[0,27:]=missing
            dataOut.ElecTempFinal[0,25:]=missing
            dataOut.EElecTempFinal[0,25:]=missing
            dataOut.IonTempFinal[0,25:]=missing
            dataOut.EIonTempFinal[0,25:]=missing
            dataOut.PhyFinal[0,25:]=missing
            dataOut.EPhyFinal[0, 25:] = missing
        '''    
        '''# 8 Sep 24
        if True: #06-18 LT
            #dataOut.DensityFinal[0,27:]=missing
            #dataOut.EDensityFinal[0,27:]=missing
            dataOut.ElecTempFinal[0,36:]=missing
            dataOut.EElecTempFinal[0,36:]=missing
            dataOut.IonTempFinal[0,36:]=missing
            dataOut.EIonTempFinal[0,36:]=missing
            dataOut.PhyFinal[0,36:]=missing
            dataOut.EPhyFinal[0, 36:] = missing'''
        '''# 24 Jan 25
        if (time_text.hour >= 5 ) and (time_text.hour <= 7):
            #dataOut.DensityFinal[0,27:]=missing
            #dataOut.EDensityFinal[0,27:]=missing
            dataOut.ElecTempFinal[0,:]=missing
            dataOut.EElecTempFinal[0,:]=missing
            dataOut.IonTempFinal[0,:]=missing
            dataOut.EIonTempFinal[0,:]=missing
            dataOut.PhyFinal[0,:]=missing
            dataOut.EPhyFinal[0, :] = missing'''
        start = time()
        flagcleandata = True

        print('dataOut.EDensityFinal.shape', dataOut.EDensityFinal.shape)
        print(dataOut.EDensityFinal)

        if flagcleandata:
            #print("Final Cleaning Process", time_text.hour, time_text.minute)
            path = os.path.join(os.path.dirname(__file__), 'clean_data.json')
            with open(path) as f:
                jsondata= json.load(f)

            corrections = {}
            

            for condition in jsondata['conditions']:
                year = condition['year']
                doy = condition['doy']
                init = condition['initial_time']
                final = condition['final_time']
                aux_index = condition['aux_index']

                only = condition.get('only')

                input_time_obj = datetime.time(hour=time_text.hour, minute=time_text.minute)
                init_time_obj = datetime.time(hour=init[0], minute=init[1])
                final_time_obj = datetime.time(hour=final[0], minute=final[1])

                is_between = init_time_obj <= input_time_obj <= final_time_obj
                
                if (year != time_text.year) or (DOY != doy) or (is_between == False):
                    #print("NON valid condition:", condition)
                    continue
                
                print("valid condition:", condition)
                indexi, indexf = aux_index[0], aux_index[1]
                #index = slice(indexi, indexf)
                #print(index, "index")
                print(only)
                if only == "te":
                    print("entered only")
                    dataOut.ElecTempFinal[0,indexi:indexf]=missing
                    dataOut.EElecTempFinal[0,indexi:indexf]=missing
                    dataOut.IonTempFinal[0,indexi:indexf]=missing
                    dataOut.EIonTempFinal[0,indexi:indexf]=missing
                    dataOut.PhyFinal[0,indexi:indexf]=missing
                    dataOut.EPhyFinal[0,indexi:indexf]=missing
                else:
                    dataOut.DensityFinal[0,indexi:indexf]=missing
                    print('dataOut.DensityFinal.shape', dataOut.DensityFinal.shape)
                    print(dataOut.DensityFinal)
                    dataOut.EDensityFinal[0,indexi:indexf]=missing
                    print('dataOut.EDensityFinal.shape', dataOut.EDensityFinal.shape)
                    print(dataOut.EDensityFinal)
                    dataOut.ElecTempFinal[0,indexi:indexf]=missing
                    dataOut.EElecTempFinal[0,indexi:indexf]=missing
                    dataOut.IonTempFinal[0,indexi:indexf]=missing
                    dataOut.EIonTempFinal[0,indexi:indexf]=missing
                    dataOut.PhyFinal[0,indexi:indexf]=missing
                    dataOut.EPhyFinal[0,indexi:indexf]=missing
                
                print(f"** Cleaning applied ** Data eliminated at {time_text} from heigh index {indexi} to {indexf}")

        end = time() ########
        #spend_clean_time += end - start
        print("clean data time: ", end - start) 

        #print("den_final",dataOut.DensityFinal)

        # for MAD
        dataOut.flagNoData = numpy.all(numpy.isnan(dataOut.DensityFinal)) #Si todos los valores son NaN no se prosigue

        '''Save of clean data information in temp csv file for den correction'''
        if not dataOut.flagNoData:
            if savecfclean:
                try:
                    import pandas as pd
                    if self.csv_flag:
                        if not os.path.exists("./cfclean"):
                            os.makedirs("./cfclean")
                        self.doy_csv = datetime.datetime.fromtimestamp(dataOut.utctime).strftime('%j')
                        self.year_csv = datetime.datetime.fromtimestamp(dataOut.utctime).strftime('%Y')
                    file = open("./cfclean/cfclean{0}{1}.csv".format(self.year_csv,self.doy_csv), "x")
                    f = csv.writer(file)
                    f.writerow(numpy.array(["timestamp",'cf']))
                    self.csv_flag = 0
                    print("Creating cf clean File")
                    print("Writing cf clean File")
                except:
                    file = open("./cfclean/cfclean{0}{1}.csv".format(self.year_csv,self.doy_csv), "a")
                    f = csv.writer(file)
                    print("Writing cf clean File")
                cf = numpy.array([dataOut.utctime,dataOut.cf])
                f.writerow(cf)
                file.close()
        # for plot
        #dataOut.flagNoData = False #Descomentar solo para ploteo #Comentar para MADWriter

        dataOut.DensityFinal *= 1.e6 #Convert units to m^⁻3
        dataOut.EDensityFinal *= 1.e6 #Convert units to m^⁻3
        print("Save Cleaner: ", dataOut.flagNoData)


        #print("den: ", dataOut.DensityFinal[0,27])
        return dataOut

class DataSaveCleanerHP(Operation):
    '''
    Written by R. Flores
    '''
    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

    def run(self,dataOut):

        dataOut.Density_DP=numpy.zeros(dataOut.cut)
        dataOut.EDensity_DP=numpy.zeros(dataOut.cut)
        dataOut.ElecTemp_DP=numpy.zeros(dataOut.cut)
        dataOut.EElecTemp_DP=numpy.zeros(dataOut.cut)
        dataOut.IonTemp_DP=numpy.zeros(dataOut.cut)
        dataOut.EIonTemp_DP=numpy.zeros(dataOut.cut)
        dataOut.Phy_DP=numpy.zeros(dataOut.cut)
        dataOut.EPhy_DP=numpy.zeros(dataOut.cut)
        dataOut.Phe_DP=numpy.empty(dataOut.cut)
        dataOut.EPhe_DP=numpy.empty(dataOut.cut)

        dataOut.Density_DP[:]=numpy.copy(dataOut.ph2[:dataOut.cut])
        dataOut.EDensity_DP[:]=numpy.copy(dataOut.sdp2[:dataOut.cut])
        dataOut.ElecTemp_DP[:]=numpy.copy(dataOut.te2[:dataOut.cut])
        dataOut.EElecTemp_DP[:]=numpy.copy(dataOut.ete2[:dataOut.cut])
        dataOut.IonTemp_DP[:]=numpy.copy(dataOut.ti2[:dataOut.cut])
        dataOut.EIonTemp_DP[:]=numpy.copy(dataOut.eti2[:dataOut.cut])
        dataOut.Phy_DP[:]=numpy.copy(dataOut.phy2[:dataOut.cut])
        dataOut.EPhy_DP[:]=numpy.copy(dataOut.ephy2[:dataOut.cut])
        dataOut.Phe_DP[:]=numpy.nan
        dataOut.EPhe_DP[:]=numpy.nan

        missing=numpy.nan
        temp_min=100.0
        temp_max_dp=3000.0

        for i in range(dataOut.cut):
            if dataOut.info2[i]!=1:
                dataOut.ElecTemp_DP[i]=dataOut.EElecTemp_DP[i]=dataOut.IonTemp_DP[i]=dataOut.EIonTemp_DP[i]=missing

            if dataOut.ElecTemp_DP[i]<=temp_min or dataOut.ElecTemp_DP[i]>temp_max_dp or dataOut.EElecTemp_DP[i]>temp_max_dp:

                dataOut.ElecTemp_DP[i]=dataOut.EElecTemp_DP[i]=missing

            if dataOut.IonTemp_DP[i]<=temp_min or dataOut.IonTemp_DP[i]>temp_max_dp or dataOut.EIonTemp_DP[i]>temp_max_dp:
                dataOut.IonTemp_DP[i]=dataOut.EIonTemp_DP[i]=missing

####################################################################################### CHECK THIS
            if dataOut.lags_to_plot[i,:][~numpy.isnan(dataOut.lags_to_plot[i,:])].shape[0]<6:
                dataOut.ElecTemp_DP[i]=dataOut.EElecTemp_DP[i]=dataOut.IonTemp_DP[i]=dataOut.EIonTemp_DP[i]=missing

            if dataOut.ut_Faraday>4 and dataOut.ut_Faraday<11:
                if numpy.nanmax(dataOut.acfs_error_to_plot[i,:])>=10:
                    dataOut.ElecTemp_DP[i]=dataOut.EElecTemp_DP[i]=dataOut.IonTemp_DP[i]=dataOut.EIonTemp_DP[i]=missing
#######################################################################################

            if dataOut.EPhy_DP[i]<0.0 or dataOut.EPhy_DP[i]>1.0:
                dataOut.Phy_DP[i]=dataOut.EPhy_DP[i]=missing
            if dataOut.EDensity_DP[i]>0.0 and dataOut.Density_DP[i]>0.0 and dataOut.Density_DP[i]<9.9e6:
                dataOut.EDensity_DP[i]=max(dataOut.EDensity_DP[i],1000.0)
            else:
                dataOut.Density_DP[i]=dataOut.EDensity_DP[i]=missing
            if dataOut.Phy_DP[i]==0 or dataOut.Phy_DP[i]>0.4:
                dataOut.Phy_DP[i]=dataOut.EPhy_DP[i]=missing
            if dataOut.ElecTemp_DP[i]==dataOut.IonTemp_DP[i]:
                dataOut.EElecTemp_DP[i]=dataOut.EIonTemp_DP[i]
            if numpy.isnan(dataOut.ElecTemp_DP[i]):
                dataOut.EElecTemp_DP[i]=missing
            if numpy.isnan(dataOut.IonTemp_DP[i]):
                dataOut.EIonTemp_DP[i]=missing
            if numpy.isnan(dataOut.ElecTemp_DP[i]) or numpy.isnan(dataOut.EElecTemp_DP[i]):
                dataOut.ElecTemp_DP[i]=dataOut.EElecTemp_DP[i]=dataOut.IonTemp_DP[i]=dataOut.EIonTemp_DP[i]=missing



        dataOut.Density_LP=numpy.zeros(dataOut.NACF-dataOut.cut)
        dataOut.EDensity_LP=numpy.zeros(dataOut.NACF-dataOut.cut)
        dataOut.ElecTemp_LP=numpy.zeros(dataOut.NACF-dataOut.cut)
        dataOut.EElecTemp_LP=numpy.zeros(dataOut.NACF-dataOut.cut)
        dataOut.IonTemp_LP=numpy.zeros(dataOut.NACF-dataOut.cut)
        dataOut.EIonTemp_LP=numpy.zeros(dataOut.NACF-dataOut.cut)
        dataOut.Phy_LP=numpy.zeros(dataOut.NACF-dataOut.cut)
        dataOut.EPhy_LP=numpy.zeros(dataOut.NACF-dataOut.cut)
        dataOut.Phe_LP=numpy.zeros(dataOut.NACF-dataOut.cut)
        dataOut.EPhe_LP=numpy.zeros(dataOut.NACF-dataOut.cut)

        dataOut.Density_LP[:]=numpy.copy(dataOut.ne[dataOut.cut:dataOut.NACF])
        dataOut.EDensity_LP[:]=numpy.copy(dataOut.ene[dataOut.cut:dataOut.NACF])
        dataOut.ElecTemp_LP[:]=numpy.copy(dataOut.te[dataOut.cut:dataOut.NACF])
        dataOut.EElecTemp_LP[:]=numpy.copy(dataOut.ete[dataOut.cut:dataOut.NACF])
        dataOut.IonTemp_LP[:]=numpy.copy(dataOut.ti[dataOut.cut:dataOut.NACF])
        dataOut.EIonTemp_LP[:]=numpy.copy(dataOut.eti[dataOut.cut:dataOut.NACF])
        dataOut.Phy_LP[:]=numpy.copy(dataOut.ph[dataOut.cut:dataOut.NACF])
        dataOut.EPhy_LP[:]=numpy.copy(dataOut.eph[dataOut.cut:dataOut.NACF])
        dataOut.Phe_LP[:]=numpy.copy(dataOut.phe[dataOut.cut:dataOut.NACF])
        dataOut.EPhe_LP[:]=numpy.copy(dataOut.ephe[dataOut.cut:dataOut.NACF])

        temp_max_lp=6000.0

        for i in range(dataOut.NACF-dataOut.cut):

            if dataOut.ElecTemp_LP[i]<=temp_min or dataOut.ElecTemp_LP[i]>temp_max_lp or dataOut.EElecTemp_LP[i]>temp_max_lp:

                dataOut.ElecTemp_LP[i]=dataOut.EElecTemp_LP[i]=missing

            if dataOut.IonTemp_LP[i]<=temp_min or dataOut.IonTemp_LP[i]>temp_max_lp or dataOut.EIonTemp_LP[i]>temp_max_lp:
                dataOut.IonTemp_LP[i]=dataOut.EIonTemp_LP[i]=missing
            if dataOut.EPhy_LP[i]<0.0 or dataOut.EPhy_LP[i]>1.0:
                dataOut.Phy_LP[i]=dataOut.EPhy_LP[i]=missing

            if dataOut.EPhe_LP[i]<0.0 or dataOut.EPhe_LP[i]>1.0:
                dataOut.Phe_LP[i]=dataOut.EPhe_LP[i]=missing
            if dataOut.EDensity_LP[i]>0.0 and dataOut.Density_LP[i]>0.0 and dataOut.Density_LP[i]<9.9e6 and dataOut.EDensity_LP[i]*dataOut.Density_LP[i]<9.9e6:
                dataOut.EDensity_LP[i]=max(dataOut.EDensity_LP[i],1000.0/dataOut.Density_LP[i])
            else:
                dataOut.Density_LP[i]=missing
                dataOut.EDensity_LP[i]=1.0

            if numpy.isnan(dataOut.Phy_LP[i]):
                dataOut.EPhy_LP[i]=missing

            if numpy.isnan(dataOut.Phe_LP[i]):
                dataOut.EPhe_LP[i]=missing


            if dataOut.ElecTemp_LP[i]==dataOut.IonTemp_LP[i]:
                dataOut.EElecTemp_LP[i]=dataOut.EIonTemp_LP[i]
            if numpy.isnan(dataOut.ElecTemp_LP[i]):
                dataOut.EElecTemp_LP[i]=missing
            if numpy.isnan(dataOut.IonTemp_LP[i]):
                dataOut.EIonTemp_LP[i]=missing
            if numpy.isnan(dataOut.ElecTemp_LP[i]) or numpy.isnan(dataOut.EElecTemp_LP[i]):
                dataOut.ElecTemp_LP[i]=dataOut.EElecTemp_LP[i]=dataOut.IonTemp_LP[i]=dataOut.EIonTemp_LP[i]=missing


        dataOut.DensityFinal=numpy.reshape(numpy.concatenate((dataOut.Density_DP,dataOut.Density_LP)),(1,-1))
        dataOut.EDensityFinal=numpy.reshape(numpy.concatenate((dataOut.EDensity_DP,dataOut.EDensity_LP)),(1,-1))
        dataOut.ElecTempFinal=numpy.reshape(numpy.concatenate((dataOut.ElecTemp_DP,dataOut.ElecTemp_LP)),(1,-1))
        dataOut.EElecTempFinal=numpy.reshape(numpy.concatenate((dataOut.EElecTemp_DP,dataOut.EElecTemp_LP)),(1,-1))
        dataOut.IonTempFinal=numpy.reshape(numpy.concatenate((dataOut.IonTemp_DP,dataOut.IonTemp_LP)),(1,-1))
        dataOut.EIonTempFinal=numpy.reshape(numpy.concatenate((dataOut.EIonTemp_DP,dataOut.EIonTemp_LP)),(1,-1))
        dataOut.PhyFinal=numpy.reshape(numpy.concatenate((dataOut.Phy_DP,dataOut.Phy_LP)),(1,-1))
        dataOut.EPhyFinal=numpy.reshape(numpy.concatenate((dataOut.EPhy_DP,dataOut.EPhy_LP)),(1,-1))
        dataOut.PheFinal=numpy.reshape(numpy.concatenate((dataOut.Phe_DP,dataOut.Phe_LP)),(1,-1))
        dataOut.EPheFinal=numpy.reshape(numpy.concatenate((dataOut.EPhe_DP,dataOut.EPhe_LP)),(1,-1))

        nan_array_2=numpy.empty(dataOut.NACF-dataOut.NDP)
        nan_array_2[:]=numpy.nan

        dataOut.acfs_DP=numpy.zeros((dataOut.NACF,dataOut.DPL),'float32')
        dataOut.acfs_error_DP=numpy.zeros((dataOut.NACF,dataOut.DPL),'float32')
        acfs_dp_aux=dataOut.acfs_to_save.transpose()
        acfs_error_dp_aux=dataOut.acfs_error_to_save.transpose()
        for i in range(dataOut.DPL):
            dataOut.acfs_DP[:,i]=numpy.concatenate((acfs_dp_aux[:,i],nan_array_2))
            dataOut.acfs_error_DP[:,i]=numpy.concatenate((acfs_error_dp_aux[:,i],nan_array_2))
        dataOut.acfs_DP=dataOut.acfs_DP.transpose()
        dataOut.acfs_error_DP=dataOut.acfs_error_DP.transpose()

        dataOut.acfs_LP=numpy.zeros((dataOut.NACF,dataOut.IBITS),'float32')
        dataOut.acfs_error_LP=numpy.zeros((dataOut.NACF,dataOut.IBITS),'float32')

        for i in range(dataOut.NACF):
            for j in range(dataOut.IBITS):
                if numpy.abs(dataOut.errors[j,i]/dataOut.output_LP_integrated.real[0,i,0])<1.0:
                    dataOut.acfs_LP[i,j]=dataOut.output_LP_integrated.real[j,i,0]/dataOut.output_LP_integrated.real[0,i,0]
                    dataOut.acfs_LP[i,j]=max(min(dataOut.acfs_LP[i,j],1.0),-1.0)

                    dataOut.acfs_error_LP[i,j]=dataOut.errors[j,i]/dataOut.output_LP_integrated.real[0,i,0]
                else:
                    dataOut.acfs_LP[i,j]=numpy.nan

                    dataOut.acfs_error_LP[i,j]=numpy.nan

        dataOut.acfs_LP=dataOut.acfs_LP.transpose()
        dataOut.acfs_error_LP=dataOut.acfs_error_LP.transpose()

        dataOut.DensityFinal *= 1.e6 #Convert units to m^⁻3
        dataOut.EDensityFinal *= 1.e6 #Convert units to m^⁻3

         ### Manual Data Cleaning
        
        time_text = datetime.datetime.utcfromtimestamp(dataOut.utctime)
        DOY = time_text.timetuple().tm_yday
        flagcleandata = True
        if flagcleandata:
            #print("Final Cleaning Process", time_text.hour, time_text.minute)
            path = os.path.join(os.path.dirname(__file__), 'clean_data.json')
            with open(path) as f:
                jsondata= json.load(f)

            corrections = {}
            

            for condition in jsondata['conditions']:
                year = condition['year']
                doy = condition['doy']
                init = condition['initial_time']
                final = condition['final_time']
                aux_index = condition['aux_index']

                only = condition.get('only')

                input_time_obj = datetime.time(hour=time_text.hour, minute=time_text.minute)
                init_time_obj = datetime.time(hour=init[0], minute=init[1])
                final_time_obj = datetime.time(hour=final[0], minute=final[1])

                is_between = init_time_obj <= input_time_obj < final_time_obj
                
                if (year != time_text.year) or (DOY != doy) or (is_between == False):
                    #print("NON valid condition:", condition)
                    continue
                
                print("valid condition:", condition)
                indexi, indexf = aux_index[0], aux_index[1]
                #index = slice(indexi, indexf)
                #print(index, "index")
                print(only)
                if only == "te":
                    print("entered only")

                    dataOut.ElecTempFinal[0,indexi:indexf]=missing
                    dataOut.EElecTempFinal[0,indexi:indexf]=missing
                    dataOut.IonTempFinal[0,indexi:indexf]=missing
                    dataOut.EIonTempFinal[0,indexi:indexf]=missing
                    dataOut.PhyFinal[0,indexi:indexf]=missing
                    dataOut.EPhyFinal[0,indexi:indexf]=missing
                    dataOut.PheFinal[0,indexi:indexf]=missing
                    dataOut.EPheFinal[0,indexi:indexf]=missing
                else:
                    dataOut.DensityFinal[0,indexi:indexf]=missing
                    dataOut.EDensityFinal[0,indexi:indexf]=missing
                    dataOut.ElecTempFinal[0,indexi:indexf]=missing
                    dataOut.EElecTempFinal[0,indexi:indexf]=missing
                    dataOut.IonTempFinal[0,indexi:indexf]=missing
                    dataOut.EIonTempFinal[0,indexi:indexf]=missing
                    dataOut.PhyFinal[0,indexi:indexf]=missing
                    dataOut.EPhyFinal[0,indexi:indexf]=missing
                    dataOut.PheFinal[0,indexi:indexf]=missing
                    dataOut.EPheFinal[0,indexi:indexf]=missing
                
                print(f"** Cleaning applied ** Data eliminated at {time_text} from heigh index {indexi} to {indexf}")




        return dataOut

class ACFs(Operation):
    '''
    Written by R. Flores
    '''
    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

        self.aux=1

    def run(self,dataOut):

        if self.aux:
            self.taup=numpy.zeros(dataOut.DPL,'float32')
            self.pacf=numpy.zeros(dataOut.DPL,'float32')
            self.sacf=numpy.zeros(dataOut.DPL,'float32')

            self.taup_full=numpy.zeros(dataOut.DPL,'float32')
            self.pacf_full=numpy.zeros(dataOut.DPL,'float32')
            self.sacf_full=numpy.zeros(dataOut.DPL,'float32')
            self.x_igcej=numpy.zeros(dataOut.DPL,'float32')
            self.y_igcej=numpy.zeros(dataOut.DPL,'float32')
            self.x_ibad=numpy.zeros(dataOut.DPL,'float32')
            self.y_ibad=numpy.zeros(dataOut.DPL,'float32')
            self.aux=0

        dataOut.acfs_to_plot=numpy.zeros((dataOut.NDP,dataOut.DPL),'float32')
        dataOut.acfs_to_save=numpy.zeros((dataOut.NDP,dataOut.DPL),'float32')
        dataOut.acfs_error_to_plot=numpy.zeros((dataOut.NDP,dataOut.DPL),'float32')
        dataOut.acfs_error_to_save=numpy.zeros((dataOut.NDP,dataOut.DPL),'float32')
        dataOut.lags_to_plot=numpy.zeros((dataOut.NDP,dataOut.DPL),'float32')
        dataOut.x_igcej_to_plot=numpy.zeros((dataOut.NDP,dataOut.DPL),'float32')
        dataOut.x_ibad_to_plot=numpy.zeros((dataOut.NDP,dataOut.DPL),'float32')
        dataOut.y_igcej_to_plot=numpy.zeros((dataOut.NDP,dataOut.DPL),'float32')
        dataOut.y_ibad_to_plot=numpy.zeros((dataOut.NDP,dataOut.DPL),'float32')

        for i in range(dataOut.NSHTS):

            acfm=dataOut.rhor[i][0]**2+dataOut.rhoi[i][0]**2

            if acfm>0:
                cc=dataOut.rhor[i][0]/acfm
                ss=dataOut.rhoi[i][0]/acfm
            else:
                cc=1.
                ss=0.

            # keep only uncontaminated data
            for l in range(dataOut.DPL):
                fact=dataOut.DH
                if (dataOut.igcej[i][l]==0 and dataOut.ibad[i][l]==0):

                    self.pacf_full[l]=min(1.0,max(-1.0,(dataOut.rhor[i][l]*cc + dataOut.rhoi[i][l]*ss)))*fact+dataOut.range1[i]
                    self.sacf_full[l]=min(1.0,numpy.sqrt(dataOut.sd[i][l]))*fact
                    self.taup_full[l]=dataOut.alag[l]
                    self.x_igcej[l]=numpy.nan
                    self.y_igcej[l]=numpy.nan
                    self.x_ibad[l]=numpy.nan
                    self.y_ibad[l]=numpy.nan

                else:
                    self.pacf_full[l]=numpy.nan
                    self.sacf_full[l]=numpy.nan
                    self.taup_full[l]=numpy.nan

                    if dataOut.igcej[i][l]:
                        self.x_igcej[l]=dataOut.alag[l]
                        self.y_igcej[l]=dataOut.range1[i]
                        self.x_ibad[l]=numpy.nan
                        self.y_ibad[l]=numpy.nan

                    if dataOut.ibad[i][l]:
                        self.x_igcej[l]=numpy.nan
                        self.y_igcej[l]=numpy.nan
                        self.x_ibad[l]=dataOut.alag[l]
                        self.y_ibad[l]=dataOut.range1[i]

            pacf_new=numpy.copy((self.pacf_full-dataOut.range1[i])/dataOut.DH)
            sacf_new=numpy.copy(self.sacf_full/dataOut.DH)
            dataOut.acfs_to_save[i,:]=numpy.copy(pacf_new)
            dataOut.acfs_error_to_save[i,:]=numpy.copy(sacf_new)
            dataOut.acfs_to_plot[i,:]=numpy.copy(self.pacf_full)
            dataOut.acfs_error_to_plot[i,:]=numpy.copy(self.sacf_full)
            dataOut.lags_to_plot[i,:]=numpy.copy(self.taup_full)
            dataOut.x_igcej_to_plot[i,:]=numpy.copy(self.x_igcej)
            dataOut.x_ibad_to_plot[i,:]=numpy.copy(self.x_ibad)
            dataOut.y_igcej_to_plot[i,:]=numpy.copy(self.y_igcej)
            dataOut.y_ibad_to_plot[i,:]=numpy.copy(self.y_ibad)

        missing=numpy.nan#-32767
        #print("dataOut.igcej",dataOut.igcej)
        #print("dataOut.ibad",dataOut.ibad)
        for i in range(dataOut.NSHTS,dataOut.NDP):
            for j in range(dataOut.DPL):
                dataOut.acfs_to_save[i,j]=missing
                dataOut.acfs_error_to_save[i,j]=missing
                dataOut.acfs_to_plot[i,j]=missing
                dataOut.acfs_error_to_plot[i,j]=missing
                dataOut.lags_to_plot[i,j]=missing
                dataOut.x_igcej_to_plot[i,j]=missing
                dataOut.x_ibad_to_plot[i,j]=missing
                dataOut.y_igcej_to_plot[i,j]=missing
                dataOut.y_ibad_to_plot[i,j]=missing

        dataOut.acfs_to_save=dataOut.acfs_to_save.transpose()
        dataOut.acfs_error_to_save=dataOut.acfs_error_to_save.transpose()
        
        return dataOut


class CohInt(Operation):

    isConfig = False
    __profIndex = 0
    __byTime = False
    __initime = None
    __lastdatatime = None
    __integrationtime = None
    __buffer = None
    __bufferStride = []
    __dataReady = False
    __profIndexStride = 0
    __dataToPutStride = False
    n = None

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

        #   self.isConfig = False

    def setup(self, n=None, timeInterval=None, stride=None, overlapping=False, byblock=False):
        """
        Set the parameters of the integration class.

        Inputs:

            n               :    Number of coherent integrations
            timeInterval    :    Time of integration. If the parameter "n" is selected this one does not work
            overlapping     :
        """

        self.__initime = None
        self.__lastdatatime = 0
        self.__buffer = None
        self.__dataReady = False
        self.byblock = byblock
        self.stride = stride

        if n == None and timeInterval == None:
            raise ValueError("n or timeInterval should be specified ...")

        if n != None:
            self.n = n
            self.__byTime = False
        else:
            self.__integrationtime = timeInterval #* 60. #if (type(timeInterval)!=integer) -> change this line
            self.n = 9999
            self.__byTime = True

        if overlapping:
            self.__withOverlapping = True
            self.__buffer = None
        else:
            self.__withOverlapping = False
            self.__buffer = 0

        self.__profIndex = 0

    def putData(self, data):

        """
        Add a profile to the __buffer and increase in one the __profileIndex

        """

        if not self.__withOverlapping:
            self.__buffer += data.copy()
            self.__profIndex += 1
            return

        #Overlapping data
        nChannels, nHeis = data.shape
        data = numpy.reshape(data, (1, nChannels, nHeis))

        #If the buffer is empty then it takes the data value
        if self.__buffer is None:
            self.__buffer = data
            self.__profIndex += 1
            return

        #If the buffer length is lower than n then stakcing the data value
        if self.__profIndex < self.n:
            self.__buffer = numpy.vstack((self.__buffer, data))
            self.__profIndex += 1
            return

        #If the buffer length is equal to n then replacing the last buffer value with the data value
        self.__buffer = numpy.roll(self.__buffer, -1, axis=0)
        self.__buffer[self.n-1] = data
        self.__profIndex = self.n
        return


    def pushData(self):
        """
        Return the sum of the last profiles and the profiles used in the sum.

        Affected:

        self.__profileIndex

        """

        if not self.__withOverlapping:
            data = self.__buffer
            n = self.__profIndex

            self.__buffer = 0
            self.__profIndex = 0

            return data, n

        #Integration with Overlapping
        data = numpy.sum(self.__buffer, axis=0)
        # print data
        # raise
        n = self.__profIndex

        return data, n

    def byProfiles(self, data):

        self.__dataReady = False
        avgdata = None
        #         n = None
        # print data
        # raise
        self.putData(data)

        if self.__profIndex == self.n:
            avgdata, n = self.pushData()
            self.__dataReady = True

        return avgdata

    def byTime(self, data, datatime):

        self.__dataReady = False
        avgdata = None
        n = None

        self.putData(data)

        if (datatime - self.__initime) >= self.__integrationtime:
            avgdata, n = self.pushData()
            self.n = n
            self.__dataReady = True

        return avgdata

    def integrateByStride(self, data, datatime):
        # print data
        if self.__profIndex == 0:
            self.__buffer = [[data.copy(), datatime]]
        else:
            self.__buffer.append([data.copy(),datatime])
        self.__profIndex += 1
        self.__dataReady = False

        if self.__profIndex == self.n * self.stride :
            self.__dataToPutStride = True
            self.__profIndexStride = 0
            self.__profIndex = 0
            self.__bufferStride = []
            for i in range(self.stride):
                current = self.__buffer[i::self.stride]
                data = numpy.sum([t[0] for t in current], axis=0)
                avgdatatime = numpy.average([t[1] for t in current])
                # print data
                self.__bufferStride.append((data, avgdatatime))

        if self.__dataToPutStride:
            self.__dataReady = True
            self.__profIndexStride += 1
            if self.__profIndexStride == self.stride:
                self.__dataToPutStride = False
            # print self.__bufferStride[self.__profIndexStride - 1]
            # raise
            return self.__bufferStride[self.__profIndexStride - 1]


        return None, None

    def integrate(self, data, datatime=None):

        if self.__initime == None:
            self.__initime = datatime

        if self.__byTime:
            avgdata = self.byTime(data, datatime)
        else:
            avgdata = self.byProfiles(data)


        self.__lastdatatime = datatime

        if avgdata is None:
            return None, None

        avgdatatime = self.__initime

        deltatime = datatime - self.__lastdatatime

        if not self.__withOverlapping:
            self.__initime = datatime
        else:
            self.__initime += deltatime

        return avgdata, avgdatatime

    def integrateByBlock(self, dataOut):

        times = int(dataOut.data.shape[1]/self.n)
        avgdata = numpy.zeros((dataOut.nChannels, times, dataOut.nHeights), dtype=numpy.complex)

        id_min = 0
        id_max = self.n

        for i in range(times):
            junk = dataOut.data[:,id_min:id_max,:]
            avgdata[:,i,:] = junk.sum(axis=1)
            id_min += self.n
            id_max += self.n

        timeInterval = dataOut.ippSeconds*self.n
        avgdatatime = (times - 1) * timeInterval + dataOut.utctime
        self.__dataReady = True
        return avgdata, avgdatatime

    def run(self, dataOut, n=None, timeInterval=None, stride=None, overlapping=False, byblock=False, **kwargs):

        if not self.isConfig:
            self.setup(n=n, stride=stride, timeInterval=timeInterval, overlapping=overlapping, byblock=byblock, **kwargs)
            self.isConfig = True
        #print("inside")
        if dataOut.flagDataAsBlock:
            """
            Si la data es leida por bloques, dimension = [nChannels, nProfiles, nHeis]
            """

            avgdata, avgdatatime = self.integrateByBlock(dataOut)
            dataOut.nProfiles /= self.n
        else:
            if stride is None:
                avgdata, avgdatatime = self.integrate(dataOut.data, dataOut.utctime)
            else:
                avgdata, avgdatatime = self.integrateByStride(dataOut.data, dataOut.utctime)


        #   dataOut.timeInterval *= n
        dataOut.flagNoData = True

        if self.__dataReady:
            dataOut.data = avgdata
            if not dataOut.flagCohInt:
                dataOut.nCohInt *= self.n
                dataOut.flagCohInt = True
            dataOut.utctime = avgdatatime
            # print avgdata, avgdatatime
            # raise
            #   dataOut.timeInterval = dataOut.ippSeconds * dataOut.nCohInt
            dataOut.flagNoData = False
        return dataOut

class TimesCode(Operation):
    '''
    Written by R. Flores
    '''
    """

    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

    def run(self,dataOut,code):

        #code = numpy.repeat(code, repeats=osamp, axis=1)
        nCodes = numpy.shape(code)[1]
        #nprofcode = dataOut.nProfiles//nCodes
        code = numpy.array(code)
        #print("nHeights",dataOut.nHeights)
        #print("nheicode",nheicode)
        #print("Code.Shape",numpy.shape(code))
        #print("Code",code[0,:])
        nheicode = dataOut.nHeights//nCodes
        res = dataOut.nHeights%nCodes
        '''
        buffer = numpy.zeros((dataOut.nChannels,
                               nprofcode,
                               nCodes,
                               ndataOut.nHeights),
                              dtype='complex')
                              '''
        #exit(1)
        #for ipr in range(dataOut.nProfiles):
        #print(dataOut.nHeights)
        #print(dataOut.data[0,384-2:])
        #print(dataOut.profileIndex)
        #print(dataOut.data[0,:2])
        #print(dataOut.data[0,0:64])
        #print(dataOut.data[0,64:64+64])
        #exit(1)
        for ich in range(dataOut.nChannels):
            for ihe in range(nheicode):
                #print(ihe*nCodes)
                #print((ihe+1)*nCodes)
                #dataOut.data[ich,ipr,ihe*nCodes:nCodes*(ihe+1)]
                #code[ipr,:]
                #print("before",dataOut.data[ich,ipr,ihe*nCodes:nCodes*(ihe+1)])
                #dataOut.data[ich,ipr,ihe*nCodes:nCodes*(ihe+1)] = numpy.prod([dataOut.data[ich,ipr,ihe*nCodes:nCodes*(ihe+1)],code[ipr,:]],axis=0)
                dataOut.data[ich,ihe*nCodes:nCodes*(ihe+1)] = numpy.prod([dataOut.data[ich,ihe*nCodes:nCodes*(ihe+1)],code[dataOut.profileIndex,:]],axis=0)

                #print("after",dataOut.data[ich,ipr,ihe*nCodes:nCodes*(ihe+1)])
                #exit(1)
        #print(dataOut.data[0,:2])
        #exit(1)
        #print(nheicode)
        #print((nheicode)*nCodes)
        #print(((nheicode)*nCodes)+res)
        if res != 0:
            for ich in range(dataOut.nChannels):
                dataOut.data[ich,nheicode*nCodes:] = numpy.prod([dataOut.data[ich,nheicode*nCodes:],code[dataOut.profileIndex,:res]],axis=0)

            #pass
        #print(dataOut.data[0,384-2:])
        #exit(1)
        #dataOut.data = numpy.mean(buffer,axis=1)
        #print(numpy.shape(dataOut.data))
        #print(dataOut.nHeights)
        #dataOut.heightList = dataOut.heightList[0:nheicode]
        #print(dataOut.nHeights)
        #dataOut.nHeights = numpy.shape(dataOut.data)[2]
        #print(numpy.shape(dataOut.data))
        #exit(1)

        return dataOut


class RemoveDcHae(Operation):
    '''
    Written by R. Flores
    '''
    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)
        self.DcCounter = 0

    def run(self, dataOut):

        if self.DcCounter == 0:
            dataOut.DcHae = numpy.zeros((dataOut.data.shape[0],320),dtype='complex')
            #dataOut.DcHae = []
            self.DcCounter = 1

        dataOut.dataaux = numpy.copy(dataOut.data)

        #dataOut.DcHae += dataOut.dataaux[:,1666:1666+320]
        dataOut.DcHae += dataOut.dataaux[:,0:0+320]
        hei = 1666
        hei = 2000
        hei = 1000
        hei = 0
        #dataOut.DcHae = numpy.concatenate([dataOut.DcHae,dataOut.dataaux[0,hei]],axis = None)



        return dataOut


class SSheightProfiles(Operation):

    step          = None
    nsamples      = None
    bufferShape   = None
    profileShape  = None
    sshProfiles   = None
    profileIndex  = None

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)
        self.isConfig = False

    def setup(self,dataOut ,step = None , nsamples = None):

        if step == None and nsamples == None:
            raise ValueError("step or nheights should be specified ...")

        self.step         = step
        self.nsamples     = nsamples
        self.__nChannels  = dataOut.nChannels
        self.__nProfiles  = dataOut.nProfiles
        self.__nHeis      = dataOut.nHeights
        shape             = dataOut.data.shape #nchannels, nprofiles, nsamples

        residue     =  (shape[1] - self.nsamples) % self.step
        if residue != 0:
            print("The residue is %d, step=%d should be multiple of %d to avoid loss of %d samples"%(residue,step,shape[1] - self.nsamples,residue))

        deltaHeight      =  dataOut.heightList[1] - dataOut.heightList[0]
        numberProfile    =  self.nsamples
        numberSamples    =  (shape[1] - self.nsamples)/self.step

        self.bufferShape  = int(shape[0]), int(numberSamples), int(numberProfile)  # nchannels, nsamples , nprofiles
        self.profileShape = int(shape[0]), int(numberProfile), int(numberSamples)  # nchannels, nprofiles, nsamples
        print("buffer shape: ", self.bufferShape)
        self.buffer       = numpy.zeros(self.bufferShape , dtype=numpy.complex)
        self.sshProfiles  = numpy.zeros(self.profileShape, dtype=numpy.complex)

    def run(self, dataOut, step, nsamples, code = None, repeat = None):
        dataOut.flagNoData      = True

        profileIndex            = None
        #print(dataOut.getFreqRange(1)/1000.)
        #exit(1)
        '''
        if dataOut.flagDataAsBlock:
            dataOut.data = numpy.average(dataOut.data,axis=1)
            #print("jee")
            '''
        dataOut.flagDataAsBlock = False
        if not self.isConfig:
            self.setup(dataOut, step=step , nsamples=nsamples)
            self.isConfig = True

        #DC_Hae = numpy.array([0.398+0.588j, -0.926+0.306j, -0.536-0.682j, -0.072+0.53j, 0.368-0.356j, 0.996+0.362j])
        #DC_Hae = numpy.array([ 0.001025  +0.0516375j,   0.03485   +0.20923125j, -0.168     -0.02720625j,
         #-0.1105375 +0.0707125j,  -0.20309375-0.09670625j,  0.189775  +0.02716875j])*(-3.5)

        #DC_Hae = numpy.array([ -32.26  +8.66j,   -32.26  +8.66j])

        #DC_Hae = numpy.array([-2.78500000e-01  -1.39175j,   -6.63237294e+02+210.4268625j])

        #print(dataOut.data[0,13:15])
        #dataOut.data = dataOut.data -  DC_Hae[:,None]
        #print(dataOut.data[0,13:15])
        #exit(1)
        if code is not None:
            code = numpy.array(code)
            code_block = code
            '''
            roll = 0
            code = numpy.roll(code,roll,axis=0)
            code = numpy.reshape(code,(5,100,64))
            block = dataOut.CurrentBlock%5

            day_dif = 0 #day_19_Oct_2021: 3
            code_block = code[block-1-day_dif,:,:]
            '''
            #print("hey i'm here!!!!!!!!!!!!!",code_block)
            if repeat is not None:
                code_block = numpy.repeat(code_block, repeats=repeat, axis=1)
                #print("REPEAT!!!!!!!!!!!!!",code_block)

        for i in range(self.buffer.shape[1]):
            #exit(1)
            if code is not None:
                #self.buffer[:,i] = dataOut.data[:,i*self.step:i*self.step + self.nsamples]*code_block[dataOut.profileIndex,:]
                #print("BEFORE: ", dataOut.data[:,i*self.step:i*self.step + self.nsamples])
                self.buffer[:,i] = dataOut.data[:,i*self.step:i*self.step + self.nsamples]*code_block
                #print("AFTER: ", self.buffer[:,i])

            else:

                self.buffer[:,i]    = dataOut.data[:,i*self.step:i*self.step + self.nsamples]#*code[dataOut.profileIndex,:]

            #self.buffer[:,j,self.__nHeis-j*self.step - self.nheights:self.__nHeis-j*self.step] = numpy.flip(dataOut.data[:,j*self.step:j*self.step + self.nheights])

        for j in range(self.buffer.shape[0]):
            self.sshProfiles[j] = numpy.transpose(self.buffer[j])

        profileIndex  =  self.nsamples
        deltaHeight   =  dataOut.heightList[1] - dataOut.heightList[0]
        ippSeconds    =  (deltaHeight*1.0e-6)/(0.15)
	#print "ippSeconds",ippSeconds
        try:
            if dataOut.concat_m  is not None:
                ippSeconds= ippSeconds/float(dataOut.concat_m)
                #print "Profile concat %d"%dataOut.concat_m
        except:
            pass

        dataOut.data            = self.sshProfiles
        dataOut.flagNoData      = False
        dataOut.heightList      = numpy.arange(self.buffer.shape[1]) *self.step*deltaHeight + dataOut.heightList[0]
        dataOut.nProfiles       = int(dataOut.nProfiles*self.nsamples)

        dataOut.profileIndex    = profileIndex
        dataOut.flagDataAsBlock = True
        dataOut.ippSeconds      = ippSeconds
        dataOut.step            = self.step


        #print("SSH")
        #print(numpy.shape(dataOut.data))
        #exit(1)
        #print(dataOut.data[0,:,150])
        #exit(1)
        #print(dataOut.data[0,:,0]*numpy.conjugate(dataOut.data[0,0,0]))
        #exit(1)

        return dataOut

class removeDCHAE(Operation):

    def run(self, dataOut, minHei, maxHei):

        heights = dataOut.heightList

        inda = numpy.where(heights >= minHei)
        indb = numpy.where(heights <= maxHei)

        minIndex = inda[0][0]
        maxIndex = indb[0][-1]

        dc = numpy.average(dataOut.data[:,minIndex:maxIndex],axis=1)
        #print(dc.shape)
        dataOut.data = dataOut.data - dc[:,None]
        #print(aux.shape)
        #exit(1)

        return dataOut

class Decoder_Original(Operation):

    isConfig = False
    __profIndex = 0

    code = None

    nCode = None
    nBaud = None

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

        self.times = None
        self.osamp = None
    #         self.__setValues = False
        self.isConfig = False
        self.setupReq = False
    def setup(self, code, osamp, dataOut):

        self.__profIndex = 0

        self.code = code

        self.nCode = len(code)
        self.nBaud = len(code[0])

        if (osamp != None) and (osamp >1):
            self.osamp = osamp
            self.code = numpy.repeat(code, repeats=self.osamp, axis=1)
            self.nBaud = self.nBaud*self.osamp

        self.__nChannels = dataOut.nChannels
        self.__nProfiles = dataOut.nProfiles
        self.__nHeis = dataOut.nHeights
        #print("self.__nHeis: ", self.__nHeis)
        #print("self.nBaud: ", self.nBaud)
        #exit(1)
        if self.__nHeis < self.nBaud:
            raise ValueError('Number of heights (%d) should be greater than number of bauds (%d)' %(self.__nHeis, self.nBaud))
        #print("JJE")
        #exit(1)
        #Frequency
        __codeBuffer = numpy.zeros((self.nCode, self.__nHeis), dtype=numpy.complex)

        __codeBuffer[:,0:self.nBaud] = self.code

        self.fft_code = numpy.conj(numpy.fft.fft(__codeBuffer, axis=1))

        if dataOut.flagDataAsBlock:

            self.ndatadec = self.__nHeis #- self.nBaud + 1

            self.datadecTime = numpy.zeros((self.__nChannels, self.__nProfiles, self.ndatadec), dtype=numpy.complex)

        else:

            #Time
            self.ndatadec = self.__nHeis #- self.nBaud + 1


            self.datadecTime = numpy.zeros((self.__nChannels, self.ndatadec), dtype=numpy.complex)

    def __convolutionInFreq(self, data):

        fft_code = self.fft_code[self.__profIndex].reshape(1,-1)

        fft_data = numpy.fft.fft(data, axis=1)

        conv = fft_data*fft_code

        data = numpy.fft.ifft(conv,axis=1)

        return data

    def __convolutionInFreqOpt(self, data):

        raise NotImplementedError

    def __convolutionInTime(self, data):

        code = self.code[self.__profIndex]
        for i in range(self.__nChannels):
            #aux=numpy.correlate(data[i,:], code, mode='full')
            #print(numpy.shape(aux))
            #print(numpy.shape(data[i,:]))
            #print(numpy.shape(code))
            #exit(1)
            self.datadecTime[i,:] = numpy.correlate(data[i,:], code, mode='full')[self.nBaud-1:]

        return self.datadecTime

    def __convolutionByBlockInTime(self, data):

        repetitions = int(self.__nProfiles / self.nCode)
        junk = numpy.lib.stride_tricks.as_strided(self.code, (repetitions, self.code.size), (0, self.code.itemsize))
        junk = junk.flatten()
        code_block = numpy.reshape(junk, (self.nCode*repetitions, self.nBaud))
        profilesList = range(self.__nProfiles)
        #print(numpy.shape(self.datadecTime))
        #print(numpy.shape(data))
        #print(profilesList)
        for i in range(self.__nChannels):
            for j in profilesList:
                #print("data.shape: ", data.shape)
                #print("code_block: ", code_block.shape)
                #print("corr: ", numpy.shape(numpy.correlate(data[i,j,:], code_block[j,:], mode='full')))
                #exit(1)
                self.datadecTime[i,j,:] = numpy.correlate(data[i,j,:], code_block[j,:], mode='full')[self.nBaud-1:]
        return self.datadecTime

    def __convolutionByBlockInFreq(self, data):

        raise NotImplementedError("Decoder by frequency fro Blocks not implemented")


        fft_code = self.fft_code[self.__profIndex].reshape(1,-1)

        fft_data = numpy.fft.fft(data, axis=2)

        conv = fft_data*fft_code

        data = numpy.fft.ifft(conv,axis=2)

        return data


    def run(self, dataOut, code=None, nCode=None, nBaud=None, mode = 0, osamp=None, times=None):

        if dataOut.flagDecodeData:
            print("This data is already decoded, recoding again ...")
        #print("code: ", numpy.shape(code))
        #exit(1)
        if not self.isConfig:

            if code is None:
                if dataOut.code is None:
                    raise ValueError("Code could not be read from %s instance. Enter a value in Code parameter" %dataOut.type)

                code = dataOut.code
            else:
                code = numpy.array(code).reshape(nCode,nBaud)
            self.setup(code, osamp, dataOut)

            self.isConfig = True

            if mode == 3:
                sys.stderr.write("Decoder Warning: mode=%d is not valid, using mode=0\n" %mode)

            if times != None:
                sys.stderr.write("Decoder Warning: Argument 'times' in not used anymore\n")

        if self.code is None:
            print("Fail decoding: Code is not defined.")
            return

        self.__nProfiles = dataOut.nProfiles
        datadec = None

        if mode == 3:
            mode = 0

        if dataOut.flagDataAsBlock:
            """
            Decoding when data have been read as block,
            """

            if mode == 0:
                datadec = self.__convolutionByBlockInTime(dataOut.data)
            if mode == 1:
                datadec = self.__convolutionByBlockInFreq(dataOut.data)
        else:
            """
            Decoding when data have been read profile by profile
            """
            if mode == 0:
                datadec = self.__convolutionInTime(dataOut.data)

            if mode == 1:
                datadec = self.__convolutionInFreq(dataOut.data)

            if mode == 2:
                datadec = self.__convolutionInFreqOpt(dataOut.data)

        if datadec is None:
            raise ValueError("Codification mode selected is not valid: mode=%d. Try selecting 0 or 1" %mode)

        dataOut.code = self.code
        dataOut.nCode = self.nCode
        dataOut.nBaud = self.nBaud

        dataOut.data = datadec#/self.nBaud

        #print("before",dataOut.heightList)
        dataOut.heightList = dataOut.heightList[0:datadec.shape[-1]]
        #print("after",dataOut.heightList)

        dataOut.flagDecodeData = True #asumo q la data esta decodificada

        if self.__profIndex == self.nCode-1:
            self.__profIndex = 0
            return dataOut

        self.__profIndex += 1

        #print("SHAPE",numpy.shape(dataOut.data))

        return dataOut
    #        dataOut.flagDeflipData = True #asumo q la data no esta sin flip

class Decoder(Operation):

    isConfig = False
    __profIndex = 0

    code = None

    nCode = None
    nBaud = None

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

        self.times = None
        self.osamp = None
    #         self.__setValues = False
        self.isConfig = False
        self.setupReq = False
    def setup(self, code, osamp, dataOut):

        self.__profIndex = 0

        self.code = code

        self.nCode = len(code)
        #self.nBaud = len(code[0])
        self.nBaud = int(numpy.shape(code)[-1])

        if (osamp != None) and (osamp >1):
            self.osamp = osamp
            self.code = numpy.repeat(code, repeats=self.osamp, axis=1)
            self.nBaud = self.nBaud*self.osamp

        self.__nChannels = dataOut.nChannels
        self.__nProfiles = dataOut.nProfiles
        self.__nHeis = dataOut.nHeights
        #print("self.__nHeis: ", self.__nHeis)
        #print("self.nBaud: ", self.nBaud)
        #exit(1)
        if self.__nHeis < self.nBaud:
            raise ValueError('Number of heights (%d) should be greater than number of bauds (%d)' %(self.__nHeis, self.nBaud))
        #print("JJE")
        #exit(1)
        '''
        #Frequency
        __codeBuffer = numpy.zeros((self.nCode, self.__nHeis), dtype=numpy.complex)

        __codeBuffer[:,0:self.nBaud] = self.code

        self.fft_code = numpy.conj(numpy.fft.fft(__codeBuffer, axis=1))
        '''
        if dataOut.flagDataAsBlock:

            self.ndatadec = self.__nHeis #- self.nBaud + 1

            self.datadecTime = numpy.zeros((self.__nChannels, self.__nProfiles, self.ndatadec), dtype=numpy.complex)

        else:

            #Time
            self.ndatadec = self.__nHeis #- self.nBaud + 1


            self.datadecTime = numpy.zeros((self.__nChannels, self.ndatadec), dtype=numpy.complex)

    def __convolutionInFreq(self, data):

        fft_code = self.fft_code[self.__profIndex].reshape(1,-1)

        fft_data = numpy.fft.fft(data, axis=1)

        conv = fft_data*fft_code

        data = numpy.fft.ifft(conv,axis=1)

        return data

    def __convolutionInFreqOpt(self, data):

        raise NotImplementedError

    def __convolutionInTime(self, data):

        code = self.code[self.__profIndex]
        for i in range(self.__nChannels):
            #aux=numpy.correlate(data[i,:], code, mode='full')
            #print(numpy.shape(aux))
            #print(numpy.shape(data[i,:]))
            #print(numpy.shape(code))
            #exit(1)
            self.datadecTime[i,:] = numpy.correlate(data[i,:], code, mode='full')[self.nBaud-1:]

        return self.datadecTime

    def __convolutionByBlockInTime(self, data, AutoDecod):

        if not AutoDecod:
            repetitions = int(self.__nProfiles / self.nCode)
            junk = numpy.lib.stride_tricks.as_strided(self.code, (repetitions, self.code.size), (0, self.code.itemsize))
            junk = junk.flatten()
            code_block = numpy.reshape(junk, (self.nCode*repetitions, self.nBaud))
        else:
            code_block = self.code
        profilesList = range(self.__nProfiles)
        #print(numpy.shape(self.datadecTime))
        #print(numpy.shape(data))
        #print(profilesList)
        for i in range(self.__nChannels):
            for j in profilesList:
                #print("data.shape: ", data.shape)
                #print("code_block: ", code_block.shape)
                #print("corr: ", numpy.shape(numpy.correlate(data[i,j,:], code_block[j,:], mode='full')))
                #exit(1)
                if not AutoDecod:
                    self.datadecTime[i,j,:] = numpy.correlate(data[i,j,:], code_block[j,:], mode='full')[self.nBaud-1:]
                else:
                    if i%2:
                        self.datadecTime[i,j,:] = numpy.correlate(data[i,j,:], code_block[0,j,:], mode='full')[self.nBaud-1:]
                    else:
                        self.datadecTime[i,j,:] = numpy.correlate(data[i,j,:], code_block[1,j,:], mode='full')[self.nBaud-1:]

        return self.datadecTime

    def __convolutionByBlockInFreq(self, data):

        raise NotImplementedError("Decoder by frequency fro Blocks not implemented")


        fft_code = self.fft_code[self.__profIndex].reshape(1,-1)

        fft_data = numpy.fft.fft(data, axis=2)

        conv = fft_data*fft_code

        data = numpy.fft.ifft(conv,axis=2)

        return data


    def run(self, dataOut, code=None, nCode=None, nBaud=None, mode = 0, osamp=None, times=None, AutoDecod = 0):
        if dataOut.flagDecodeData:
            print("This data is already decoded, recoding again ...")
        #print("code: ", code, numpy.shape(code))

        if not self.isConfig:

            if code is None and not AutoDecod:
                if dataOut.code is None:
                    raise ValueError("Code could not be read from %s instance. Enter a value in Code parameter" %dataOut.type)

                code = dataOut.code
            else:
                if not AutoDecod:
                    code = numpy.array(code).reshape(nCode,nBaud)
                else:
                    po = 0
                    print("***********AutoDecod***********")
                    code = dataOut.data[:2,:,po:po+64]
                    #print("AutoDecod Shape: ", numpy.shape(code))
                    #exit(1)

            self.setup(code, osamp, dataOut)
            self.isConfig = True

            if mode == 3:
                sys.stderr.write("Decoder Warning: mode=%d is not valid, using mode=0\n" %mode)

            if times != None:
                sys.stderr.write("Decoder Warning: Argument 'times' in not used anymore\n")

        if self.code is None:
            print("Fail decoding: Code is not defined.")
            return

        self.__nProfiles = dataOut.nProfiles
        datadec = None

        if mode == 3:
            mode = 0

        if dataOut.flagDataAsBlock:
            """
            Decoding when data have been read as block,
            """

            if mode == 0:
                datadec = self.__convolutionByBlockInTime(dataOut.data, AutoDecod)
            if mode == 1:
                datadec = self.__convolutionByBlockInFreq(dataOut.data)
        else:
            """
            Decoding when data have been read profile by profile
            """
            if mode == 0:
                datadec = self.__convolutionInTime(dataOut.data)

            if mode == 1:
                datadec = self.__convolutionInFreq(dataOut.data)

            if mode == 2:
                datadec = self.__convolutionInFreqOpt(dataOut.data)

        if datadec is None:
            raise ValueError("Codification mode selected is not valid: mode=%d. Try selecting 0 or 1" %mode)

        dataOut.code = self.code
        dataOut.nCode = self.nCode
        dataOut.nBaud = self.nBaud

        dataOut.data = datadec#/self.nBaud

        #print("before",dataOut.heightList)
        dataOut.heightList = dataOut.heightList[0:datadec.shape[-1]]
        #print("after",dataOut.heightList)

        dataOut.flagDecodeData = True #asumo q la data esta decodificada

        if self.__profIndex == self.nCode-1:
            self.__profIndex = 0
            return dataOut

        self.__profIndex += 1

        return dataOut


class DecoderRoll(Operation):

    isConfig = False
    __profIndex = 0

    code = None

    nCode = None
    nBaud = None

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

        self.times = None
        self.osamp = None
    #         self.__setValues = False
        self.isConfig = False
        self.setupReq = False
    def setup(self, code, osamp, dataOut):

        self.__profIndex = 0


        self.code = code

        self.nCode = len(code)
        self.nBaud = len(code[0])

        if (osamp != None) and (osamp >1):
            self.osamp = osamp
            self.code = numpy.repeat(code, repeats=self.osamp, axis=1)
            self.nBaud = self.nBaud*self.osamp

        self.__nChannels = dataOut.nChannels
        self.__nProfiles = dataOut.nProfiles
        self.__nHeis = dataOut.nHeights

        if self.__nHeis < self.nBaud:
            raise ValueError('Number of heights (%d) should be greater than number of bauds (%d)' %(self.__nHeis, self.nBaud))

        #Frequency
        __codeBuffer = numpy.zeros((self.nCode, self.__nHeis), dtype=numpy.complex)

        __codeBuffer[:,0:self.nBaud] = self.code

        self.fft_code = numpy.conj(numpy.fft.fft(__codeBuffer, axis=1))

        if dataOut.flagDataAsBlock:

            self.ndatadec = self.__nHeis #- self.nBaud + 1

            self.datadecTime = numpy.zeros((self.__nChannels, self.__nProfiles, self.ndatadec), dtype=numpy.complex)

        else:

            #Time
            self.ndatadec = self.__nHeis #- self.nBaud + 1


            self.datadecTime = numpy.zeros((self.__nChannels, self.ndatadec), dtype=numpy.complex)

    def __convolutionInFreq(self, data):

        fft_code = self.fft_code[self.__profIndex].reshape(1,-1)

        fft_data = numpy.fft.fft(data, axis=1)

        conv = fft_data*fft_code

        data = numpy.fft.ifft(conv,axis=1)

        return data

    def __convolutionInFreqOpt(self, data):

        raise NotImplementedError

    def __convolutionInTime(self, data):

        code = self.code[self.__profIndex]
        #print("code",code[0,0])
        for i in range(self.__nChannels):
            #aux=numpy.correlate(data[i,:], code, mode='full')
            #print(numpy.shape(aux))
            #print(numpy.shape(data[i,:]))
            #print(numpy.shape(code))
            #exit(1)
            self.datadecTime[i,:] = numpy.correlate(data[i,:], code, mode='full')[self.nBaud-1:]

        return self.datadecTime

    def __convolutionByBlockInTime(self, data):

        repetitions = int(self.__nProfiles / self.nCode)
        junk = numpy.lib.stride_tricks.as_strided(self.code, (repetitions, self.code.size), (0, self.code.itemsize))
        junk = junk.flatten()
        code_block = numpy.reshape(junk, (self.nCode*repetitions, self.nBaud))
        profilesList = range(self.__nProfiles)
        #print(numpy.shape(self.datadecTime))
        #print(numpy.shape(data))
        for i in range(self.__nChannels):
            for j in profilesList:
                self.datadecTime[i,j,:] = numpy.correlate(data[i,j,:], code_block[j,:], mode='full')[self.nBaud-1:]
        return self.datadecTime

    def __convolutionByBlockInFreq(self, data):

        raise NotImplementedError("Decoder by frequency fro Blocks not implemented")


        fft_code = self.fft_code[self.__profIndex].reshape(1,-1)

        fft_data = numpy.fft.fft(data, axis=2)

        conv = fft_data*fft_code

        data = numpy.fft.ifft(conv,axis=2)

        return data


    def run(self, dataOut, code=None, nCode=None, nBaud=None, mode = 0, osamp=None, times=None):

        if dataOut.flagDecodeData:
            print("This data is already decoded, recoding again ...")


        roll = 0

        if self.isConfig:
            code = numpy.array(code)

            code = numpy.roll(code,roll,axis=0)
            code = numpy.reshape(code,(5,100,64))
            block = dataOut.CurrentBlock%5
            #code = code[block-1,:,:] #NormalizeDPPower
            code = code[block-1-1,:,:] #Next Day
            self.code = numpy.repeat(code, repeats=self.osamp, axis=1)


        if not self.isConfig:

            if code is None:
                if dataOut.code is None:
                    raise ValueError("Code could not be read from %s instance. Enter a value in Code parameter" %dataOut.type)

                code = dataOut.code
            else:
                code = numpy.array(code)

                #roll = 29
                code = numpy.roll(code,roll,axis=0)
                code = numpy.reshape(code,(5,100,64))
                block = dataOut.CurrentBlock%5
                code = code[block-1-1,:,:]
                #print(code.shape())
                #exit(1)

                code = numpy.array(code).reshape(nCode,nBaud)
            self.setup(code, osamp, dataOut)

            self.isConfig = True

            if mode == 3:
                sys.stderr.write("Decoder Warning: mode=%d is not valid, using mode=0\n" %mode)

            if times != None:
                sys.stderr.write("Decoder Warning: Argument 'times' in not used anymore\n")

        if self.code is None:
            print("Fail decoding: Code is not defined.")
            return

        self.__nProfiles = dataOut.nProfiles
        datadec = None

        if mode == 3:
            mode = 0

        if dataOut.flagDataAsBlock:
            """
            Decoding when data have been read as block,
            """

            if mode == 0:
                datadec = self.__convolutionByBlockInTime(dataOut.data)
            if mode == 1:
                datadec = self.__convolutionByBlockInFreq(dataOut.data)
        else:
            """
            Decoding when data have been read profile by profile
            """
            if mode == 0:
                datadec = self.__convolutionInTime(dataOut.data)

            if mode == 1:
                datadec = self.__convolutionInFreq(dataOut.data)

            if mode == 2:
                datadec = self.__convolutionInFreqOpt(dataOut.data)

        if datadec is None:
            raise ValueError("Codification mode selected is not valid: mode=%d. Try selecting 0 or 1" %mode)

        dataOut.code = self.code
        dataOut.nCode = self.nCode
        dataOut.nBaud = self.nBaud

        dataOut.data = datadec
        #print("before",dataOut.heightList)
        dataOut.heightList = dataOut.heightList[0:datadec.shape[-1]]
        #print("after",dataOut.heightList)

        dataOut.flagDecodeData = True #asumo q la data esta decodificada

        if self.__profIndex == self.nCode-1:
            self.__profIndex = 0
            return dataOut

        self.__profIndex += 1

        #print("SHAPE",numpy.shape(dataOut.data))

        return dataOut


class ProfileConcat(Operation):

    isConfig = False
    buffer = None

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)
        self.profileIndex = 0

    def reset(self):
        self.buffer = numpy.zeros_like(self.buffer)
        self.start_index = 0
        self.times = 1

    def setup(self, data, m, n=1):
        self.buffer = numpy.zeros((data.shape[0],data.shape[1]*m),dtype=type(data[0,0]))
        self.nHeights = data.shape[1]#.nHeights
        self.start_index = 0
        self.times = 1

    def concat(self, data):

        self.buffer[:,self.start_index:self.nHeights*self.times] = data.copy()
        self.start_index = self.start_index + self.nHeights

    def run(self, dataOut, m):
        dataOut.flagNoData = True

        if not self.isConfig:
            self.setup(dataOut.data, m, 1)
            self.isConfig = True

        if dataOut.flagDataAsBlock:
            raise ValueError("ProfileConcat can only be used when voltage have been read profile by profile, getBlock = False")

        else:
            self.concat(dataOut.data)
            self.times += 1
            if self.times > m:
                dataOut.data = self.buffer
                self.reset()
                dataOut.flagNoData = False
                # se deben actualizar mas propiedades del header y del objeto dataOut, por ejemplo, las alturas
                deltaHeight = dataOut.heightList[1] - dataOut.heightList[0]
                xf = dataOut.heightList[0] + dataOut.nHeights * deltaHeight * m
                dataOut.heightList = numpy.arange(dataOut.heightList[0], xf, deltaHeight)
                dataOut.ippSeconds *= m
        return dataOut

class ProfileSelector(Operation):

    profileIndex = None
    # Tamanho total de los perfiles
    nProfiles = None

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)
        self.profileIndex = 0

    def incProfileIndex(self):

        self.profileIndex += 1
        
        if self.profileIndex >= self.nProfiles:
            self.profileIndex = 0

    def isThisProfileInRange(self, profileIndex, minIndex, maxIndex):

        if profileIndex < minIndex:
            return False

        if profileIndex > maxIndex:
            return False

        return True

    def isThisProfileInList(self, profileIndex, profileList):

        if profileIndex not in profileList:
            return False

        return True

    def run(self, dataOut, profileList=None, profileRangeList=None, beam=None, byblock=False, rangeList = None, nProfiles=None):

        """
        ProfileSelector:

        Inputs:
            profileList        :    Index of profiles selected. Example: profileList = (0,1,2,7,8)

            profileRangeList    :    Minimum and maximum profile indexes. Example: profileRangeList = (4, 30)

            rangeList            :    List of profile ranges. Example: rangeList = ((4, 30), (32, 64), (128, 256))

        """
        #print("nProfiles,", self.nProfiles)
        #print("dataOut.flagDataAsBlock", dataOut.flagDataAsBlock)
        #print("dataOut.data", numpy.shape(dataOut.data))
        if rangeList is not None:
            if type(rangeList[0]) not in (tuple, list):
                rangeList = [rangeList]

        dataOut.flagNoData = True

        if dataOut.flagDataAsBlock:
            """
            data dimension  = [nChannels, nProfiles, nHeis]
            """
            if profileList != None:
                dataOut.data = dataOut.data[:,profileList,:]

            if profileRangeList != None:
                minIndex = profileRangeList[0]
                maxIndex = profileRangeList[1]
                profileList = list(range(minIndex, maxIndex+1))

                dataOut.data = dataOut.data[:,minIndex:maxIndex+1,:]

            if rangeList != None:

                profileList = []

                for thisRange in rangeList:
                    minIndex = thisRange[0]
                    maxIndex = thisRange[1]

                    profileList.extend(list(range(minIndex, maxIndex+1)))

                dataOut.data = dataOut.data[:,profileList,:]

            dataOut.nProfiles = len(profileList)
            dataOut.profileIndex = dataOut.nProfiles - 1
            dataOut.flagNoData = False
            #print("Shape after prof select: ", numpy.shape(dataOut.data))
            #print(dataOut.heightList)
            #exit(1)
            '''
            po = 5
            import matplotlib.pyplot as plt
            this_data = dataOut.data[0,0,:]
            plt.plot(this_data*numpy.conjugate(this_data),marker='*',linestyle='-')
            plt.axvline(po)
            plt.axvline(po+64-1)
            plt.grid()
            plt.show()
            '''
            return dataOut

        """
        data dimension  = [nChannels, nHeis]
        """

        if profileList != None:

            if self.isThisProfileInList(dataOut.profileIndex, profileList):

                self.nProfiles = len(profileList)
                dataOut.nProfiles = self.nProfiles
                dataOut.profileIndex = self.profileIndex
                dataOut.flagNoData = False

                self.incProfileIndex()
            return dataOut

        if profileRangeList != None:

            minIndex = profileRangeList[0]
            maxIndex = profileRangeList[1]

            if self.isThisProfileInRange(dataOut.profileIndex, minIndex, maxIndex):

                self.nProfiles = maxIndex - minIndex + 1
                dataOut.nProfiles = self.nProfiles
                dataOut.profileIndex = self.profileIndex
                dataOut.flagNoData = False

                self.incProfileIndex()
            return dataOut

        if rangeList != None:

            nProfiles = 0

            for thisRange in rangeList:
                minIndex = thisRange[0]
                maxIndex = thisRange[1]

                nProfiles += maxIndex - minIndex + 1

            for thisRange in rangeList:

                minIndex = thisRange[0]
                maxIndex = thisRange[1]

                if self.isThisProfileInRange(dataOut.profileIndex, minIndex, maxIndex):

                    self.nProfiles = nProfiles
                    dataOut.nProfiles = self.nProfiles
                    dataOut.profileIndex = self.profileIndex
                    dataOut.flagNoData = False

                    self.incProfileIndex()

                    break

            return dataOut


        if beam != None: #beam is only for AMISR data
            if self.isThisProfileInList(dataOut.profileIndex, dataOut.beamRangeDict[beam]):
                dataOut.flagNoData = False
                dataOut.profileIndex = self.profileIndex

                self.incProfileIndex()

            return dataOut

        raise ValueError("ProfileSelector needs profileList, profileRangeList or rangeList parameter")

        #return False
        return dataOut

class Reshaper(Operation):

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

        self.__buffer = None
        self.__nitems = 0

    def __appendProfile(self, dataOut, nTxs):

        if self.__buffer is None:
            shape = (dataOut.nChannels, int(dataOut.nHeights/nTxs) )
            self.__buffer = numpy.empty(shape, dtype = dataOut.data.dtype)

        ini = dataOut.nHeights * self.__nitems
        end = ini + dataOut.nHeights

        self.__buffer[:, ini:end] = dataOut.data

        self.__nitems += 1

        return int(self.__nitems*nTxs)

    def __getBuffer(self):

        if self.__nitems == int(1./self.__nTxs):

            self.__nitems = 0

            return self.__buffer.copy()

        return None

    def __checkInputs(self, dataOut, shape, nTxs):

        if shape is None and nTxs is None:
            raise ValueError("Reshaper: shape of factor should be defined")

        if nTxs:
            if nTxs < 0:
                raise ValueError("nTxs should be greater than 0")

            if nTxs < 1 and dataOut.nProfiles % (1./nTxs) != 0:
                raise ValueError("nProfiles= %d is not divisibled by (1./nTxs) = %f" %(dataOut.nProfiles, (1./nTxs)))

            shape = [dataOut.nChannels, dataOut.nProfiles*nTxs, dataOut.nHeights/nTxs]

            return shape, nTxs

        if len(shape) != 2 and len(shape) !=  3:
            raise ValueError("shape dimension should be equal to 2 or 3. shape = (nProfiles, nHeis) or (nChannels, nProfiles, nHeis). Actually shape = (%d, %d, %d)" %(dataOut.nChannels, dataOut.nProfiles, dataOut.nHeights))

        if len(shape) == 2:
            shape_tuple = [dataOut.nChannels]
            shape_tuple.extend(shape)
        else:
            shape_tuple = list(shape)
        #print(shape_tuple)
        #print(dataOut.nProfiles)
        nTxs = 1.0*shape_tuple[1]/dataOut.nProfiles

        return shape_tuple, nTxs

    def run(self, dataOut, shape=None, nTxs=None):
        #print(numpy.shape(dataOut.data))
        #exit(1)
        shape_tuple, self.__nTxs = self.__checkInputs(dataOut, shape, nTxs)
        #print(self.__nTxs)
        #print(dataOut.flagDataAsBlock)
        dataOut.flagNoData = True
        profileIndex = None

        if dataOut.flagDataAsBlock:

            dataOut.data = numpy.reshape(dataOut.data, shape_tuple)
            dataOut.flagNoData = False

            profileIndex = int(dataOut.nProfiles*self.__nTxs) - 1

        else:


            if self.__nTxs < 1:

                self.__appendProfile(dataOut, self.__nTxs)
                new_data = self.__getBuffer()

                if new_data is not None:
                    dataOut.data = new_data
                    dataOut.flagNoData = False

                    profileIndex = dataOut.profileIndex*nTxs

            else:
                raise ValueError("nTxs should be greater than 0 and lower than 1, or use VoltageReader(..., getblock=True)")

        deltaHeight = dataOut.heightList[1] - dataOut.heightList[0]

        dataOut.heightList = numpy.arange(dataOut.nHeights/self.__nTxs) * deltaHeight + dataOut.heightList[0]

        dataOut.nProfiles = int(dataOut.nProfiles*self.__nTxs)

        dataOut.profileIndex = profileIndex

        dataOut.ippSeconds /= self.__nTxs
        #print(numpy.shape(dataOut.data))
        #exit(1)

        return dataOut

class SplitProfiles(Operation):

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

    def run(self, dataOut, n):

        dataOut.flagNoData = True
        profileIndex = None

        if dataOut.flagDataAsBlock:

            #nchannels, nprofiles, nsamples
            shape = dataOut.data.shape

            if shape[2] % n != 0:
                raise ValueError("Could not split the data, n=%d has to be multiple of %d" %(n, shape[2]))

            new_shape = shape[0], shape[1]*n, int(shape[2]/n)

            dataOut.data = numpy.reshape(dataOut.data, new_shape)
            dataOut.flagNoData = False

            profileIndex = int(dataOut.nProfiles/n) - 1

        else:

            raise ValueError("Could not split the data when is read Profile by Profile. Use VoltageReader(..., getblock=True)")

        deltaHeight = dataOut.heightList[1] - dataOut.heightList[0]

        dataOut.heightList = numpy.arange(dataOut.nHeights/n) * deltaHeight + dataOut.heightList[0]

        dataOut.nProfiles = int(dataOut.nProfiles*n)

        dataOut.profileIndex = profileIndex

        dataOut.ippSeconds /= n

        return dataOut

class CombineProfiles(Operation):
    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

        self.__remData = None
        self.__profileIndex = 0

    def run(self, dataOut, n):

        dataOut.flagNoData = True
        profileIndex = None

        if dataOut.flagDataAsBlock:

            #nchannels, nprofiles, nsamples
            shape = dataOut.data.shape
            new_shape = shape[0], shape[1]/n, shape[2]*n

            if shape[1] % n != 0:
                raise ValueError("Could not split the data, n=%d has to be multiple of %d" %(n, shape[1]))

            dataOut.data = numpy.reshape(dataOut.data, new_shape)
            dataOut.flagNoData = False

            profileIndex = int(dataOut.nProfiles*n) - 1

        else:

            #nchannels, nsamples
            if self.__remData is None:
                newData = dataOut.data
            else:
                newData = numpy.concatenate((self.__remData, dataOut.data), axis=1)

            self.__profileIndex += 1

            if self.__profileIndex < n:
                self.__remData = newData
                #continue
                return

            self.__profileIndex = 0
            self.__remData = None

            dataOut.data = newData
            dataOut.flagNoData = False

            profileIndex = dataOut.profileIndex/n


        deltaHeight = dataOut.heightList[1] - dataOut.heightList[0]

        dataOut.heightList = numpy.arange(dataOut.nHeights*n) * deltaHeight + dataOut.heightList[0]

        dataOut.nProfiles = int(dataOut.nProfiles/n)

        dataOut.profileIndex = profileIndex

        dataOut.ippSeconds *= n

        return dataOut
# import collections
# from scipy.stats import mode
#
# class Synchronize(Operation):
#
#     isConfig = False
#     __profIndex = 0
#
#     def __init__(self, **kwargs):
#
#         Operation.__init__(self, **kwargs)
# #         self.isConfig = False
#         self.__powBuffer = None
#         self.__startIndex = 0
#         self.__pulseFound = False
#
#     def __findTxPulse(self, dataOut, channel=0, pulse_with = None):
#
#         #Read data
#
#         powerdB = dataOut.getPower(channel = channel)
#         noisedB = dataOut.getNoise(channel = channel)[0]
#
#         self.__powBuffer.extend(powerdB.flatten())
#
#         dataArray = numpy.array(self.__powBuffer)
#
#         filteredPower = numpy.correlate(dataArray, dataArray[0:self.__nSamples], "same")
#
#         maxValue = numpy.nanmax(filteredPower)
#
#         if maxValue < noisedB + 10:
#             #No se encuentra ningun pulso de transmision
#             return None
#
#         maxValuesIndex = numpy.where(filteredPower > maxValue - 0.1*abs(maxValue))[0]
#
#         if len(maxValuesIndex) < 2:
#             #Solo se encontro un solo pulso de transmision de un baudio, esperando por el siguiente TX
#             return None
#
#         phasedMaxValuesIndex = maxValuesIndex - self.__nSamples
#
#         #Seleccionar solo valores con un espaciamiento de nSamples
#         pulseIndex = numpy.intersect1d(maxValuesIndex, phasedMaxValuesIndex)
#
#         if len(pulseIndex) < 2:
#             #Solo se encontro un pulso de transmision con ancho mayor a 1
#             return None
#
#         spacing = pulseIndex[1:] - pulseIndex[:-1]
#
#         #remover senales que se distancien menos de 10 unidades o muestras
#         #(No deberian existir IPP menor a 10 unidades)
#
#         realIndex = numpy.where(spacing > 10 )[0]
#
#         if len(realIndex) < 2:
#             #Solo se encontro un pulso de transmision con ancho mayor a 1
#             return None
#
#         #Eliminar pulsos anchos (deja solo la diferencia entre IPPs)
#         realPulseIndex = pulseIndex[realIndex]
#
#         period = mode(realPulseIndex[1:] - realPulseIndex[:-1])[0][0]
#
#         print "IPP = %d samples" %period
#
#         self.__newNSamples = dataOut.nHeights #int(period)
#         self.__startIndex = int(realPulseIndex[0])
#
#         return 1
#
#
#     def setup(self, nSamples, nChannels, buffer_size = 4):
#
#         self.__powBuffer = collections.deque(numpy.zeros( buffer_size*nSamples,dtype=numpy.float),
#                                           maxlen = buffer_size*nSamples)
#
#         bufferList = []
#
#         for i in range(nChannels):
#             bufferByChannel = collections.deque(numpy.zeros( buffer_size*nSamples, dtype=numpy.complex) +  numpy.NAN,
#                                           maxlen = buffer_size*nSamples)
#
#             bufferList.append(bufferByChannel)
#
#         self.__nSamples = nSamples
#         self.__nChannels = nChannels
#         self.__bufferList = bufferList
#
#     def run(self, dataOut, channel = 0):
#
#         if not self.isConfig:
#             nSamples = dataOut.nHeights
#             nChannels = dataOut.nChannels
#             self.setup(nSamples, nChannels)
#             self.isConfig = True
#
#         #Append new data to internal buffer
#         for thisChannel in range(self.__nChannels):
#             bufferByChannel = self.__bufferList[thisChannel]
#             bufferByChannel.extend(dataOut.data[thisChannel])
#
#         if self.__pulseFound:
#             self.__startIndex -= self.__nSamples
#
#         #Finding Tx Pulse
#         if not self.__pulseFound:
#             indexFound = self.__findTxPulse(dataOut, channel)
#
#             if indexFound == None:
#                 dataOut.flagNoData = True
#                 return
#
#             self.__arrayBuffer = numpy.zeros((self.__nChannels, self.__newNSamples), dtype = numpy.complex)
#             self.__pulseFound = True
#             self.__startIndex = indexFound
#
#         #If pulse was found ...
#         for thisChannel in range(self.__nChannels):
#             bufferByChannel = self.__bufferList[thisChannel]
#             #print self.__startIndex
#             x = numpy.array(bufferByChannel)
#             self.__arrayBuffer[thisChannel] = x[self.__startIndex:self.__startIndex+self.__newNSamples]
#
#         deltaHeight = dataOut.heightList[1] - dataOut.heightList[0]
#         dataOut.heightList = numpy.arange(self.__newNSamples)*deltaHeight
# #             dataOut.ippSeconds = (self.__newNSamples / deltaHeight)/1e6
#
#         dataOut.data = self.__arrayBuffer
#
#         self.__startIndex += self.__newNSamples
#
#         return







##############################LONG PULSE##############################



class CrossProdHybrid(CrossProdDP):
    """Operation to calculate cross products of the Hybrid Experiment.

    Parameters:
    -----------
    NLAG : int
        Number of lags for Long Pulse.
    NRANGE : int
        Number of samples (heights) for Long Pulse.
    NCAL : int
        .*
    DPL : int
        Number of lags for Double Pulse.
    NDN : int
        .*
    NDT : int
        Number of heights for Double Pulse.*
    NDP : int
        Number of heights for Double Pulse.*
    NSCAN : int
        Number of profiles when the transmitter is on.
    lagind : intlist
        .*
    lagfirst : intlist
        .*
    NAVG : int
        Number of blocks to be "averaged".
    nkill : int
        Number of blocks not to be considered when averaging.

    Example
    --------

    op = proc_unit.addOperation(name='CrossProdHybrid', optype='other')
    op.addParameter(name='NLAG', value='16', format='int')
    op.addParameter(name='NRANGE', value='200', format='int')
    op.addParameter(name='NCAL', value='0', format='int')
    op.addParameter(name='DPL', value='11', format='int')
    op.addParameter(name='NDN', value='0', format='int')
    op.addParameter(name='NDT', value='67', format='int')
    op.addParameter(name='NDP', value='67', format='int')
    op.addParameter(name='NSCAN', value='128', format='int')
    op.addParameter(name='lagind', value='(0,1,2,3,4,5,6,7,0,3,4,5,6,8,9,10)', format='intlist')
    op.addParameter(name='lagfirst', value='(1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1)', format='intlist')
    op.addParameter(name='NAVG', value='16', format='int')
    op.addParameter(name='nkill', value='6', format='int')

    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)
        self.bcounter=0
        self.aux=1
        self.aux_cross_lp=1
        self.lag_products_LP_median_estimates_aux=1

    def get_products_cabxys_HP(self,dataOut):

        if self.aux==1:
            self.set_header_output(dataOut)
            self.aux=0

        self.cax=numpy.zeros((dataOut.NDP,dataOut.DPL,2))# hp:67x11x2  dp: 66x11x2
        self.cay=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.cbx=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.cby=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.cax2=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.cay2=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.cbx2=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.cby2=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.caxbx=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.caxby=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.caybx=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.cayby=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.caxay=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        self.cbxby=numpy.zeros((dataOut.NDP,dataOut.DPL,2))
        for i in range(2):   # flipped and unflipped
            for j in range(dataOut.NDP): # loop over true ranges # 67
                for k in range(int(dataOut.NSCAN)): # 128

                    n=dataOut.lagind[k%dataOut.NLAG] # 128=16x8

                    ax=dataOut.data[0,k,dataOut.NRANGE+dataOut.NCAL+j+i*dataOut.NDT].real#-dataOut.dc.real[0]
                    ay=dataOut.data[0,k,dataOut.NRANGE+dataOut.NCAL+j+i*dataOut.NDT].imag#-dataOut.dc.imag[0]

                    if dataOut.NRANGE+dataOut.NCAL+j+i*dataOut.NDT+2*n<dataOut.read_samples:

                        bx=dataOut.data[1,k,dataOut.NRANGE+dataOut.NCAL+j+i*dataOut.NDT+2*n].real#-dataOut.dc.real[1]
                        by=dataOut.data[1,k,dataOut.NRANGE+dataOut.NCAL+j+i*dataOut.NDT+2*n].imag#-dataOut.dc.imag[1]
                        #print(bx)
                        #print(by)
                        #exit(1)
                    else:
                        #print(i,j,n)
                        #exit(1)

                        if k+1<int(dataOut.NSCAN):
                            bx=dataOut.data[1,k+1,(dataOut.NRANGE+dataOut.NCAL+j+i*dataOut.NDT+2*n)%dataOut.NDP].real
                            by=dataOut.data[1,k+1,(dataOut.NRANGE+dataOut.NCAL+j+i*dataOut.NDT+2*n)%dataOut.NDP].imag
                            #print(bx)
                            #print(by)
                            #exit(1)
                        if k+1==int(dataOut.NSCAN):## ESTO ES UN PARCHE PUES NO SE TIENE EL SIGUIENTE BLOQUE
                            bx=dataOut.data[1,k,(dataOut.NRANGE+dataOut.NCAL+j+i*dataOut.NDT+2*n)%dataOut.NDP].real
                            by=dataOut.data[1,k,(dataOut.NRANGE+dataOut.NCAL+j+i*dataOut.NDT+2*n)%dataOut.NDP].imag
                    #if n == 7 and j == 65:
                        #print(k)
                        #print(bx)
                        #print(by)
                            #exit(1)
                    if(k<dataOut.NLAG and dataOut.lagfirst[k%dataOut.NLAG]==1):# if(k<16 && lagfirst[k%16]==1)
                        self.cax[j][n][i]=ax
                        self.cay[j][n][i]=ay
                        self.cbx[j][n][i]=bx
                        self.cby[j][n][i]=by
                        self.cax2[j][n][i]=ax*ax
                        self.cay2[j][n][i]=ay*ay
                        self.cbx2[j][n][i]=bx*bx
                        self.cby2[j][n][i]=by*by
                        self.caxbx[j][n][i]=ax*bx
                        self.caxby[j][n][i]=ax*by
                        self.caybx[j][n][i]=ay*bx
                        self.cayby[j][n][i]=ay*by
                        self.caxay[j][n][i]=ax*ay
                        self.cbxby[j][n][i]=bx*by
                    else:
                        self.cax[j][n][i]+=ax
                        self.cay[j][n][i]+=ay
                        self.cbx[j][n][i]+=bx
                        self.cby[j][n][i]+=by
                        self.cax2[j][n][i]+=ax*ax
                        self.cay2[j][n][i]+=ay*ay
                        self.cbx2[j][n][i]+=bx*bx
                        self.cby2[j][n][i]+=by*by
                        self.caxbx[j][n][i]+=ax*bx
                        self.caxby[j][n][i]+=ax*by
                        self.caybx[j][n][i]+=ay*bx
                        self.cayby[j][n][i]+=ay*by
                        self.caxay[j][n][i]+=ax*ay
                        self.cbxby[j][n][i]+=bx*by
                    '''
                    if j == 20 and n == 10:
                        print("i: ",i)
                        print(ax+ay*1.j)
                        print("b",bx+by*1.j)
                        '''
                    #if n ==7:
                        #print()


        # FindMe
        pa1 = 20
        pa2 = 0
        #for pa1 in range(67):
        '''
        print(self.cax[pa1,pa2,0]+self.cax[pa1,pa2,1]+self.cay[pa1,pa2,0]*1.j+self.cay[pa1,pa2,1]*1.j)
        print("b",self.cbx[pa1,pa2,0]+self.cbx[pa1,pa2,1]+self.cby[pa1,pa2,0]*1.j+self.cby[pa1,pa2,1]*1.j)
        '''
        '''
        print(self.caxbx[pa1,pa2,0]+self.caxbx[pa1,pa2,1]+self.cayby[pa1,pa2,0]+self.cayby[pa1,pa2,1])
        print("c",self.caybx[pa1,pa2,0]+self.caybx[pa1,pa2,1]-self.caxby[pa1,pa2,0]-self.caxby[pa1,pa2,1])

        exit(1)
        '''


    def lag_products_LP(self,dataOut):


        buffer=dataOut.data
        if self.aux_cross_lp==1:

            #self.dataOut.nptsfft2=150
            self.cnorm=float((dataOut.nProfiles-dataOut.NSCAN)/dataOut.NSCAN)
            self.lagp0=numpy.zeros((dataOut.NLAG,dataOut.NRANGE,dataOut.NAVG),'complex128')
            ww=numpy.zeros((dataOut.NLAG,dataOut.NRANGE,dataOut.NSCAN,dataOut.NAVG),'complex128')
            self.lagp1=numpy.zeros((dataOut.NLAG,dataOut.NRANGE,dataOut.NAVG),'complex128')
            self.lagp2=numpy.zeros((dataOut.NLAG,dataOut.NRANGE,dataOut.NAVG),'complex128')
            self.lagp3=numpy.zeros((dataOut.NLAG,dataOut.NRANGE,dataOut.NAVG),'complex128')

            #self.lagp4=numpy.zeros((dataOut.NLAG,dataOut.NRANGE,dataOut.NAVG),'complex64')
            self.aux_cross_lp=0

        for i in range(dataOut.NR):
            buffer_dc=dataOut.dc[i]
            for j in range(dataOut.NRANGE):

                range_for_n=numpy.min((dataOut.NRANGE-j,dataOut.NLAG))

                buffer_aux=numpy.conj(buffer[i,:dataOut.nProfiles,j])#-buffer_dc)
                for n in range(range_for_n):

                    c=(buffer_aux)*(buffer[i,:dataOut.nProfiles,j+n])#-buffer_dc)
                    #print(c.dtype)
                    #print(self.lagp0.dtype)
                    #exit(1)
                    if i==0:
                        self.lagp0[n][j][self.bcounter-1]=numpy.sum(c[:dataOut.NSCAN])
                        #ww[n,j,:,self.bcounter-1]=c[:dataOut.NSCAN]
                        self.lagp3[n][j][self.bcounter-1]=numpy.sum(c[dataOut.NSCAN:]/self.cnorm)
                    elif i==1:
                        self.lagp1[n][j][self.bcounter-1]=numpy.sum(c[:dataOut.NSCAN])
                    elif i==2:
                        self.lagp2[n][j][self.bcounter-1]=numpy.sum(c[:dataOut.NSCAN])
                    #if i==2 and j==30 and n==0:
                        #print(c[:dataOut.NSCAN])
                        #print(numpy.sum(c[:dataOut.NSCAN]))

        self.lagp0[:,:,self.bcounter-1]=numpy.conj(self.lagp0[:,:,self.bcounter-1])
        self.lagp1[:,:,self.bcounter-1]=numpy.conj(self.lagp1[:,:,self.bcounter-1])
        self.lagp2[:,:,self.bcounter-1]=numpy.conj(self.lagp2[:,:,self.bcounter-1])
        self.lagp3[:,:,self.bcounter-1]=numpy.conj(self.lagp3[:,:,self.bcounter-1])
        #print(self.lagp3[2,197,0])
        #print(self.lagp0[0,0,self.bcounter-1])
        #print(sum(self.buffer[3,:,199,2]))
        #print(self.cnorm)
        #exit(1)
        #print("self,lagp0: ", self.lagp0[0,0,self.bcounter-1])
        #print(ww[:,0,0,self.bcounter-1])
        #exit(1)


    def LP_median_estimates_original(self,dataOut):

        if self.bcounter==dataOut.NAVG:

            if self.lag_products_LP_median_estimates_aux==1:
                self.output=numpy.zeros((dataOut.NLAG,dataOut.NRANGE,dataOut.NR),'complex128')
                self.lag_products_LP_median_estimates_aux=0

                #print("self,lagp0: ", numpy.sum(self.lagp0[0,0,:]))
            for i in range(dataOut.NLAG):
                for j in range(dataOut.NRANGE):
                    for l in range(4): #four outputs
                        for k in range(dataOut.NAVG):
                            if k==0:
                                self.output[i,j,l]=0.0+0.j
                                if l==0:
                                    self.lagp0[i,j,:]=sorted(self.lagp0[i,j,:], key=lambda x: x.real)  #sorted(self.lagp0[i,j,:].real)
                                if l==1:
                                    self.lagp1[i,j,:]=sorted(self.lagp1[i,j,:], key=lambda x: x.real)  #sorted(self.lagp1[i,j,:].real)
                                if l==2:
                                    self.lagp2[i,j,:]=sorted(self.lagp2[i,j,:], key=lambda x: x.real)  #sorted(self.lagp2[i,j,:].real)
                                if l==3:
                                    self.lagp3[i,j,:]=sorted(self.lagp3[i,j,:], key=lambda x: x.real)  #sorted(self.lagp3[i,j,:].real)
                                #print(self.lagp0[2,100,1]);exit(1)
                            if k>=dataOut.nkill/2 and k<dataOut.NAVG-dataOut.nkill/2:
                                if l==0:
                                    self.output[i,j,l]=self.output[i,j,l]+((float(dataOut.NAVG)/(float)(dataOut.NAVG-dataOut.nkill))*self.lagp0[i,j,k])
                                if l==1:
                                    self.output[i,j,l]=self.output[i,j,l]+((float(dataOut.NAVG)/(float)(dataOut.NAVG-dataOut.nkill))*self.lagp1[i,j,k])
                                if l==2:
                                    self.output[i,j,l]=self.output[i,j,l]+((float(dataOut.NAVG)/(float)(dataOut.NAVG-dataOut.nkill))*self.lagp2[i,j,k])
                                if l==3:
                                    self.output[i,j,l]=self.output[i,j,l]+((float(dataOut.NAVG)/(float)(dataOut.NAVG-dataOut.nkill))*self.lagp3[i,j,k])
                                #if j== 30 and i==2 and l==0:
                                    #print(self.output[2,30,0])
            dataOut.output_LP=self.output
            dataOut.data_for_RTI_LP=numpy.zeros((4,dataOut.NRANGE))
            dataOut.data_for_RTI_LP[0],dataOut.data_for_RTI_LP[1],dataOut.data_for_RTI_LP[2],dataOut.data_for_RTI_LP[3]=self.RTI_LP(dataOut.output_LP,dataOut.NRANGE)
            #print("output:",self.output[2,30,0])
            #exit(1)

    def LP_median_estimates(self,dataOut):

        if self.bcounter==dataOut.NAVG:

            if self.lag_products_LP_median_estimates_aux==1:
                self.output=numpy.zeros((dataOut.NLAG,dataOut.NRANGE,dataOut.NR),'complex128')
                self.lag_products_LP_median_estimates_aux=0
                #print("self,lagp0: ", numpy.sum(self.lagp0[0,0,:]))

            for i in range(dataOut.NLAG):
            #my_list = ([0,1,2,3,4,5,6,7]) #hasta 7 funciona, en 6 ya no
            #for i in my_list:
                for j in range(dataOut.NRANGE):
                    for l in range(4): #four outputs
                        for k in range(dataOut.NAVG):
                            if k==0:
                                self.output[i,j,l]=0.0+0.j
                                if l==0:
                                    self.lagp0[i,j,:]=sorted(self.lagp0[i,j,:], key=lambda x: x.real)  #sorted(self.lagp0[i,j,:].real)
                                if l==1:
                                    self.lagp1[i,j,:]=sorted(self.lagp1[i,j,:], key=lambda x: x.real)  #sorted(self.lagp1[i,j,:].real)
                                if l==2:
                                    self.lagp2[i,j,:]=sorted(self.lagp2[i,j,:], key=lambda x: x.real)  #sorted(self.lagp2[i,j,:].real)
                                if l==3:
                                    self.lagp3[i,j,:]=sorted(self.lagp3[i,j,:], key=lambda x: x.real)  #sorted(self.lagp3[i,j,:].real)
                            '''
                            x = 2
                            if k>=x and k<dataOut.NAVG-dataOut.nkill/2:
                                if l==0:
                                    self.output[i,j,l]=self.output[i,j,l]+((float(dataOut.NAVG)/(float)(dataOut.NAVG-x-dataOut.nkill/2))*self.lagp0[i,j,k])
                                if l==1:
                                    self.output[i,j,l]=self.output[i,j,l]+((float(dataOut.NAVG)/(float)(dataOut.NAVG-x-dataOut.nkill/2))*self.lagp1[i,j,k])
                                if l==2:
                                    self.output[i,j,l]=self.output[i,j,l]+((float(dataOut.NAVG)/(float)(dataOut.NAVG-x-dataOut.nkill/2))*self.lagp2[i,j,k])
                                if l==3:
                                    self.output[i,j,l]=self.output[i,j,l]+((float(dataOut.NAVG)/(float)(dataOut.NAVG-x-dataOut.nkill/2))*self.lagp3[i,j,k])
                                    '''
                            #'''
                            if k>=dataOut.nkill/2 and k<dataOut.NAVG-dataOut.nkill/2:
                                if l==0:
                                    self.output[i,j,l]=self.output[i,j,l]+((float(dataOut.NAVG)/(float)(dataOut.NAVG-dataOut.nkill))*self.lagp0[i,j,k])
                                if l==1:
                                    self.output[i,j,l]=self.output[i,j,l]+((float(dataOut.NAVG)/(float)(dataOut.NAVG-dataOut.nkill))*self.lagp1[i,j,k])
                                if l==2:
                                    self.output[i,j,l]=self.output[i,j,l]+((float(dataOut.NAVG)/(float)(dataOut.NAVG-dataOut.nkill))*self.lagp2[i,j,k])
                                if l==3:
                                    self.output[i,j,l]=self.output[i,j,l]+((float(dataOut.NAVG)/(float)(dataOut.NAVG-dataOut.nkill))*self.lagp3[i,j,k])
                                    #'''


            dataOut.output_LP=self.output
            dataOut.data_for_RTI_LP=numpy.zeros((4,dataOut.NRANGE))
            dataOut.data_for_RTI_LP[0],dataOut.data_for_RTI_LP[1],dataOut.data_for_RTI_LP[2],dataOut.data_for_RTI_LP[3]=self.RTI_LP(dataOut.output_LP,dataOut.NRANGE)

    def get_dc(self,dataOut):

        if self.bcounter==0:
            dataOut.dc=numpy.zeros(dataOut.NR,dtype='complex64')

        dataOut.dc+=numpy.sum(dataOut.data[:,:,2*dataOut.NLAG:dataOut.NRANGE],axis=(1,2))
        dataOut.dc=dataOut.dc/float(dataOut.nProfiles*(dataOut.NRANGE-2*dataOut.NLAG))

    def noise_estimation4x_HP(self,dataOut):
        if self.bcounter==dataOut.NAVG:
            dataOut.noise_final=numpy.zeros(dataOut.NR,'float32')
            #snoise=numpy.zeros((NR,NAVG),'float32')
            #nvector1=numpy.zeros((NR,NAVG,MAXNRANGENDT),'float32')
            sorted_data=numpy.zeros((dataOut.MAXNRANGENDT,dataOut.NR,dataOut.NAVG),'float32')
            for i in range(dataOut.NR):
                dataOut.noise_final[i]=0.0
                for j in range(dataOut.MAXNRANGENDT):
                    sorted_data[j,i,:]=numpy.copy(sorted(dataOut.noisevector[j,i,:]))
                    l=dataOut.MAXNRANGENDT-2
                    for k in range(dataOut.NAVG):
                        if k>=dataOut.nkill/2 and k<dataOut.NAVG-dataOut.nkill/2:
                            dataOut.noise_final[i]+=sorted_data[min(j,l),i,k]*float(dataOut.NAVG)/float(dataOut.NAVG-dataOut.nkill)

    def noisevectorizer(self,NSCAN,nProfiles,NR,MAXNRANGENDT,noisevector,data,dc):

        #rnormalizer= 1./(float(nProfiles - NSCAN))
        rnormalizer= float(NSCAN)/((float(nProfiles - NSCAN))*float(MAXNRANGENDT))
        for i in range(NR):
            for j in range(MAXNRANGENDT):
                for k in range(NSCAN,nProfiles):
                    #TODO:integrate just 2nd quartile gates
                    if k==NSCAN:
                        noisevector[j][i][self.bcounter]=(abs(data[i][k][j]-dc[i])**2)*rnormalizer
                        ##noisevector[j][i][iavg]=(abs(cdata[k][j][i])**2)*rnormalizer
                    else:
                        noisevector[j][i][self.bcounter]+=(abs(data[i][k][j]-dc[i])**2)*rnormalizer


    def RTI_LP(self,output,NRANGE):
        x00=numpy.zeros(NRANGE,dtype='float32')
        x01=numpy.zeros(NRANGE,dtype='float32')
        x02=numpy.zeros(NRANGE,dtype='float32')
        x03=numpy.zeros(NRANGE,dtype='float32')

        for i in range(2): #first couple of lags
            for j in range(NRANGE): #
                #fx=numpy.sqrt((kaxbx[i,j,k]+kayby[i,j,k])**2+(kaybx[i,j,k]-kaxby[i,j,k])**2)
                x00[j]+=numpy.abs(output[i,j,0]) #Ch0
                x01[j]+=numpy.abs(output[i,j,1]) #Ch1
                x02[j]+=numpy.abs(output[i,j,2]) #Ch2
                x03[j]+=numpy.abs(output[i,j,3]) #Ch3
                #x02[i]=x02[i]+fx

                x00[j]=10.0*numpy.log10(x00[j]/2.)
                x01[j]=10.0*numpy.log10(x01[j]/2.)
                x02[j]=10.0*numpy.log10(x02[j]/2.)
                x03[j]=10.0*numpy.log10(x03[j]/2.)
                #x02[i]=10.0*numpy.log10(x02[i])
        return x00,x01,x02,x03

    def run(self, dataOut, NLAG=16, NRANGE=200, NCAL=0, DPL=11,
        NDN=0, NDT=67, NDP=67, NSCAN=128,
        lagind=(0,1,2,3,4,5,6,7,0,3,4,5,6,8,9,10), lagfirst=(1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1),
        NAVG=16, nkill=6):
        #print(dataOut.data[1,:12,200:200+15])
        #exit(1)
        dataOut.NLAG=NLAG
        dataOut.NR=len(dataOut.channelList)
        dataOut.NRANGE=NRANGE
        dataOut.NCAL=NCAL
        dataOut.DPL=DPL
        dataOut.NDN=NDN
        dataOut.NDT=NDT
        dataOut.NDP=NDP
        dataOut.NSCAN=NSCAN
        dataOut.DH=dataOut.heightList[1]-dataOut.heightList[0]
        dataOut.H0=int(dataOut.heightList[0])
        dataOut.lagind=lagind
        dataOut.lagfirst=lagfirst
        dataOut.NAVG=NAVG
        dataOut.nkill=nkill

        dataOut.flagNoData =  True

        self.get_dc(dataOut)
        self.get_products_cabxys_HP(dataOut)
        self.cabxys_navg(dataOut)
        self.lag_products_LP(dataOut)
        self.LP_median_estimates(dataOut)
        self.noise_estimation4x_HP(dataOut)
        self.kabxys(dataOut)

        return dataOut


class RemoveDebris(Operation):
    """Operation to remove blocks where an outlier is found for Double (Long) Pulse.

    Parameters:
    -----------
    None

    Example
    --------

    op = proc_unit.addOperation(name='RemoveDebris', optype='other')

    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

    def run(self,dataOut):
        debris=numpy.zeros(dataOut.NRANGE,'float32')

        for j in range(0,3):
            for i in range(dataOut.NRANGE):
                if j==0:
                    debris[i]=10*numpy.log10(numpy.abs(dataOut.output_LP[j,i,0]))
                else:
                    debris[i]+=10*numpy.log10(numpy.abs(dataOut.output_LP[j,i,0]))

        thresh=8.0+4+4+4
        for i in range(47,100):
            if ((debris[i-2]+debris[i-1]+debris[i]+debris[i+1])>
                ((debris[i-12]+debris[i-11]+debris[i-10]+debris[i-9]+
                debris[i+12]+debris[i+11]+debris[i+10]+debris[i+9])/2.0+
                thresh)):

                dataOut.flagNoData=True
                print("LP Debris detected at",i*15,"km")

        debris=numpy.zeros(dataOut.NDP,dtype='float32')
        Range=numpy.arange(0,3000,15)
        for k in range(2): #flip
            for i in range(dataOut.NDP): #
                debris[i]+=numpy.sqrt((dataOut.kaxbx[i,0,k]+dataOut.kayby[i,0,k])**2+(dataOut.kaybx[i,0,k]-dataOut.kaxby[i,0,k])**2)

        if gmtime(dataOut.utctime).tm_hour > 11:
            for i in range(2,dataOut.NDP-2):
                if (debris[i]>3.0*debris[i-2] and
                    debris[i]>3.0*debris[i+2] and
                    Range[i]>200.0 and Range[i]<=540.0):
                    dataOut.flagNoData=True
                    print("DP Debris detected at",i*15,"km")

        return dataOut


class IntegrationHP(IntegrationDP):
    """Operation to integrate Double Pulse and Long Pulse data.

    Parameters:
    -----------
    nint : int
        Number of integrations.

    Example
    --------

    op = proc_unit.addOperation(name='IntegrationHP', optype='other')
    op.addParameter(name='nint', value='30', format='int')

    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

        self.counter = 0
        self.aux = 0

    def integration_noise(self,dataOut):

        if self.counter == 0:
            dataOut.tnoise=numpy.zeros((dataOut.NR),dtype='float32')

        dataOut.tnoise+=dataOut.noise_final

    def integration_for_long_pulse(self,dataOut):

        if self.counter == 0:
            dataOut.output_LP_integrated=numpy.zeros((dataOut.NLAG,dataOut.NRANGE,dataOut.NR),order='F',dtype='complex128')

        dataOut.output_LP_integrated+=dataOut.output_LP

    def run(self,dataOut,nint=None):

        dataOut.flagNoData=True

        dataOut.nint=nint
        dataOut.paramInterval=0#int(dataOut.nint*dataOut.header[7][0]*2 )
        dataOut.lat=-11.95
        dataOut.lon=-76.87

        self.integration_for_long_pulse(dataOut)

        self.integration_noise(dataOut)

        if self.counter==dataOut.nint-1:
            dataOut.nis=dataOut.NSCAN*dataOut.NAVG*dataOut.nint*10
            #print(dataOut.tnoise)
            #exit(1)
            dataOut.tnoise[0]*=0.995
            dataOut.tnoise[1]*=0.995
            dataOut.pan=dataOut.tnoise[0]/float(dataOut.NSCAN*dataOut.nint*dataOut.NAVG)
            dataOut.pbn=dataOut.tnoise[1]/float(dataOut.NSCAN*dataOut.nint*dataOut.NAVG)
        #print(self.counter) ToDo: Fix when nint = 1
            '''
            print("pan: ",dataOut.pan)
            print("pbn: ",dataOut.pbn)
            #print("tnoise: ",dataOut.tnoise)
            exit(1)
            '''
            #print(dataOut.output_LP_integrated[0,30,2])
            #exit(1)
        self.integration_for_double_pulse(dataOut)
        #if self.counter==dataOut.nint:
        #print(dataOut.kabxys_integrated[8][53,6,0]+dataOut.kabxys_integrated[11][53,6,0])
        #print(dataOut.kabxys_integrated[8][53,9,0]+dataOut.kabxys_integrated[11][53,9,0])
        #exit(1)
        #print(dataOut.flagNoData)
        return dataOut

class SumFlipsHP(SumFlips):
    """Operation to sum the flip and unflip part of certain cross products of the Double Pulse.

    Parameters:
    -----------
    None

    Example
    --------

    op = proc_unit.addOperation(name='SumFlipsHP', optype='other')

    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)

    def rint2HP(self,dataOut):

        dataOut.rnint2=numpy.zeros(dataOut.DPL,'float32')
        #print(dataOut.nint,dataOut.NAVG)
        for l in range(dataOut.DPL):
            if(l==0 or (l>=3 and l <=6)):
                dataOut.rnint2[l]=0.5/float(dataOut.nint*dataOut.NAVG*16.0)
            else:
                dataOut.rnint2[l]=0.5/float(dataOut.nint*dataOut.NAVG*8.0)

    def run(self,dataOut):

        self.rint2HP(dataOut)
        self.SumLags(dataOut)

        hei = 2
        lag = 0
        '''
        for hei in range(67):
            print("hei",hei)
            print(dataOut.kabxys_integrated[8][hei,:,0]+dataOut.kabxys_integrated[11][hei,:,0])
            print(dataOut.kabxys_integrated[10][hei,:,0]-dataOut.kabxys_integrated[9][hei,:,0])
        exit(1)
        '''
        '''
        print("b",(dataOut.kabxys_integrated[4][hei,lag,0]+dataOut.kabxys_integrated[5][hei,lag,0]))
        print((dataOut.kabxys_integrated[6][hei,lag,0]+dataOut.kabxys_integrated[7][hei,lag,0]))
        print("c",(dataOut.kabxys_integrated[8][hei,lag,0]+dataOut.kabxys_integrated[11][hei,lag,0]))
        print((dataOut.kabxys_integrated[10][hei,lag,0]-dataOut.kabxys_integrated[9][hei,lag,0]))
        exit(1)
        '''
        #print(dataOut.rnint2)
        #print(numpy.sum(dataOut.kabxys_integrated[4][:,1,0]+dataOut.kabxys_integrated[5][:,1,0]))
        #print(dataOut.nis)
        #exit(1)
        return dataOut


class LongPulseAnalysis(Operation):
    """Operation to estimate ACFs, temperatures, total electron density and Hydrogen/Helium fractions from the Long Pulse data.

    Parameters:
    -----------
    NACF : int
        .*

    Example
    --------

    op = proc_unit.addOperation(name='LongPulseAnalysis', optype='other')
    op.addParameter(name='NACF', value='16', format='int')

    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)
        self.aux=1

    def run(self,dataOut,NACF):

        dataOut.NACF=NACF
        dataOut.heightList=dataOut.DH*(numpy.arange(dataOut.NACF))
        anoise0=dataOut.tnoise[0]
        anoise1=anoise0*0.0       #seems to be noise in 1st lag 0.015 before '14
        #print(anoise0)
        #exit(1)
        if self.aux:
            #dataOut.cut=31#26#height=31*15=465
            self.cal=numpy.zeros((dataOut.NLAG),'float32')
            self.drift=numpy.zeros((200),'float32')
            self.rdrift=numpy.zeros((200),'float32')
            self.ddrift=numpy.zeros((200),'float32')
            self.sigma=numpy.zeros((dataOut.NRANGE),order='F',dtype='float32')
            self.powera=numpy.zeros((dataOut.NRANGE),order='F',dtype='float32')
            self.powerb=numpy.zeros((dataOut.NRANGE),order='F',dtype='float32')
            self.perror=numpy.zeros((dataOut.NRANGE),order='F',dtype='float32')
            dataOut.ene=numpy.zeros((dataOut.NRANGE),'float32')
            self.dpulse=numpy.zeros((dataOut.NACF),'float32')
            self.lpulse=numpy.zeros((dataOut.NACF),'float32')
            dataOut.lags_LP=numpy.zeros((dataOut.IBITS),order='F',dtype='float32')
            self.lagp=numpy.zeros((dataOut.NACF),'float32')
            self.u=numpy.zeros((2*dataOut.NACF,2*dataOut.NACF),'float32')
            dataOut.ne=numpy.zeros((dataOut.NRANGE),order='F',dtype='float32')
            dataOut.te=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.ete=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.ti=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.eti=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.ph=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.eph=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.phe=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.ephe=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.errors=numpy.zeros((dataOut.IBITS,max(dataOut.NRANGE,dataOut.NSHTS)),order='F',dtype='float32')
            dataOut.fit_array_real=numpy.zeros((max(dataOut.NRANGE,dataOut.NSHTS),dataOut.NLAG),order='F',dtype='float32')
            dataOut.status=numpy.zeros(1,'float32')
            dataOut.tx=240.0 #debería provenir del header #hybrid
            
            for i in range(dataOut.IBITS):
                dataOut.lags_LP[i]=float(i)*(dataOut.tx/150.0)/float(dataOut.IBITS) # (float)i*(header.tx/150.0)/(float)IBITS;

            self.aux=0

        ## Juntion point between DP and LP
        dataOut.cut=30
        for i in range(30,15,-1):
            if numpy.nanmax(dataOut.acfs_error_to_plot[i,:])>=10 or dataOut.info2[i]==0:
                dataOut.cut=i-1

        ##
        for i in range(dataOut.NLAG):
            self.cal[i]=sum(dataOut.output_LP_integrated[i,:,3].real) #Lag x Height x Channel
        #print('*', dataOut.output_LP_integrated, numpy.shape(dataOut.output_LP_integrated), dataOut.output_LP_integrated[0,:,0]); exit(1)
        #print(numpy.sum(self.cal)) #Coinciden
        #exit(1)
        self.cal/=float(dataOut.NRANGE)
        #print("anoise0", anoise0)
        #print("anoise1", anoise1)


        #exit(1)
        #print("nis: ", dataOut.nis)
        #print("pan: ", dataOut.pan)
        #print("pbn: ", dataOut.pbn)
        #print(numpy.sum(dataOut.output_LP_integrated[0,:,0]))
        
        '''import matplotlib.pyplot as plt
        plt.plot(dataOut.output_LP_integrated[:,40,0])
        plt.show()'''
        
        #print(dataOut.output_LP_integrated[0,40,0])
        #print(numpy.sum(dataOut.output_LP_integrated[:,0,0]))
        #exit(1)

        #################### PROBAR MÁS INTEGRACIÓN, SINO MODIFICAR VALOR DE "NIS" ####################
                                    # VER dataOut.nProfiles_LP #

        
        #PLOTEAR POTENCIA VS RUIDO, QUIZA SE ESTA REMOVIENDO MUCHA SEÑAL

        '''import matplotlib.pyplot as plt
        plt.plot(10*numpy.log10(dataOut.output_LP_integrated.real[0,:,0]),dataOut.range1)
        #plt.plot(10*numpy.log10(dataOut.output_LP_integrated.real[0,:,0]/dataOut.nProfiles_LP),dataOut.range1)
        plt.axvline(10*numpy.log10(anoise0),color='k',linestyle='dashed')
        plt.grid()
        plt.xlim(20,100)
        plt.show()'''
        

        ## Noise subtract 0th and 1st lag, power profile ==> powera

        for j in range(dataOut.NACF+2*dataOut.IBITS+2):

            dataOut.output_LP_integrated.real[0,j,0]-=anoise0   #lag0 ch0
            dataOut.output_LP_integrated.real[1,j,0]-=anoise1   #lag1 ch0

            for i in range(1,dataOut.NLAG):  #remove cal data from certain lags
                 dataOut.output_LP_integrated.real[i,j,0]-=self.cal[i]
            k=max(j,26)   #constant power below range 26
            self.powera[j]=dataOut.output_LP_integrated.real[0,k,0] #Lag0 and Channel 0
        #print(numpy.sum(self.powera))
        print("***", dataOut.NACF, dataOut.IBITS, dataOut.NACF+2*dataOut.IBITS+2)
        print(numpy.shape(dataOut.output_LP_integrated))


        ## Examine drifts here - based on 60 'indep.' estimates

        #nis=dataOut.NSCAN*dataOut.NAVG*dataOut.nint*10
        nis = dataOut.nis
        #print("nis",nis)
        alpha=beta=delta=0.0
        nest=0
        gamma=3.0/(2.0*numpy.pi*dataOut.lags_LP[1]*1.0e-3)
        beta=gamma*(math.atan2(dataOut.output_LP_integrated.imag[14,0,2],dataOut.output_LP_integrated.real[14,0,2])-math.atan2(dataOut.output_LP_integrated.imag[1,0,2],dataOut.output_LP_integrated.real[1,0,2]))/13.0
        #print(gamma,beta)

        for i in range(1,3):
            gamma=3.0/(2.0*numpy.pi*dataOut.lags_LP[i]*1.0e-3)
            #print("gamma",gamma)
            for j in range(34,44):
                rho2=numpy.abs(dataOut.output_LP_integrated[i,j,0])/numpy.abs(dataOut.output_LP_integrated[0,j,0])
                dataOut.dphi2=(1.0/rho2-1.0)/(float(2*nis))
                dataOut.dphi2*=gamma**2
                pest=gamma*math.atan(dataOut.output_LP_integrated.imag[i,j,0]/dataOut.output_LP_integrated.real[i,j,0])
                #print("1",dataOut.output_LP_integrated.imag[i,j,0])
                #print("2",dataOut.output_LP_integrated.real[i,j,0])
                self.drift[nest]=pest
                self.ddrift[nest]=dataOut.dphi2
                self.rdrift[nest]=float(nest)
                nest+=1

        ## Rank the drift estimates
        sorted(self.drift[:nest])
        for j in range(int(nest/4),int(3*nest/4)):
            #i=int(self.rdrift[j])
            alpha+=self.drift[j]/self.ddrift[j]
            delta+=1.0/self.ddrift[j]

        alpha/=delta
        delta=1./numpy.sqrt(delta)
        vdrift=alpha-beta
        dvdrift=delta

        '''
        Hysell note: need to develop estimate of complete density profile using all
        available data
        '''

        ## Estimatate sample variances for long-pulse power profile
        #nis=dataOut.NSCAN*dataOut.NAVG*dataOut.nint
        nis = dataOut.nis/10
        self.sigma[:dataOut.NACF+2*dataOut.IBITS+2]=((anoise0+self.powera[:dataOut.NACF+2*dataOut.IBITS+2])**2)/float(nis)
        ioff=1
        #print("nis",nis)
        #print(self.sigma)
        #exit(1)
        
        '''if(ut>16.0) // needs to be re-evaluated  # ask Roberto why not included
        ioff=1;
        else if(ut>14.0)
        ioff=2;  
        else
        ioff=1;'''


        #deconvolve rectangular pulse shape from profile ==> powerb, perror
        ############# START nnlswrap#############
        
        if dataOut.ut_Faraday>14.0:
            alpha_nnlswrap=20.0
        else:
            alpha_nnlswrap=30.0

        range1_nnls=dataOut.NACF
        range2_nnls=dataOut.NACF+dataOut.IBITS-1

        g_nnlswrap=numpy.zeros((range1_nnls,range2_nnls),'float32')
        a_nnlswrap=numpy.zeros((range2_nnls,range2_nnls),'float64')

        for i in range(range1_nnls):
            for j in range(range2_nnls):
                if j>=i and j<i+dataOut.IBITS:
                    g_nnlswrap[i,j]=1.0
                else:
                    g_nnlswrap[i,j]=0.0

        a_nnlswrap[:]=numpy.matmul(numpy.transpose(g_nnlswrap),g_nnlswrap)

        numpy.fill_diagonal(a_nnlswrap,a_nnlswrap.diagonal()+alpha_nnlswrap**2)

                    #ERROR ANALYSIS#

        self.perror[:range2_nnls]=0.0
        self.perror[:range2_nnls]=numpy.matmul(1./(self.sigma[dataOut.IBITS+ioff:range1_nnls+dataOut.IBITS+ioff]),g_nnlswrap**2)
        self.perror[:range1_nnls]+=(alpha_nnlswrap**2)/(self.sigma[dataOut.IBITS+ioff:range1_nnls+dataOut.IBITS+ioff])
        self.perror[:range2_nnls]=1.00/self.perror[:range2_nnls]

        b_nnlswrap=numpy.zeros(range2_nnls,'float64')
        b_nnlswrap[:]=numpy.matmul(self.powera[dataOut.IBITS+ioff:range1_nnls+dataOut.IBITS+ioff],g_nnlswrap) #match filter alturas

        x_nnlswrap=numpy.zeros(range2_nnls,'float64')
        x_nnlswrap[:]=nnls(a_nnlswrap,b_nnlswrap)[0]

        self.powerb[:range2_nnls]=x_nnlswrap
        #print(self.powerb[40])
        #print(self.powerb[66])
        #exit(1)
        #############END nnlswrap#############

        #print(numpy.sum(numpy.sqrt(self.perror[0:dataOut.NACF])))
        #print(self.powerb[0:dataOut.NACF])
        #print(dataOut.NACF)

        ## Estimate relative error for deconvolved profile (scaling irrelevant)

        dataOut.ene[0:dataOut.NACF]=numpy.sqrt(self.perror[0:dataOut.NACF])/self.powerb[0:dataOut.NACF]
        #print(numpy.sum(dataOut.ene))

        ## Reconvolve long- and double-pulse power profiles with rectangular pulse
        aux=0
        for i in range(dataOut.IBITS,dataOut.NACF):
            self.dpulse[i]=self.lpulse[i]=0.0
            for j in range(dataOut.IBITS):
                k=int(i-j)
                if k<36-aux and k>16:
                    self.dpulse[i]+=dataOut.ph2[k]/dataOut.h2[k]
                elif k>=36-aux:
                    self.lpulse[i]+=self.powerb[k]
            self.lagp[i]=self.powera[i]

        ## Find scale factor that best merges profiles
        
        qi=sum(self.dpulse[32:dataOut.NACF]**2/(self.lagp[32:dataOut.NACF]+anoise0)**2)
        ri=sum((self.dpulse[32:dataOut.NACF]*self.lpulse[32:dataOut.NACF])/(self.lagp[32:dataOut.NACF]+anoise0)**2)
        si=sum((self.dpulse[32:dataOut.NACF]*self.lagp[32:dataOut.NACF])/(self.lagp[32:dataOut.NACF]+anoise0)**2)
        ui=sum(self.lpulse[32:dataOut.NACF]**2/(self.lagp[32:dataOut.NACF]+anoise0)**2)
        vi=sum((self.lpulse[32:dataOut.NACF]*self.lagp[32:dataOut.NACF])/(self.lagp[32:dataOut.NACF]+anoise0)**2)

        alpha=(si*ui-vi*ri)/(qi*ui-ri*ri)
        beta=(qi*vi-ri*si)/(qi*ui-ri*ri)

        ## Form density profile estimate, merging rescaled power profiles

        self.powerb[16:36-aux]=alpha*dataOut.ph2[16:36-aux]/dataOut.h2[16:36-aux]
        self.powerb[36-aux:dataOut.NACF]*=beta
        #print(dataOut.h2)
        #print(numpy.sum(alpha))
        #print(numpy.sum(dataOut.ph2))

        # Form Ne estimate, fill in error estimate at low altitudes

        dataOut.ene[0:36-aux]=dataOut.sdp2[0:36-aux]/dataOut.ph2[0:36-aux]
        dataOut.ne[:dataOut.NACF]=self.powerb[:dataOut.NACF]*dataOut.h2[:dataOut.NACF]/alpha
        #print(numpy.sum(self.powerb))
        #print(numpy.sum(dataOut.ene))
        #print(numpy.sum(dataOut.ne))

        ## Error propagation: store zero lag error covariance in u

        nis=dataOut.NSCAN*dataOut.NAVG*dataOut.nint/1   # DLH serious debris removal

        for i in range(dataOut.NACF):
            for j in range(i,dataOut.NACF):
                if j-i>=dataOut.IBITS:
                    self.u[i,j]=0.0
                else:
                    self.u[i,j]=dataOut.output_LP_integrated.real[j-i,i,0]**2/float(nis)
                    self.u[i,j]*=(anoise0+dataOut.output_LP_integrated.real[0,i,0])/dataOut.output_LP_integrated.real[0,i,0]
                    self.u[i,j]*=(anoise0+dataOut.output_LP_integrated.real[0,j,0])/dataOut.output_LP_integrated.real[0,j,0]

                self.u[j,i]=self.u[i,j]

        ## Error analyis for lag product matrix (diag), place in acf_err

        for i in range(dataOut.NACF):
            for j in range(dataOut.IBITS):
                if j==0:
                    dataOut.errors[0,i]=numpy.sqrt(self.u[i,i])
                else:
                    dataOut.errors[j,i]=numpy.sqrt(((dataOut.output_LP_integrated.real[0,i,0]+anoise0)*(dataOut.output_LP_integrated.real[0,i+j,0]+anoise0)+dataOut.output_LP_integrated.real[j,i,0]**2)/float(2*nis))
        '''
        print(numpy.sum(dataOut.output_LP_integrated))
        print(numpy.sum(dataOut.errors))
        print(numpy.sum(self.powerb))
        print(numpy.sum(dataOut.ne))
        print(numpy.sum(dataOut.lags_LP))
        print(numpy.sum(dataOut.thb))
        print(numpy.sum(dataOut.bfm))
        print(numpy.sum(dataOut.te))
        print(numpy.sum(dataOut.ete))
        print(numpy.sum(dataOut.ti))
        print(numpy.sum(dataOut.eti))
        print(numpy.sum(dataOut.ph))
        print(numpy.sum(dataOut.eph))
        print(numpy.sum(dataOut.phe))
        print(numpy.sum(dataOut.ephe))
        print(numpy.sum(dataOut.range1))
        print(numpy.sum(dataOut.ut))
        print(numpy.sum(dataOut.NACF))
        print(numpy.sum(dataOut.fit_array_real))
        print(numpy.sum(dataOut.status))
        print(numpy.sum(dataOut.NRANGE))
        print(numpy.sum(dataOut.IBITS))
        exit(1)
        '''
        '''
        print(dataOut.te2[13:16])
        print(numpy.sum(dataOut.te2))
        exit(1)
        '''

        ###################Correlation pulse and itself

        print("LP Estimation")
        with suppress_stdout_stderr():
            #pass
            full_profile_profile.profile(numpy.transpose(dataOut.output_LP_integrated,(2,1,0)),
                                        numpy.transpose(dataOut.errors),self.powerb,dataOut.ne,dataOut.lags_LP,dataOut.thb,
                                        dataOut.bfm,dataOut.te,dataOut.ete,dataOut.ti,dataOut.eti,dataOut.ph,dataOut.eph,
                                        dataOut.phe,dataOut.ephe,dataOut.range1,dataOut.ut,dataOut.NACF,dataOut.fit_array_real,
                                        dataOut.status,dataOut.NRANGE,dataOut.IBITS)

        print("status: ",dataOut.status)

        if dataOut.status>=3.5:
            dataOut.te[:]=numpy.nan
            dataOut.ete[:]=numpy.nan
            dataOut.ti[:]=numpy.nan
            dataOut.eti[:]=numpy.nan
            dataOut.ph[:]=numpy.nan
            dataOut.eph[:]=numpy.nan
            dataOut.phe[:]=numpy.nan
            dataOut.ephe[:]=numpy.nan

        return dataOut

class LongPulseAnalysisSpectra(Operation):
    """Operation to estimate ACFs, temperatures, total electron density and Hydrogen/Helium fractions from the Long Pulse data.

    Parameters:
    -----------
    NACF : int
        .*

    Example
    --------

    op = proc_unit.addOperation(name='LongPulseAnalysis', optype='other')
    op.addParameter(name='NACF', value='16', format='int')

    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)
        self.aux=1

    def run(self,dataOut,NACF):

        dataOut.NACF=NACF
        dataOut.heightList=dataOut.DH*(numpy.arange(dataOut.NACF))
        anoise0=dataOut.tnoise[0]
        anoise1=anoise0*0.0       #seems to be noise in 1st lag 0.015 before '14
        #print(anoise0)
        #exit(1)
        if self.aux:
            #dataOut.cut=31#26#height=31*15=465
            self.cal=numpy.zeros((dataOut.NLAG),'float32')
            self.drift=numpy.zeros((200),'float32')
            self.rdrift=numpy.zeros((200),'float32')
            self.ddrift=numpy.zeros((200),'float32')
            self.sigma=numpy.zeros((dataOut.NRANGE),order='F',dtype='float32')
            self.powera=numpy.zeros((dataOut.NRANGE),order='F',dtype='float32')
            self.powerb=numpy.zeros((dataOut.NRANGE),order='F',dtype='float32')
            self.perror=numpy.zeros((dataOut.NRANGE),order='F',dtype='float32')
            dataOut.ene=numpy.zeros((dataOut.NRANGE),'float32')
            self.dpulse=numpy.zeros((dataOut.NACF),'float32')
            self.lpulse=numpy.zeros((dataOut.NACF),'float32')
            dataOut.lags_LP=numpy.zeros((dataOut.IBITS),order='F',dtype='float32')
            self.lagp=numpy.zeros((dataOut.NACF),'float32')
            self.u=numpy.zeros((2*dataOut.NACF,2*dataOut.NACF),'float32')
            dataOut.ne=numpy.zeros((dataOut.NRANGE),order='F',dtype='float32')
            dataOut.te=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.ete=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.ti=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.eti=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.ph=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.eph=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.phe=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.ephe=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.errors=numpy.zeros((dataOut.IBITS,max(dataOut.NRANGE,dataOut.NSHTS)),order='F',dtype='float32')
            dataOut.fit_array_real=numpy.zeros((max(dataOut.NRANGE,dataOut.NSHTS),dataOut.NLAG),order='F',dtype='float32')
            dataOut.status=numpy.zeros(1,'float32')
            dataOut.tx=240.0 #debería provenir del header #hybrid

            for i in range(dataOut.IBITS):
                dataOut.lags_LP[i]=float(i)*(dataOut.tx/150.0)/float(dataOut.IBITS) # (float)i*(header.tx/150.0)/(float)IBITS;

            self.aux=0

        dataOut.cut=30
        for i in range(30,15,-1): #Aquí se calcula en donde se unirá DP y LP en la parte final
            if numpy.nanmax(dataOut.acfs_error_to_plot[i,:])>=10 or dataOut.info2[i]==0:
                dataOut.cut=i-1

        for i in range(dataOut.NLAG):
            self.cal[i]=sum(dataOut.output_LP_integrated[i,:,3].real) #Lag x Height x Channel

        #print(numpy.sum(self.cal)) #Coinciden
        #exit(1)
        self.cal/=float(dataOut.NRANGE)


        #################### PROBAR MÁS INTEGRACIÓN, SINO MODIFICAR VALOR DE "NIS" ####################
                                    # VER dataOut.nProfiles_LP #

        '''
        #PLOTEAR POTENCIA VS RUIDO, QUIZA SE ESTA REMOVIENDO MUCHA SEÑAL
        #print(dataOut.heightList)
        import matplotlib.pyplot as plt
        plt.plot(10*numpy.log10(dataOut.output_LP_integrated.real[0,:,0]),dataOut.range1)
        #plt.plot(10*numpy.log10(dataOut.output_LP_integrated.real[0,:,0]/dataOut.nProfiles_LP),dataOut.range1)
        plt.axvline(10*numpy.log10(anoise0),color='k',linestyle='dashed')
        plt.grid()
        plt.xlim(20,100)
        plt.show()
        '''


        for j in range(dataOut.NACF+2*dataOut.IBITS+2):

            dataOut.output_LP_integrated.real[0,j,0]-=anoise0   #lag0 ch0
            dataOut.output_LP_integrated.real[1,j,0]-=anoise1   #lag1 ch0

            for i in range(1,dataOut.NLAG):  #remove cal data from certain lags
                 dataOut.output_LP_integrated.real[i,j,0]-=self.cal[i]
            k=max(j,26)   #constant power below range 26
            self.powera[j]=dataOut.output_LP_integrated.real[0,k,0] #Lag0 and Channel 0

            ## examine drifts here - based on 60 'indep.' estimates
        #print(numpy.sum(self.powera))
        #exit(1)
        #nis=dataOut.NSCAN*dataOut.NAVG*dataOut.nint*10
        nis = dataOut.nis
        #print("nis",nis)
        alpha=beta=delta=0.0
        nest=0
        gamma=3.0/(2.0*numpy.pi*dataOut.lags_LP[1]*1.0e-3)
        beta=gamma*(math.atan2(dataOut.output_LP_integrated.imag[14,0,2],dataOut.output_LP_integrated.real[14,0,2])-math.atan2(dataOut.output_LP_integrated.imag[1,0,2],dataOut.output_LP_integrated.real[1,0,2]))/13.0
        #print(gamma,beta)
        #exit(1)
        for i in range(1,3):
            gamma=3.0/(2.0*numpy.pi*dataOut.lags_LP[i]*1.0e-3)
            #print("gamma",gamma)
            for j in range(34,44):
                rho2=numpy.abs(dataOut.output_LP_integrated[i,j,0])/numpy.abs(dataOut.output_LP_integrated[0,j,0])
                dataOut.dphi2=(1.0/rho2-1.0)/(float(2*nis))
                dataOut.dphi2*=gamma**2
                pest=gamma*math.atan(dataOut.output_LP_integrated.imag[i,j,0]/dataOut.output_LP_integrated.real[i,j,0])
                #print("1",dataOut.output_LP_integrated.imag[i,j,0])
                #print("2",dataOut.output_LP_integrated.real[i,j,0])
                self.drift[nest]=pest
                self.ddrift[nest]=dataOut.dphi2
                self.rdrift[nest]=float(nest)
                nest+=1

        sorted(self.drift[:nest])

        #print(dataOut.dphi2)
        #exit(1)

        for j in range(int(nest/4),int(3*nest/4)):
            #i=int(self.rdrift[j])
            alpha+=self.drift[j]/self.ddrift[j]
            delta+=1.0/self.ddrift[j]

        alpha/=delta
        delta=1./numpy.sqrt(delta)
        vdrift=alpha-beta
        dvdrift=delta

        #need to develop estimate of complete density profile using all
        #available data

        #estimate sample variances for long-pulse power profile

        #nis=dataOut.NSCAN*dataOut.NAVG*dataOut.nint
        nis = dataOut.nis/10
        #print("nis",nis)

        self.sigma[:dataOut.NACF+2*dataOut.IBITS+2]=((anoise0+self.powera[:dataOut.NACF+2*dataOut.IBITS+2])**2)/float(nis)
        #print(self.sigma)
        #exit(1)
        ioff=1

        #deconvolve rectangular pulse shape from profile ==> powerb, perror

        '''
        ############# START nnlswrap#############

        if dataOut.ut_Faraday>14.0:
            alpha_nnlswrap=20.0
        else:
            alpha_nnlswrap=30.0

        range1_nnls=dataOut.NACF
        range2_nnls=dataOut.NACF+dataOut.IBITS-1

        g_nnlswrap=numpy.zeros((range1_nnls,range2_nnls),'float32')
        a_nnlswrap=numpy.zeros((range2_nnls,range2_nnls),'float64')

        for i in range(range1_nnls):
            for j in range(range2_nnls):
                if j>=i and j<i+dataOut.IBITS:
                    g_nnlswrap[i,j]=1.0
                else:
                    g_nnlswrap[i,j]=0.0

        a_nnlswrap[:]=numpy.matmul(numpy.transpose(g_nnlswrap),g_nnlswrap)

        numpy.fill_diagonal(a_nnlswrap,a_nnlswrap.diagonal()+alpha_nnlswrap**2)

                    #ERROR ANALYSIS#

        self.perror[:range2_nnls]=0.0
        self.perror[:range2_nnls]=numpy.matmul(1./(self.sigma[dataOut.IBITS+ioff:range1_nnls+dataOut.IBITS+ioff]),g_nnlswrap**2)
        self.perror[:range1_nnls]+=(alpha_nnlswrap**2)/(self.sigma[dataOut.IBITS+ioff:range1_nnls+dataOut.IBITS+ioff])
        self.perror[:range2_nnls]=1.00/self.perror[:range2_nnls]

        b_nnlswrap=numpy.zeros(range2_nnls,'float64')
        b_nnlswrap[:]=numpy.matmul(self.powera[dataOut.IBITS+ioff:range1_nnls+dataOut.IBITS+ioff],g_nnlswrap)

        x_nnlswrap=numpy.zeros(range2_nnls,'float64')
        x_nnlswrap[:]=nnls(a_nnlswrap,b_nnlswrap)[0]

        self.powerb[:range2_nnls]=x_nnlswrap
        #print(self.powerb[40])
        #print(self.powerb[66])
        #exit(1)
        #############END nnlswrap#############
        '''
        self.powerb[:] = self.powera
        self.perror[:] = 0.
        #print(numpy.sum(numpy.sqrt(self.perror[0:dataOut.NACF])))
        #print(self.powerb[0:dataOut.NACF])
        #exit(1)
        #estimate relative error for deconvolved profile (scaling irrelevant)
        #print(dataOut.NACF)
        dataOut.ene[0:dataOut.NACF]=numpy.sqrt(self.perror[0:dataOut.NACF])/self.powerb[0:dataOut.NACF]
        #print(numpy.sum(dataOut.ene))
        #exit(1)
        aux=0

        for i in range(dataOut.IBITS,dataOut.NACF):
            self.dpulse[i]=self.lpulse[i]=0.0
            for j in range(dataOut.IBITS):
                k=int(i-j)
                if k<36-aux and k>16:
                    self.dpulse[i]+=dataOut.ph2[k]/dataOut.h2[k]
                elif k>=36-aux:
                    self.lpulse[i]+=self.powerb[k]
            self.lagp[i]=self.powera[i]

        #find scale factor that best merges profiles

        qi=sum(self.dpulse[32:dataOut.NACF]**2/(self.lagp[32:dataOut.NACF]+anoise0)**2)
        ri=sum((self.dpulse[32:dataOut.NACF]*self.lpulse[32:dataOut.NACF])/(self.lagp[32:dataOut.NACF]+anoise0)**2)
        si=sum((self.dpulse[32:dataOut.NACF]*self.lagp[32:dataOut.NACF])/(self.lagp[32:dataOut.NACF]+anoise0)**2)
        ui=sum(self.lpulse[32:dataOut.NACF]**2/(self.lagp[32:dataOut.NACF]+anoise0)**2)
        vi=sum((self.lpulse[32:dataOut.NACF]*self.lagp[32:dataOut.NACF])/(self.lagp[32:dataOut.NACF]+anoise0)**2)

        alpha=(si*ui-vi*ri)/(qi*ui-ri*ri)
        beta=(qi*vi-ri*si)/(qi*ui-ri*ri)

        #form density profile estimate, merging rescaled power profiles
        #print(dataOut.h2)
        #print(numpy.sum(alpha))
        #print(numpy.sum(dataOut.ph2))
        self.powerb[16:36-aux]=alpha*dataOut.ph2[16:36-aux]/dataOut.h2[16:36-aux]
        self.powerb[36-aux:dataOut.NACF]*=beta

        #form Ne estimate, fill in error estimate at low altitudes

        dataOut.ene[0:36-aux]=dataOut.sdp2[0:36-aux]/dataOut.ph2[0:36-aux]
        dataOut.ne[:dataOut.NACF]=self.powerb[:dataOut.NACF]*dataOut.h2[:dataOut.NACF]/alpha
        #print(numpy.sum(self.powerb))
        #print(numpy.sum(dataOut.ene))
        #print(numpy.sum(dataOut.ne))
        #exit(1)
        #now do error propagation: store zero lag error covariance in u

        nis=dataOut.NSCAN*dataOut.NAVG*dataOut.nint/1   # DLH serious debris removal

        for i in range(dataOut.NACF):
            for j in range(i,dataOut.NACF):
                if j-i>=dataOut.IBITS:
                    self.u[i,j]=0.0
                else:
                    self.u[i,j]=dataOut.output_LP_integrated.real[j-i,i,0]**2/float(nis)
                    self.u[i,j]*=(anoise0+dataOut.output_LP_integrated.real[0,i,0])/dataOut.output_LP_integrated.real[0,i,0]
                    self.u[i,j]*=(anoise0+dataOut.output_LP_integrated.real[0,j,0])/dataOut.output_LP_integrated.real[0,j,0]

                self.u[j,i]=self.u[i,j]

        #now error analyis for lag product matrix (diag), place in acf_err

        for i in range(dataOut.NACF):
            for j in range(dataOut.IBITS):
                if j==0:
                    dataOut.errors[0,i]=numpy.sqrt(self.u[i,i])
                else:
                    dataOut.errors[j,i]=numpy.sqrt(((dataOut.output_LP_integrated.real[0,i,0]+anoise0)*(dataOut.output_LP_integrated.real[0,i+j,0]+anoise0)+dataOut.output_LP_integrated.real[j,i,0]**2)/float(2*nis))

        print("Success")
        #print(dataOut.NRANGE)
        with suppress_stdout_stderr():
            pass
            #full_profile_profile.profile(numpy.transpose(dataOut.output_LP_integrated,(2,1,0)),numpy.transpose(dataOut.errors),self.powerb,dataOut.ne,dataOut.lags_LP,dataOut.thb,dataOut.bfm,dataOut.te,dataOut.ete,dataOut.ti,dataOut.eti,dataOut.ph,dataOut.eph,dataOut.phe,dataOut.ephe,dataOut.range1,dataOut.ut,dataOut.NACF,dataOut.fit_array_real,dataOut.status,dataOut.NRANGE,dataOut.IBITS)

        print("status: ",dataOut.status)

        if dataOut.status>=3.5:
            dataOut.te[:]=numpy.nan
            dataOut.ete[:]=numpy.nan
            dataOut.ti[:]=numpy.nan
            dataOut.eti[:]=numpy.nan
            dataOut.ph[:]=numpy.nan
            dataOut.eph[:]=numpy.nan
            dataOut.phe[:]=numpy.nan
            dataOut.ephe[:]=numpy.nan

        return dataOut

class LongPulseAnalysis_V2(Operation):
    """Operation to estimate ACFs, temperatures, total electron density and Hydrogen/Helium fractions from the Long Pulse data.

    Parameters:
    -----------
    NACF : int
        .*

    Example
    --------

    op = proc_unit.addOperation(name='LongPulseAnalysis', optype='other')
    op.addParameter(name='NACF', value='16', format='int')

    """

    def __init__(self, **kwargs):

        Operation.__init__(self, **kwargs)
        self.aux=1

    def run(self,dataOut,NACF):

        dataOut.NACF=NACF
        dataOut.heightList=dataOut.DH*(numpy.arange(dataOut.NACF))
        anoise0=dataOut.tnoise[0]
        anoise1=anoise0*0.0       #seems to be noise in 1st lag 0.015 before '14
        #print(anoise0)
        #exit(1)
        if self.aux:
            #dataOut.cut=31#26#height=31*15=465
            self.cal=numpy.zeros((dataOut.NLAG),'float32')
            self.drift=numpy.zeros((200),'float32')
            self.rdrift=numpy.zeros((200),'float32')
            self.ddrift=numpy.zeros((200),'float32')
            self.sigma=numpy.zeros((dataOut.NRANGE),order='F',dtype='float32')
            self.powera=numpy.zeros((dataOut.NRANGE),order='F',dtype='float32')
            self.powerb=numpy.zeros((dataOut.NRANGE),order='F',dtype='float32')
            self.perror=numpy.zeros((dataOut.NRANGE),order='F',dtype='float32')
            dataOut.ene=numpy.zeros((dataOut.NRANGE),'float32')
            self.dpulse=numpy.zeros((dataOut.NACF),'float32')
            self.lpulse=numpy.zeros((dataOut.NACF),'float32')
            dataOut.lags_LP=numpy.zeros((dataOut.IBITS),order='F',dtype='float32')
            self.lagp=numpy.zeros((dataOut.NACF),'float32')
            self.u=numpy.zeros((2*dataOut.NACF,2*dataOut.NACF),'float32')
            dataOut.ne=numpy.zeros((dataOut.NRANGE),order='F',dtype='float32')
            dataOut.te=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.ete=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.ti=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.eti=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.ph=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.eph=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.phe=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.ephe=numpy.zeros((dataOut.NACF),order='F',dtype='float32')
            dataOut.errors=numpy.zeros((dataOut.IBITS,max(dataOut.NRANGE,dataOut.NSHTS)),order='F',dtype='float32')
            dataOut.fit_array_real=numpy.zeros((max(dataOut.NRANGE,dataOut.NSHTS),dataOut.NLAG),order='F',dtype='float32')
            dataOut.status=numpy.zeros(1,'float32')
            dataOut.tx=240.0 #debería provenir del header #hybrid

            for i in range(dataOut.IBITS):
                dataOut.lags_LP[i]=float(i)*(dataOut.tx/150.0)/float(dataOut.IBITS) # (float)i*(header.tx/150.0)/(float)IBITS;

            self.aux=0

        dataOut.cut=30
        for i in range(30,15,-1):
            if numpy.nanmax(dataOut.acfs_error_to_plot[i,:])>=10 or dataOut.info2[i]==0:
                dataOut.cut=i-1

        for i in range(dataOut.NLAG):
            self.cal[i]=sum(dataOut.output_LP_integrated[i,:,3].real)

        #print(numpy.sum(self.cal)) #Coinciden
        #exit(1)
        self.cal/=float(dataOut.NRANGE)
        #print(anoise0)
        #print(anoise1)
        #exit(1)

        for j in range(dataOut.NACF+2*dataOut.IBITS+2):

            dataOut.output_LP_integrated.real[0,j,0]-=anoise0   #lag0 ch0
            dataOut.output_LP_integrated.real[1,j,0]-=anoise1   #lag1 ch0

            for i in range(1,dataOut.NLAG):  #remove cal data from certain lags
                 dataOut.output_LP_integrated.real[i,j,0]-=self.cal[i]
            k=max(j,26)   #constant power below range 26
            self.powera[j]=dataOut.output_LP_integrated.real[0,k,0]

            ## examine drifts here - based on 60 'indep.' estimates
        #print(numpy.sum(self.powera))
        #exit(1)
        #nis=dataOut.NSCAN*dataOut.NAVG*dataOut.nint*10
        nis = dataOut.nis
        #print("nis",nis)
        alpha=beta=delta=0.0
        nest=0
        gamma=3.0/(2.0*numpy.pi*dataOut.lags_LP[1]*1.0e-3)
        beta=gamma*(math.atan2(dataOut.output_LP_integrated.imag[14,0,2],dataOut.output_LP_integrated.real[14,0,2])-math.atan2(dataOut.output_LP_integrated.imag[1,0,2],dataOut.output_LP_integrated.real[1,0,2]))/13.0
        #print(gamma,beta)
        #exit(1)
        for i in range(1,3):
            gamma=3.0/(2.0*numpy.pi*dataOut.lags_LP[i]*1.0e-3)
            #print("gamma",gamma)
            for j in range(34,44):
                rho2=numpy.abs(dataOut.output_LP_integrated[i,j,0])/numpy.abs(dataOut.output_LP_integrated[0,j,0])
                dataOut.dphi2=(1.0/rho2-1.0)/(float(2*nis))
                dataOut.dphi2*=gamma**2
                pest=gamma*math.atan(dataOut.output_LP_integrated.imag[i,j,0]/dataOut.output_LP_integrated.real[i,j,0])
                #print("1",dataOut.output_LP_integrated.imag[i,j,0])
                #print("2",dataOut.output_LP_integrated.real[i,j,0])
                self.drift[nest]=pest
                self.ddrift[nest]=dataOut.dphi2
                self.rdrift[nest]=float(nest)
                nest+=1

        sorted(self.drift[:nest])

        #print(dataOut.dphi2)
        #exit(1)

        for j in range(int(nest/4),int(3*nest/4)):
            #i=int(self.rdrift[j])
            alpha+=self.drift[j]/self.ddrift[j]
            delta+=1.0/self.ddrift[j]

        alpha/=delta
        delta=1./numpy.sqrt(delta)
        vdrift=alpha-beta
        dvdrift=delta

        #need to develop estimate of complete density profile using all
        #available data

        #estimate sample variances for long-pulse power profile

        #nis=dataOut.NSCAN*dataOut.NAVG*dataOut.nint
        nis = dataOut.nis/10
        #print("nis",nis)

        self.sigma[:dataOut.NACF+2*dataOut.IBITS+2]=((anoise0+self.powera[:dataOut.NACF+2*dataOut.IBITS+2])**2)/float(nis)
        #print(self.sigma)
        #exit(1)
        ioff=1

        #deconvolve rectangular pulse shape from profile ==> powerb, perror


        ############# START nnlswrap#############

        if dataOut.ut_Faraday>14.0:
            alpha_nnlswrap=20.0
        else:
            alpha_nnlswrap=30.0

        range1_nnls=dataOut.NACF
        range2_nnls=dataOut.NACF+dataOut.IBITS-1

        g_nnlswrap=numpy.zeros((range1_nnls,range2_nnls),'float32')
        a_nnlswrap=numpy.zeros((range2_nnls,range2_nnls),'float64')

        for i in range(range1_nnls):
            for j in range(range2_nnls):
                if j>=i and j<i+dataOut.IBITS:
                    g_nnlswrap[i,j]=1.0
                else:
                    g_nnlswrap[i,j]=0.0

        a_nnlswrap[:]=numpy.matmul(numpy.transpose(g_nnlswrap),g_nnlswrap)

        numpy.fill_diagonal(a_nnlswrap,a_nnlswrap.diagonal()+alpha_nnlswrap**2)

                    #ERROR ANALYSIS#

        self.perror[:range2_nnls]=0.0
        self.perror[:range2_nnls]=numpy.matmul(1./(self.sigma[dataOut.IBITS+ioff:range1_nnls+dataOut.IBITS+ioff]),g_nnlswrap**2)
        self.perror[:range1_nnls]+=(alpha_nnlswrap**2)/(self.sigma[dataOut.IBITS+ioff:range1_nnls+dataOut.IBITS+ioff])
        self.perror[:range2_nnls]=1.00/self.perror[:range2_nnls]

        b_nnlswrap=numpy.zeros(range2_nnls,'float64')
        b_nnlswrap[:]=numpy.matmul(self.powera[dataOut.IBITS+ioff:range1_nnls+dataOut.IBITS+ioff],g_nnlswrap)

        x_nnlswrap=numpy.zeros(range2_nnls,'float64')
        x_nnlswrap[:]=nnls(a_nnlswrap,b_nnlswrap)[0]

        self.powerb[:range2_nnls]=x_nnlswrap
        #print(self.powerb[40])
        #print(self.powerb[66])
        #exit(1)
        #############END nnlswrap#############
        #print(numpy.sum(numpy.sqrt(self.perror[0:dataOut.NACF])))
        #print(self.powerb[0:dataOut.NACF])
        #exit(1)
        #estimate relative error for deconvolved profile (scaling irrelevant)
        #print(dataOut.NACF)
        dataOut.ene[0:dataOut.NACF]=numpy.sqrt(self.perror[0:dataOut.NACF])/self.powerb[0:dataOut.NACF]
        #print(numpy.sum(dataOut.ene))
        #exit(1)
        aux=0

        for i in range(dataOut.IBITS,dataOut.NACF):
            self.dpulse[i]=self.lpulse[i]=0.0
            for j in range(dataOut.IBITS):
                k=int(i-j)
                if k<36-aux and k>16:
                    self.dpulse[i]+=dataOut.ph2[k]/dataOut.h2[k]
                elif k>=36-aux:
                    self.lpulse[i]+=self.powerb[k]
            self.lagp[i]=self.powera[i]

        #find scale factor that best merges profiles

        qi=sum(self.dpulse[32:dataOut.NACF]**2/(self.lagp[32:dataOut.NACF]+anoise0)**2)
        ri=sum((self.dpulse[32:dataOut.NACF]*self.lpulse[32:dataOut.NACF])/(self.lagp[32:dataOut.NACF]+anoise0)**2)
        si=sum((self.dpulse[32:dataOut.NACF]*self.lagp[32:dataOut.NACF])/(self.lagp[32:dataOut.NACF]+anoise0)**2)
        ui=sum(self.lpulse[32:dataOut.NACF]**2/(self.lagp[32:dataOut.NACF]+anoise0)**2)
        vi=sum((self.lpulse[32:dataOut.NACF]*self.lagp[32:dataOut.NACF])/(self.lagp[32:dataOut.NACF]+anoise0)**2)

        alpha=(si*ui-vi*ri)/(qi*ui-ri*ri)
        beta=(qi*vi-ri*si)/(qi*ui-ri*ri)

        #form density profile estimate, merging rescaled power profiles
        #print(dataOut.h2)
        #print(numpy.sum(alpha))
        #print(numpy.sum(dataOut.ph2))
        self.powerb[16:36-aux]=alpha*dataOut.ph2[16:36-aux]/dataOut.h2[16:36-aux]
        self.powerb[36-aux:dataOut.NACF]*=beta

        #form Ne estimate, fill in error estimate at low altitudes

        dataOut.ene[0:36-aux]=dataOut.sdp2[0:36-aux]/dataOut.ph2[0:36-aux]
        dataOut.ne[:dataOut.NACF]=self.powerb[:dataOut.NACF]*dataOut.h2[:dataOut.NACF]/alpha
        #print(numpy.sum(self.powerb))
        #print(numpy.sum(dataOut.ene))
        #print(numpy.sum(dataOut.ne))
        #exit(1)
        #now do error propagation: store zero lag error covariance in u

        nis=dataOut.NSCAN*dataOut.NAVG*dataOut.nint/1   # DLH serious debris removal

        for i in range(dataOut.NACF):
            for j in range(i,dataOut.NACF):
                if j-i>=dataOut.IBITS:
                    self.u[i,j]=0.0
                else:
                    self.u[i,j]=dataOut.output_LP_integrated.real[j-i,i,0]**2/float(nis)
                    self.u[i,j]*=(anoise0+dataOut.output_LP_integrated.real[0,i,0])/dataOut.output_LP_integrated.real[0,i,0]
                    self.u[i,j]*=(anoise0+dataOut.output_LP_integrated.real[0,j,0])/dataOut.output_LP_integrated.real[0,j,0]

                self.u[j,i]=self.u[i,j]

        #now error analyis for lag product matrix (diag), place in acf_err

        for i in range(dataOut.NACF):
            for j in range(dataOut.IBITS):
                if j==0:
                    dataOut.errors[0,i]=numpy.sqrt(self.u[i,i])
                else:
                    dataOut.errors[j,i]=numpy.sqrt(((dataOut.output_LP_integrated.real[0,i,0]+anoise0)*(dataOut.output_LP_integrated.real[0,i+j,0]+anoise0)+dataOut.output_LP_integrated.real[j,i,0]**2)/float(2*nis))

        print("Success")
        with suppress_stdout_stderr():
            #pass
            full_profile_profile.profile(numpy.transpose(dataOut.output_LP_integrated,(2,1,0)),numpy.transpose(dataOut.errors),self.powerb,dataOut.ne,dataOut.lags_LP,dataOut.thb,dataOut.bfm,dataOut.te,dataOut.ete,dataOut.ti,dataOut.eti,dataOut.ph,dataOut.eph,dataOut.phe,dataOut.ephe,dataOut.range1,dataOut.ut,dataOut.NACF,dataOut.fit_array_real,dataOut.status,dataOut.NRANGE,dataOut.IBITS)

        if dataOut.status>=3.5:
            dataOut.te[:]=numpy.nan
            dataOut.ete[:]=numpy.nan
            dataOut.ti[:]=numpy.nan
            dataOut.eti[:]=numpy.nan
            dataOut.ph[:]=numpy.nan
            dataOut.eph[:]=numpy.nan
            dataOut.phe[:]=numpy.nan
            dataOut.ephe[:]=numpy.nan

        return dataOut

class PulsePairVoltage(Operation):
    '''
    Function PulsePair(Signal Power, Velocity)
    The real component of Lag[0] provides Intensity Information
    The imag component of Lag[1] Phase provides Velocity Information

    Configuration Parameters:
    nPRF = Number of Several PRF
    theta = Degree Azimuth angel Boundaries

    Input:
          self.dataOut
          lag[N]
    Affected:
          self.dataOut.spc
    '''
    isConfig       = False
    __profIndex    = 0
    __initime      = None
    __lastdatatime = None
    __buffer       = None
    noise          = None
    __dataReady    = False
    n              = None
    __nch          = 0
    __nHeis        = 0
    removeDC       = False
    ipp            = None
    lambda_        = 0

    def __init__(self,**kwargs):
        Operation.__init__(self,**kwargs)

    def setup(self, dataOut, n = None, removeDC=False):
        '''
        n= Numero de PRF's de entrada
        '''
        self.__initime        = None
        self.__lastdatatime   = 0
        self.__dataReady      = False
        self.__buffer         = 0
        self.__profIndex      = 0
        self.noise            = None
        self.__nch            = dataOut.nChannels
        self.__nHeis          = dataOut.nHeights
        self.removeDC         = removeDC
        self.lambda_          = 3.0e8/(9345.0e6)
        self.ippSec           = dataOut.ippSeconds
        self.nCohInt          = dataOut.nCohInt
        print("IPPseconds",dataOut.ippSeconds)

        print("ELVALOR DE n es:", n)
        if n == None:
            raise ValueError("n should be specified.")

        if n != None:
            if n<2:
                raise ValueError("n should be greater than 2")

        self.n       = n
        self.__nProf = n

        self.__buffer = numpy.zeros((dataOut.nChannels,
                                           n,
                                           dataOut.nHeights),
                                          dtype='complex')

    def putData(self,data):
        '''
        Add a profile to he __buffer and increase in one the __profiel Index
        '''
        self.__buffer[:,self.__profIndex,:]= data
        self.__profIndex      += 1
        return

    def pushData(self,dataOut):
        '''
        Return the PULSEPAIR and the profiles used in the operation
        Affected :  self.__profileIndex
        '''
        #----------------- Remove DC-----------------------------------
        if self.removeDC==True:
            mean    = numpy.mean(self.__buffer,1)
            tmp     = mean.reshape(self.__nch,1,self.__nHeis)
            dc= numpy.tile(tmp,[1,self.__nProf,1])
            self.__buffer = self.__buffer -  dc
        #------------------Calculo de Potencia ------------------------
        pair0       = self.__buffer*numpy.conj(self.__buffer)
        pair0       = pair0.real
        lag_0       = numpy.sum(pair0,1)
        #------------------Calculo de Ruido x canal--------------------
        self.noise  = numpy.zeros(self.__nch)
        for i in range(self.__nch):
            daux         = numpy.sort(pair0[i,:,:],axis= None)
            self.noise[i]=hildebrand_sekhon( daux ,self.nCohInt)

        self.noise       = self.noise.reshape(self.__nch,1)
        self.noise       = numpy.tile(self.noise,[1,self.__nHeis])
        noise_buffer     = self.noise.reshape(self.__nch,1,self.__nHeis)
        noise_buffer     = numpy.tile(noise_buffer,[1,self.__nProf,1])
        #------------------ Potencia recibida= P , Potencia senal = S , Ruido= N--
        #------------------   P= S+N  ,P=lag_0/N ---------------------------------
        #-------------------- Power --------------------------------------------------
        data_power       = lag_0/(self.n*self.nCohInt)
        #------------------  Senal  ---------------------------------------------------
        data_intensity   = pair0 - noise_buffer
        data_intensity   = numpy.sum(data_intensity,axis=1)*(self.n*self.nCohInt)#*self.nCohInt)
        #data_intensity   = (lag_0-self.noise*self.n)*(self.n*self.nCohInt)
        for i in range(self.__nch):
            for j in range(self.__nHeis):
                if data_intensity[i][j]  < 0:
                    data_intensity[i][j] = numpy.min(numpy.absolute(data_intensity[i][j]))

        #----------------- Calculo de Frecuencia y Velocidad doppler--------
        pair1            = self.__buffer[:,:-1,:]*numpy.conjugate(self.__buffer[:,1:,:])
        lag_1            = numpy.sum(pair1,1)
        data_freq        = (-1/(2.0*math.pi*self.ippSec*self.nCohInt))*numpy.angle(lag_1)
        data_velocity    = (self.lambda_/2.0)*data_freq

        #---------------- Potencia promedio estimada de la Senal-----------
        lag_0            = lag_0/self.n
        S                = lag_0-self.noise

        #---------------- Frecuencia Doppler promedio ---------------------
        lag_1            = lag_1/(self.n-1)
        R1               = numpy.abs(lag_1)

        #---------------- Calculo del SNR----------------------------------
        data_snrPP       = S/self.noise
        for i in range(self.__nch):
            for j in range(self.__nHeis):
                if data_snrPP[i][j]  < 1.e-20:
                    data_snrPP[i][j] = 1.e-20

        #----------------- Calculo del ancho espectral ----------------------
        L                = S/R1
        L                = numpy.where(L<0,1,L)
        L                = numpy.log(L)
        tmp              = numpy.sqrt(numpy.absolute(L))
        data_specwidth   = (self.lambda_/(2*math.sqrt(2)*math.pi*self.ippSec*self.nCohInt))*tmp*numpy.sign(L)
        n                = self.__profIndex

        self.__buffer    = numpy.zeros((self.__nch, self.__nProf,self.__nHeis),  dtype='complex')
        self.__profIndex = 0
        return data_power,data_intensity,data_velocity,data_snrPP,data_specwidth,n


    def pulsePairbyProfiles(self,dataOut):

        self.__dataReady     =  False
        data_power           =  None
        data_intensity       =  None
        data_velocity        =  None
        data_specwidth       =  None
        data_snrPP           =  None
        self.putData(data=dataOut.data)
        if self.__profIndex  == self.n:
            data_power,data_intensity, data_velocity,data_snrPP,data_specwidth, n   = self.pushData(dataOut=dataOut)
            self.__dataReady                   = True

        return data_power, data_intensity, data_velocity, data_snrPP, data_specwidth


    def pulsePairOp(self, dataOut, datatime= None):

        if self.__initime == None:
            self.__initime = datatime
        data_power, data_intensity, data_velocity, data_snrPP, data_specwidth = self.pulsePairbyProfiles(dataOut)
        self.__lastdatatime           = datatime

        if data_power is None:
            return None, None, None,None,None,None

        avgdatatime    = self.__initime
        deltatime      = datatime - self.__lastdatatime
        self.__initime = datatime

        return data_power, data_intensity, data_velocity, data_snrPP, data_specwidth, avgdatatime

    def run(self, dataOut,n = None,removeDC= False, overlapping= False,**kwargs):

        if not self.isConfig:
            self.setup(dataOut = dataOut, n    = n , removeDC=removeDC , **kwargs)
            self.isConfig   = True
        data_power, data_intensity, data_velocity,data_snrPP,data_specwidth, avgdatatime = self.pulsePairOp(dataOut, dataOut.utctime)
        dataOut.flagNoData                         = True

        if self.__dataReady:
            dataOut.nCohInt        *= self.n
            dataOut.dataPP_POW      = data_intensity # S
            dataOut.dataPP_POWER    = data_power     # P
            dataOut.dataPP_DOP      = data_velocity
            dataOut.dataPP_SNR      = data_snrPP
            dataOut.dataPP_WIDTH    = data_specwidth
            dataOut.PRFbyAngle      = self.n         #numero de PRF*cada angulo rotado que equivale a un tiempo.
            dataOut.utctime         = avgdatatime
            dataOut.flagNoData      = False
        return dataOut



# import collections
# from scipy.stats import mode
#
# class Synchronize(Operation):
#
#     isConfig = False
#     __profIndex = 0
#
#     def __init__(self, **kwargs):
#
#         Operation.__init__(self, **kwargs)
# #         self.isConfig = False
#         self.__powBuffer = None
#         self.__startIndex = 0
#         self.__pulseFound = False
#
#     def __findTxPulse(self, dataOut, channel=0, pulse_with = None):
#
#         #Read data
#
#         powerdB = dataOut.getPower(channel = channel)
#         noisedB = dataOut.getNoise(channel = channel)[0]
#
#         self.__powBuffer.extend(powerdB.flatten())
#
#         dataArray = numpy.array(self.__powBuffer)
#
#         filteredPower = numpy.correlate(dataArray, dataArray[0:self.__nSamples], "same")
#
#         maxValue = numpy.nanmax(filteredPower)
#
#         if maxValue < noisedB + 10:
#             #No se encuentra ningun pulso de transmision
#             return None
#
#         maxValuesIndex = numpy.where(filteredPower > maxValue - 0.1*abs(maxValue))[0]
#
#         if len(maxValuesIndex) < 2:
#             #Solo se encontro un solo pulso de transmision de un baudio, esperando por el siguiente TX
#             return None
#
#         phasedMaxValuesIndex = maxValuesIndex - self.__nSamples
#
#         #Seleccionar solo valores con un espaciamiento de nSamples
#         pulseIndex = numpy.intersect1d(maxValuesIndex, phasedMaxValuesIndex)
#
#         if len(pulseIndex) < 2:
#             #Solo se encontro un pulso de transmision con ancho mayor a 1
#             return None
#
#         spacing = pulseIndex[1:] - pulseIndex[:-1]
#
#         #remover senales que se distancien menos de 10 unidades o muestras
#         #(No deberian existir IPP menor a 10 unidades)
#
#         realIndex = numpy.where(spacing > 10 )[0]
#
#         if len(realIndex) < 2:
#             #Solo se encontro un pulso de transmision con ancho mayor a 1
#             return None
#
#         #Eliminar pulsos anchos (deja solo la diferencia entre IPPs)
#         realPulseIndex = pulseIndex[realIndex]
#
#         period = mode(realPulseIndex[1:] - realPulseIndex[:-1])[0][0]
#
#         print "IPP = %d samples" %period
#
#         self.__newNSamples = dataOut.nHeights #int(period)
#         self.__startIndex = int(realPulseIndex[0])
#
#         return 1
#
#
#     def setup(self, nSamples, nChannels, buffer_size = 4):
#
#         self.__powBuffer = collections.deque(numpy.zeros( buffer_size*nSamples,dtype=numpy.float),
#                                           maxlen = buffer_size*nSamples)
#
#         bufferList = []
#
#         for i in range(nChannels):
#             bufferByChannel = collections.deque(numpy.zeros( buffer_size*nSamples, dtype=numpy.complex) +  numpy.NAN,
#                                           maxlen = buffer_size*nSamples)
#
#             bufferList.append(bufferByChannel)
#
#         self.__nSamples = nSamples
#         self.__nChannels = nChannels
#         self.__bufferList = bufferList
#
#     def run(self, dataOut, channel = 0):
#
#         if not self.isConfig:
#             nSamples = dataOut.nHeights
#             nChannels = dataOut.nChannels
#             self.setup(nSamples, nChannels)
#             self.isConfig = True
#
#         #Append new data to internal buffer
#         for thisChannel in range(self.__nChannels):
#             bufferByChannel = self.__bufferList[thisChannel]
#             bufferByChannel.extend(dataOut.data[thisChannel])
#
#         if self.__pulseFound:
#             self.__startIndex -= self.__nSamples
#
#         #Finding Tx Pulse
#         if not self.__pulseFound:
#             indexFound = self.__findTxPulse(dataOut, channel)
#
#             if indexFound == None:
#                 dataOut.flagNoData = True
#                 return
#
#             self.__arrayBuffer = numpy.zeros((self.__nChannels, self.__newNSamples), dtype = numpy.complex)
#             self.__pulseFound = True
#             self.__startIndex = indexFound
#
#         #If pulse was found ...
#         for thisChannel in range(self.__nChannels):
#             bufferByChannel = self.__bufferList[thisChannel]
#             #print self.__startIndex
#             x = numpy.array(bufferByChannel)
#             self.__arrayBuffer[thisChannel] = x[self.__startIndex:self.__startIndex+self.__newNSamples]
#
#         deltaHeight = dataOut.heightList[1] - dataOut.heightList[0]
#         dataOut.heightList = numpy.arange(self.__newNSamples)*deltaHeight
# #             dataOut.ippSeconds = (self.__newNSamples / deltaHeight)/1e6
#
#         dataOut.data = self.__arrayBuffer
#
#         self.__startIndex += self.__newNSamples
#
#         return
