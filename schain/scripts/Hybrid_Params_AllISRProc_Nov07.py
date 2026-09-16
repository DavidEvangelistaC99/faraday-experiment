'''
#!/usr/bin/python3.7.4
'''

from schainpy.controller import Project

import json

#  dpath = '/home/cportilla/jars/ALL_ISR/hybrid'
# dpath = '/home/cportilla/Minotaur/2025_11/Faraday/main_radar/rawdata'
dpath = '/mnt/compartido2'
#dpath= '/media/cportilla/HDD/Faraday/Hybrid'
#dpath = '/home/cportilla/PC_DATA/Hybrid'

controller = Project()
controller.setup(id = '001',
                 name='Hybrid',
                 description='DP+LP')

figpath='/home/idi/Documents/DATA/All-ISR'
figpath_server=figpath
procpath = figpath
Show_plots = 1
Lag = '0'
navg = 16
nint_LP = 5 #6
nint_DP = int(navg*nint_LP)


read_unit = controller.addReadUnit(datatype='VoltageReader',
                                   path=dpath,
                                   startDate='2025/11/07',
                                   endDate='2025/11/07',
                                   startTime='00:00:00',
                                   endTime='23:59:59',
                                   online=0,
                                   getByBlock='True',
                                   nTries=20,#120,
                                   nTries_file =20,#9,
                                   walk=1,
                                   delay=30)

############################################# LP #############################################

proc_unit_LP = controller.addProcUnit(datatype='VoltageProc',inputId=read_unit.getId())
proc_unit_LP.addParameter(name='runNextUnit', value=True)

op1 = proc_unit_LP.addOperation(name='selectChannels')
op1.addParameter(name='channelList', value='0,1,2,3')

op121 = proc_unit_LP.addOperation(name='ProfileSelector')
op121.addParameter(name='profileRangeList', value='0,149')

op2 = proc_unit_LP.addOperation(name='filterByHeights')
op2.addParameter(name='window', value='10')

op1 = proc_unit_LP.addOperation(name='saturatedBlock')
op1.addParameter(name='limit', value='70')
op1.addParameter(name='mode', value='1')

op3 = proc_unit_LP.addOperation(name='CombineChannels_V2')
op3.addParameter(name='ch_list',   value='[[0,1],[0,1],[2],[3]]')
op3.addParameter(name='comb_list', value='sub,sum,none,none')

op3 = proc_unit_LP.addOperation(name='CrossProdHybrid')
op3.addParameter(name='nkill', value='12')

op4 = proc_unit_LP.addOperation(name='IntegrationHP')
op4.addParameter(name='nint', value=nint_LP)

############################################# DP #############################################
proc_unit_DP = controller.addProcUnit(datatype='VoltageProc',inputId=read_unit.getId())

op1 = proc_unit_DP.addOperation(name='selectChannels')
op1.addParameter(name='channelList', value='0,1')

op3 = proc_unit_DP.addOperation(name='CombineChannels_V2')
op3.addParameter(name='ch_list',   value='[[0,1],[0,1]]')
op3.addParameter(name='comb_list', value='sub,sum')

op121 = proc_unit_DP.addOperation(name='ProfileSelector')
op121.addParameter(name='profileRangeList', value='0,127')

op4 = proc_unit_DP.addOperation(name='deFlipHP')
op4.addParameter(name='channelList', value='(1,)')
op4.addParameter(name='byHeights', value='True')
op4.addParameter(name='HeiRangeList', value='(2670, 3340)') # 1068,1336 for dh 3.75 km # 2670, 3340 for dh 1.5km

op2 = proc_unit_DP.addOperation(name='filterByHeights')
op2.addParameter(name='window', value='10')

op5 = proc_unit_DP.addOperation(name='LagsReshapeDP_V2')

proc_spectra_DP = controller.addProcUnit(datatype='SpectraLagProc',inputId=proc_unit_DP.getId())
proc_spectra_DP.addParameter(name='nFFTPoints', value='32')
proc_spectra_DP.addParameter(name='nProfiles', value='32')
proc_spectra_DP.addParameter(name='ByLags', value=True)
proc_spectra_DP.addParameter(name='nLags', value='11')
proc_spectra_DP.addParameter(name='LagPlot', value=Lag)

op44 = proc_spectra_DP.addOperation(name='selectHeights')
op44.addParameter(name='minIndex', value='0')
op44.addParameter(name='maxIndex', value='66') #66

op08 = proc_spectra_DP.addOperation(name='IntegrationFaradaySpectra2')
op08.addParameter(name='n', value=nint_DP) #16*6 = navg*nint

op6 = proc_spectra_DP.addOperation(name='removeDCLagFlip')

merge = controller.addProcUnit(datatype='MergeProc',inputId=[proc_unit_LP.getId(),proc_spectra_DP.getId()])

merge.addParameter(name='attr_data', value='dataLag_spc') #DP
merge.addParameter(name='attr_data_2', value='dataLag_cspc') #DP
merge.addParameter(name='attr_data_3', value='output_LP_integrated') #LP
merge.addParameter(name='mode', value='6')

op10 = merge.addOperation(name='SpcVoltageDataToHybrid')

op6 = merge.addOperation(name='FlagBadHeightsSpectra')

#op11 = merge.addOperation(name='CleanCohEchoes')

op8 = merge.addOperation(name='DoublePulseACFs')

op9 = merge.addOperation(name='FaradayAngleAndDPPower')

op10 = merge.addOperation(name='IGRFModel')

op11 = merge.addOperation(name='ElectronDensityFaraday')
op11.addParameter(name='NSHTS', value='41')
op11.addParameter(name='RATE', value='1.8978873e-6')

op12 = merge.addOperation(name='NormalizeDPPowerRoberto_V2')

op012 = merge.addOperation(name='DPTemperaturesEstimation')
op012.addParameter(name='IBITS', value='16')

op185 = merge.addOperation(name='DenCorrection', optype='other')
op185.addParameter(name='savecf', value='1')

op022 = merge.addOperation(name='ACFs')

op0122 = merge.addOperation(name='LongPulseAnalysis')
op0122.addParameter(name='NACF', value='90')

op026 = merge.addOperation(name='ACFsPlot')
op026.addParameter(name='id', value='177')
op026.addParameter(name='wintitle', value='ACFs DP')
op026.addParameter(name='ymin', value='180')
op026.addParameter(name='ymax', value='600')
op026.addParameter(name='xmax', value='2')
op026.addParameter(name='save', value=figpath_server)

op026 = merge.addOperation(name='ACFsLPPlot')
op026.addParameter(name='id', value='177')
op026.addParameter(name='wintitle', value='ACFs LP')
op026.addParameter(name='ymin', value='550')
op026.addParameter(name='save', value=figpath_server)

op026 = merge.addOperation(name='EDensityHPPlot')
op026.addParameter(name='id', value='179')
op026.addParameter(name='wintitle', value='Electron Density')
op026.addParameter(name='ymin', value='150')
op026.addParameter(name='ymax', value='1350')
op026.addParameter(name='xmin', value='1e0')
op026.addParameter(name='xmax', value='1e9')
op026.addParameter(name='save', value=figpath_server)

op0166 = merge.addOperation(name='TempsHPPlot')
op0166.addParameter(name='id', value='1755')
op0166.addParameter(name='wintitle', value='Temperatures')
op0166.addParameter(name='ymin', value='160')
op0166.addParameter(name='save', value=figpath_server)

op0177 = merge.addOperation(name='FracsHPPlot')
op0177.addParameter(name='id', value='1700')
op0177.addParameter(name='wintitle', value='Composition')
op0177.addParameter(name='ymin', value='160')
op0177.addParameter(name='save', value=figpath_server)

op018 = merge.addOperation(name='DataSaveCleanerHP')

op181 = merge.addOperation(name='ETempRTIPlot')
op181.addParameter(name='id', value='176')
op181.addParameter(name='wintitle', value='ElectronTemp')
op181.addParameter(name='xmin', value='0')
op181.addParameter(name='xmax', value='24')
op181.addParameter(name='ymin', value='180')
op181.addParameter(name='ymax', value='1350')
op181.addParameter(name='zmin', value='500')
op181.addParameter(name='zmax', value='6500')
op181.addParameter(name='xrange', value=str(24*1))
op181.addParameter(name='save', value=figpath_server)

op182 = merge.addOperation(name='ITempRTIPlot')
op182.addParameter(name='id', value='176')
op182.addParameter(name='wintitle', value='IonTemp')
op182.addParameter(name='xmin', value='0')
op182.addParameter(name='xmax', value='24')
op182.addParameter(name='ymin', value='180')
op182.addParameter(name='ymax', value='1350')
op182.addParameter(name='zmin', value='500')
op182.addParameter(name='zmax', value='6500')
op182.addParameter(name='xrange', value=str(24*1))
op182.addParameter(name='save', value=figpath_server)

op183 = merge.addOperation(name='HFracRTIPlot')
op183.addParameter(name='id', value='176')
op183.addParameter(name='wintitle', value='H+ Frac')
op183.addParameter(name='xmin', value='0')
op183.addParameter(name='xmax', value='24')
op183.addParameter(name='ymin', value='180')
op183.addParameter(name='ymax', value='1350')
op183.addParameter(name='zmin', value='0')
op183.addParameter(name='zmax', value='1')
op183.addParameter(name='xrange', value=str(24*1))
op183.addParameter(name='save', value=figpath_server)

op184 = merge.addOperation(name='HeFracRTIPlot')
op184.addParameter(name='id', value='176')
op184.addParameter(name='wintitle', value='He+ Frac')
op184.addParameter(name='xmin', value='0')
op184.addParameter(name='xmax', value='24')
op184.addParameter(name='ymin', value='180')
op184.addParameter(name='ymax', value='1350')
op184.addParameter(name='zmin', value='0')
op184.addParameter(name='zmax', value='0.5')
op184.addParameter(name='xrange', value=str(24*1))
op184.addParameter(name='save', value=figpath_server)

op16 = merge.addOperation(name='DenRTIPlot')
op16.addParameter(name='id', value='174')
op16.addParameter(name='wintitle', value='Electron Density RTI')
op16.addParameter(name='xmin', value='0')
op16.addParameter(name='xmax', value='24')
op16.addParameter(name='ymin', value='180')
op16.addParameter(name='ymax', value='1350')
op16.addParameter(name='zmin', value='1000')
op16.addParameter(name='zmax', value='10000000')
op16.addParameter(name='xrange', value=str(24*1))
op16.addParameter(name='save', value=figpath_server)

# HDF5
#dataList = ['utctime','utctimeInit','data_spc','data_snr','data_pow','data_dop','data_width']
#metadataList = ['nIncohInt','nCohInt','heightList','timeZone','normFactor','abscissaList']
'''
dataList = ['utctime','utctimeInit','ElecTempFinal','IonTempFinal','PhyFinal','PheFinal','DensityFinal']
metadataList = ['nIncohInt','nCohInt','heightList','timeZone','normFactor','abscissaList']


writer = merge.addOperation(name='HDFWriter')
writer.addParameter(name='path', value=procpath)
writer.addParameter(name='blocksPerFile', value=1000)
writer.addParameter(name='dataList', value=dataList, format='list')
writer.addParameter(name='metadataList', value=metadataList, format='list')
'''
one = {'gdlatr': 'lat', 'gdlonr': 'lon', 'inttms': 'paramInterval',
        'azm': 'azm', 'elm': 'elm', 'pl':'pl', 'rgate':'DH', 'tfreq': 'tfreq'} #reader gdlatr-->lat only 1D
two = {
    'gdalt': 'heightList',   #<----- nmonics
    #'VIPE': ('data_output', 0),
    # 'DVIPE': ('data_output', 1),
    'NE': ('DensityFinal', 0),
    'DNE': ('EDensityFinal', 0),
    'TE': ('ElecTempFinal', 0),
    'DTE': ('EElecTempFinal', 0),
    'TI': ('IonTempFinal', 0),
    'DTI': ('EIonTempFinal', 0),
    'PH+': ('PhyFinal', 0),
    'DPH+': ('EPhyFinal', 0),
    'PHE+': ('PheFinal', 0),
    'DPHE+': ('EPheFinal', 0)



    # 'DVIPN': ('data_output', 3),
    #'snl': 'data_SNR'
    } #writer

#f=open('/home/cportilla/proyecto/Faraday/moder_test.txt','r')
#file_contents=f.read()


ind = ['gdalt']

meta = {
    'kinst': 10, #instrumnet code
    'kindat': 1802, #type of data
    'catalog': {
        'principleInvestigator': 'Danny Scipión',
        'expPurpose': 'Electron Density'#,
        #'sciRemarks': file_contents
        },
    'header': {
        'analyst': 'David Hysell'
    }
}
#f.close()

'''
op_writer = merge.addOperation(name='MADWriter')
op_writer.addParameter(name='path', value=figpath)
op_writer.addParameter(name='format', value='hdf5')
op_writer.addParameter(name='oneDDict', value=json.dumps(one))
op_writer.addParameter(name='twoDDict', value=json.dumps(two))
op_writer.addParameter(name='ind2DList', value=json.dumps(ind))
op_writer.addParameter(name='metadata', value=json.dumps(meta))
#op_writer.addParameter(name='blocks', value=60)
'''

controller.start()
