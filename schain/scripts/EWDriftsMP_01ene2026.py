
import os, sys
import json

#from controller import *
from schainpy.controller import Project

desc = "EW DRIFTS MP Experiment"
filename = "EWDrifts.xml"

controllerObj = Project()

controllerObj.setup(id = '191', name='test01', description=desc)

#Experimentos

#path = "/data/dia"
#path = '/home/pcondor/data'
#path = '/media/pcondor/DATA1/Database/ewdriftsschain2023prue/data'
path = '/home/david/Documentos/DATA/DRIFTS/rawdata'
#path = '/data/isr_julia'
#pathFigure = '/media/pcondor/DATA1/Database/ewdriftsschain2023wh5'
#pathFile = '/media/pcondor/DATA1/Database/ewdriftsdic2024MPsch/EW_Drifts_01dic'
pathFile = '/home/david/Documentos/DATA/DRIFTS/rawdata/Drifts'
pathFigure = pathFile
pathFileavg = pathFile+'/avg'
pathFiledata = pathFile+'/Drifts-data'

xmin = 0
xmax = 24
#------------------------------------------------------------------------------------------------
readUnitConfObj = controllerObj.addReadUnit(datatype='VoltageReader',
                                            path=path,
                                            startDate='2026/03/29',
                                            endDate=  '2026/03/29',
                                            startTime='03:30:00',
                                            endTime=  '03:40:59',
                                            online=0,
                                            getByBlock=1,
                                            walk=1)

#--------------------------------------------------------------------------------------------------

procUnitConfObj0 = controllerObj.addProcUnit(datatype='VoltageProc', inputId=readUnitConfObj.getId())

#opObj11 = procUnitConfObj0.addOperation(name='selectHeights')
# # opObj11.addParameter(name='minHei', value='320.0', format='float')
# # opObj11.addParameter(name='maxHei', value='350.0', format='float')
#opObj11.addParameter(name='minHei', value='0.01', format='float')
#opObj11.addParameter(name='maxHei', value='960.0', format='float')

opObj11 = procUnitConfObj0.addOperation(name='selectChannels')
opObj11.addParameter(name='channelList', value='0,0,1,1', format='intlist')

#opObj11 = procUnitConfObj0.addOperation(name='Reshaper')
#opObj11.addParameter(name='shape', value='(500,980)', format='intlist')

opObj11 = procUnitConfObj0.addOperation(name='ProfileSelector', optype='other')
opObj11.addParameter(name='profileRangeList', value='0,127', format='intlist')
 
opObj11 = procUnitConfObj0.addOperation(name='filterByHeights')
opObj11.addParameter(name='window', value='10', format='int')

#code=[[1,1,-1],[-1,-1,1]]  #hasta 12 nov
code = [[1,1,-1] ,[-1,-1,1],[1,1,-1],[1,1,-1] ,[1,1,-1] ,[1,1,-1],[-1,-1,1] ,[-1,-1,1] ,[-1,-1,1],[-1,-1,1] ,[-1,-1,1] ,[-1,-1,1],[-1,-1,1] ,[-1,-1,1] ,[1,1,-1],[1,1,-1] ,[1,1,-1] ,[1,1,-1],[1,1,-1] ,[-1,-1,1] ,[1,1,-1],[-1,-1,1] ,[-1,-1,1] ,[1,1,-1],[1,1,-1] ,[-1,-1,1] ,[1,1,-1],[-1,-1,1] ,[-1,-1,1] ,[-1,-1,1],[-1,-1,1] ,[1,1,-1] ,[-1,-1,1],[-1,-1,1],[-1,-1,1],[1,1,-1],[-1,-1,1],[1,1,-1],[1,1,-1],[1,1,-1],[-1,-1,1],[1,1,-1],[1,1,-1],[1,1,-1],[1,1,-1],[1,1,-1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[1,1,-1],[-1,-1,1],[1,1,-1],[1,1,-1],[1,1,-1],[1,1,-1],[1,1,-1],[1,1,-1],[-1,-1,1],[1,1,-1],[-1,-1,1],[1,1,-1],[-1,-1,1],[-1,-1,1],[1,1,-1],[-1,-1,1],[1,1,-1],[1,1,-1],[-1,-1,1],[-1,-1,1],[1,1,-1],[-1,-1,1],[1,1,-1],[1,1,-1],[1,1,-1],[-1,-1,1],[1,1,-1],[1,1,-1],[-1,-1,1],[1,1,-1],[-1,-1,1],[-1,-1,1],[1,1,-1],[1,1,-1],[-1,-1,1],[1,1,-1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[1,1,-1],[-1,-1,1],[1,1,-1],[1,1,-1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[1,1,-1],[1,1,-1],[1,1,-1],[-1,-1,1],[1,1,-1],[1,1,-1],[-1,-1,1],[-1,-1,1],[1,1,-1],[1,1,-1],[1,1,-1],[1,1,-1],[-1,-1,1],[1,1,-1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[1,1,-1],[1,1,-1],[1,1,-1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[1,1,-1],[1,1,-1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[-1,-1,1]]
#code = [[1,1,-1],[-1,-1,1],[1,1,-1],[1,1,-1],[1,1,-1],[-1,-1,1],[-1,-1,1],[1,1,-1],[1,1,-1],[-1,-1,1],[1,1,-1],[1,1,-1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[1,1,-1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[1,1,-1],[1,1,-1],[-1,-1,1],[-1,-1,1],[1,1,-1],[-1,-1,1],[1,1,-1],[1,1,-1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[1,1,-1],[-1,-1,1],[1,1,-1],[1,1,-1],[-1,-1,1],[-1,-1,1],[1,1,-1],[1,1,-1],[-1,-1,1],[1,1,-1],[-1,-1,1],[-1,-1,1],[1,1,-1],[1,1,-1],[-1,-1,1],[-1,-1,1],[1,1,-1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[1,1,-1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[1,1,-1],[1,1,-1],[1,1,-1],[1,1,-1],[1,1,-1],[1,1,-1],[1,1,-1],[1,1,-1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[1,1,-1],[-1,-1,1],[1,1,-1],[1,1,-1],[-1,-1,1],[-1,-1,1],[1,1,-1],[-1,-1,1],[1,1,-1],[1,1,-1],[1,1,-1],[1,1,-1],[-1,-1,1],[1,1,-1],[1,1,-1],[1,1,-1],[-1,-1,1],[1,1,-1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[1,1,-1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[1,1,-1],[1,1,-1],[1,1,-1],[1,1,-1],[1,1,-1],[-1,-1,1],[1,1,-1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[1,1,-1],[-1,-1,1],[1,1,-1],[-1,-1,1],[1,1,-1],[1,1,-1],[-1,-1,1],[1,1,-1],[-1,-1,1],[-1,-1,1],[1,1,-1],[1,1,-1],[-1,-1,1],[1,1,-1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[1,1,-1],[-1,-1,1],[-1,-1,1],[1,1,-1],[-1,-1,1],[1,1,-1],[1,1,-1],[-1,-1,1],[-1,-1,1],[1,1,-1],[-1,-1,1],[1,1,-1],[1,1,-1],[1,1,-1],[-1,-1,1],[1,1,-1],[-1,-1,1],[-1,-1,1],[1,1,-1],[1,1,-1],[1,1,-1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[1,1,-1],[-1,-1,1],[-1,-1,1],[1,1,-1],[1,1,-1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[1,1,-1],[-1,-1,1],[1,1,-1],[1,1,-1],[1,1,-1],[1,1,-1],[1,1,-1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[1,1,-1],[1,1,-1],[1,1,-1],[-1,-1,1],[-1,-1,1],[1,1,-1],[1,1,-1],[1,1,-1],[1,1,-1],[1,1,-1],[-1,-1,1],[1,1,-1],[1,1,-1],[-1,-1,1],[1,1,-1],[1,1,-1],[-1,-1,1],[-1,-1,1],[-1,-1,1],[1,1,-1],[1,1,-1],[1,1,-1],[1,1,-1],[1,1,-1],[-1,-1,1],[-1,-1,1],[1,1,-1],[-1,-1,1],[-1,-1,1],[1,1,-1],[1,1,-1]]
opObj11 = procUnitConfObj0.addOperation(name='Decoder', optype='other')
opObj11.addParameter(name='code', value=code, format='floatlist')
opObj11.addParameter(name='nCode', value='128', format='int')
opObj11.addParameter(name='nBaud', value='3', format='int')

opObj11 = procUnitConfObj0.addOperation(name='selectHeights')
opObj11.addParameter(name='minHei', value='0.0', format='float')
opObj11.addParameter(name='maxHei', value='960', format='float')

procUnitConfObj1 = controllerObj.addProcUnit(datatype='SpectraProc', inputId=procUnitConfObj0.getId())
procUnitConfObj1.addParameter(name='nFFTPoints', value='128', format='int')
procUnitConfObj1.addParameter(name='nProfiles', value='128', format='int')
#procUnitConfObj1.addParameter(name='pairsList', value='(2,3),(4,5)', format='pairsList')#,(2,3)
procUnitConfObj1.addParameter(name='pairsList', value='(0,1),(2,3)', format='pairsList')

#opObj11 = procUnitConfObj1.addOperation(name='selectHeights')
# # opObj11.addParameter(name='minHei', value='320.0', format='float') 
# # opObj11.addParameter(name='maxHei', value='350.0', format='float') 
#opObj11.addParameter(name='minHei', value='0.0', format='float') 
#opObj11.addParameter(name='maxHei', value='960.0', format='float')

#opObj11 = procUnitConfObj1.addOperation(name='selectChannels')
#opObj11.addParameter(name='channelList', value='2,3,4,5', format='intlist') 
   
opObj11 = procUnitConfObj1.addOperation(name='IncohInt', optype='other')
opObj11.addParameter(name='n', value='1', format='float')
#opObj11.addParameter(name='timeInterval', value='300.0', format='float')

#opObj13 = procUnitConfObj1.addOperation(name='removeDC')
#opObj11 = procUnitConfObj1.addOperation(name='removeInterference')

#opObj14 = procUnitConfObj1.addOperation(name='SpectraPlot', optype='other')
#opObj14.addParameter(name='id', value='65', format='int')
## # opObj14.addParameter(name='wintitle', value='Con interf', format='str')
#opObj14.addParameter(name='save', value=pathFigure, format='str')
##opObj14.addParameter(name='save_period', value=1, format='int')
#opObj14.addParameter(name='zmin', value='10', format='int')
#opObj14.addParameter(name='zmax', value='26', format='int')
# 

#opObj12 = procUnitConfObj1.addOperation(name='RTIPlot', optype='other')
#opObj12.addParameter(name='id', value='63', format='int')
#opObj12.addParameter(name='wintitle', value='RTI Plot', format='str')
#opObj12.addParameter(name='save', value=pathFigure, format='str')
#opObj12.addParameter(name='save_period', value=3600, format='int')
##opObj12.addParameter(name='figpath', value = pathFigure, format='str')
#opObj12.addParameter(name='xmin', value=xmin, format='float')
#opObj12.addParameter(name='xmax', value=xmax, format='float')
#opObj12.addParameter(name='zmin', value='20', format='int')
#opObj12.addParameter(name='zmax', value='36', format='int')

#--------------------------------------------------------------------------------------------------

procUnitConfObj2 = controllerObj.addProcUnit(datatype='ParametersProc', inputId=procUnitConfObj1.getId())
opObj20 = procUnitConfObj2.addOperation(name='SpectralFitting', optype='other')
opObj20.addParameter(name='path', value='/home/pcondor/DIR_MADRIGAL/schain/schainpy/model/proc', format='str')
opObj20.addParameter(name='file', value='modelSpectralFitting', format='str')
opObj20.addParameter(name='groupList', value='(0,1),(2,3)',format='multiList')
opObj20.addParameter(name='taver', value='5')
opObj20.addParameter(name='coh_th', value='[1]',format='multiList')
opObj20.addParameter(name='hei_th', value='[2000]',format='multiList')
opObj20.addParameter(name='filec', value='weightfit2025', format='str')
opObj20.addParameter(name='snr_coh', value='7')
opObj20.addParameter(name='snr_th', value='-12')

opObj22 = procUnitConfObj2.addOperation(name='HDFWriter', optype='other')
opObj22.addParameter(name='path', value=pathFiledata)
opObj22.addParameter(name='blocksPerFile', value='1')
opObj22.addParameter(name='metadataList',value='heightList,timeZone')
opObj22.addParameter(name='dataList',value='tmp_spectra_i,tmp_cspectra_i,tmp_spectra_c,tmp_cspectra_c,clean_num_aver,coh_num_aver,sat_spectra,sat_cspectra,index,utctime,nIncohInt,nCohInt,nProfiles,nFFTPoints,ippFactor,ippSeconds,paramInterval')
##opObj22.addParameter(name='dataList',value='tmp_spectra_i,tmp_cspectra_i,tmp_spectra_c,tmp_cspectra_c,clean_num_aver,coh_num_aver,index,utctime,nIncohInt,nCohInt,nProfiles,nFFTPoints,normFactor,channelList,ippFactor,ippSeconds')

#angles :-2.41116      3.01082
#angles = [ -2.40286, 3.00663] + 0.104139
#beam_pos= [[-11704, 8612],[13207, 8757]]
# user_angles = 1 & angles = [ -2.22796, 3.60630] + 0.116634 [-2.11133, 3.72293] desde 20-2-2025 hasta 23-6-25 -120.838 86.1597 124.463 87.1177
# user_angles = 1 & angles = [ -2.22796, 3.60630] + 0.129126 [-2.09883, 3.73543] desde 24-6-2025 hasta 08-04-26
# user_angles = 1 & angles = [ -2.22796, 3.60630] + 0.141619 [-2.08591, 3.74792] desde 09-4-2026 hasta 16-11-26
# user_angles = 1 & angles = [ -2.22796, 3.61879] + 0.141619 [-2.08591, 3.76041] desde 17-11-2026 hasta 31-12-26
opObj21 = procUnitConfObj2.addOperation(name='EWDriftsEstimation', optype='other')
opObj21.addParameter(name='zenith', value='-2.09883, 3.73543', format='floatlist')
opObj21.addParameter(name='zenithCorrection', value='0.0', format='float')
opObj21.addParameter(name='beam_pos', value='-120.83, 86.16, 124.46, 87.12',format='floatlist')
opObj21.addParameter(name='fileDrifts', value=pathFile)

# Drifts en h5
one = {'gdlatr': 'lat', 'gdlonr': 'lon', 'spcst':'spcst','pl':'pl','cbadn':'cbadn','inttms': 'inttms','azdir7':'azw','eldir7':'elw','azdir8':'aze','eldir8':'ele','jro14':'jro14','jro15':'jro15','jro16':'jro16','nwlos':'nwlos'}
two = {
    'range': ('params', 0),
    'gdalt': ('params', 1),
    'VIPN': ('params', 2),
    'dvipn': ('params', 3),
    'vipe': ('params', 4),
    'dvipe': ('params', 5),
    'vi7': ('params', 6),
    'dvi7': ('params', 7),
    'vi8': ('params', 8),
    'dvi8': ('params', 9),
    'PAIWL': ('params', 10),
    'pacwl': ('params', 11),
    'pbiwl': ('params', 12),
    'pbcwl': ('params', 13),
    'pciel': ('params', 14),
    'pccel': ('params', 15),
    'pdiel': ('params', 16),
    'pdcel': ('params', 17),
    'jro10': ('params', 18),
    'jro11': ('params', 19)
    } #writer
ind = ['gdalt']

meta = {
    'kinst': 12, #instrument code
    'kindat': 1910, #type of data
    'catalog': {
        'principleInvestigator': 'Danny Scipión',
        'expPurpose': 'Drifts version 2 Improved calculation and cleaning'#,
        #'sciRemarks': file_contents
        },
    'header': {
        'analyst': 'Danny Scipión'
    }
}

op_writer = procUnitConfObj2.addOperation(name='MADWriter')
op_writer.addParameter(name='path', value=pathFile)
op_writer.addParameter(name='format', value='hdf5')
op_writer.addParameter(name='oneDDict', value=json.dumps(one))
op_writer.addParameter(name='twoDDict', value=json.dumps(two))
op_writer.addParameter(name='ind2DList', value=json.dumps(ind))
op_writer.addParameter(name='metadata', value=json.dumps(meta))

op_writer = procUnitConfObj2.addOperation(name='setHeightDriftsavg')

# Avg Drifts
one_avg = {'gdlatr': 'lat', 'gdlonr': 'lon', 'spcst':'spcst','pl':'pl','cbadn':'cbadn','inttms': 'inttms'}
two_avg = {
    'range': ('params_avg', 4),
    'gdalt': ('params_avg', 5),
    'altav': ('params_avg', 6),
    'VIPN': ('params_avg', 0),
    'dvipn': ('params_avg', 1),
    'vipe': ('params_avg', 2),
    'dvipe': ('params_avg', 3)
    }
ind_avg = ['gdalt']
meta = {
    'kinst': 12, #instrument code
    'kindat': 1911, #type of data
    'catalog': {
        'principleInvestigator': 'Danny Scipión',
        'expPurpose': 'Drifts version 2 Improved calculation and cleaning'#,
        #'sciRemarks': file_contents
        },
    'header': {
        'analyst': 'Danny Scipión'
    }
}
#dataOut.heightList = dataOut.params_avg[4]
op_writer = procUnitConfObj2.addOperation(name='MADWriter')
op_writer.addParameter(name='path', value=pathFileavg)
op_writer.addParameter(name='format', value='hdf5')
op_writer.addParameter(name='oneDDict', value=json.dumps(one_avg))
op_writer.addParameter(name='twoDDict', value=json.dumps(two_avg))
op_writer.addParameter(name='ind2DList', value=json.dumps(ind_avg))
op_writer.addParameter(name='metadata', value=json.dumps(meta))

op_writer = procUnitConfObj2.addOperation(name='setHeightDrifts')

opObj24 = procUnitConfObj2.addOperation(name='SpectralMomentsPlot', optype='other')
opObj24.addParameter(name='id', value='1', format='int')
### # opObj14.addParameter(name='wintitle', value='Spectral Averaged', format='str')
opObj24.addParameter(name='save', value=pathFigure, format='str')
###opObj24.addParameter(name='save_period', value=1, format='int')
opObj24.addParameter(name='zmin', value='-9', format='int')
opObj24.addParameter(name='zmax', value='19', format='int')
opObj24.addParameter(name='xaxis', value='Velocity', format='str')

# 
titles=('SNR+1,Vertical Drifts,Zonal Drifts')
#titles=('Zonal Drifts,Vertical Drifts')
opObj23 = procUnitConfObj2.addOperation(name='GenericRTIPlot')
opObj23.addParameter(name='colormaps', value='jro,seismic,seismic')
opObj23.addParameter(name='attr_data', value='data_snr1,data_output')
#opObj23.addParameter(name='colormaps', value='RdBu,RdBu')
#opObj23.addParameter(name='attr_data', value='data_output')
opObj23.addParameter(name='wintitle', value='EW Drifts')
opObj23.addParameter(name='save', value=pathFigure)
opObj23.addParameter(name='titles', value=titles)
opObj23.addParameter(name='zfactors', value='1,1,1')
opObj23.addParameter(name='zlimits', value='(-4,16),(-50,50),(-300,300)')
opObj23.addParameter(name='cb_labels', value='dB,m/s,m/s')
opObj23.addParameter(name='throttle', value='1')
opObj23.addParameter(name='xmin', value=xmin)
opObj23.addParameter(name='xmax', value=xmax)
#opObj23.addParameter(name='exp_code', value='110', format='int')
#opObj23.addParameter(name='server', value='10.10.110.243:4444', format='int')
#opObj23.addParameter(name='tag', value= 'jicamarca', format='str')

#--------------------------------------------------------------------------------------------------

controllerObj.start()
