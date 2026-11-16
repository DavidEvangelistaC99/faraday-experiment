'''
#!/usr/bin/python3.7.4
'''

from schainpy.controller import Project
import json

controller = Project()
controller.setup(id = '001',
                 name='Faraday',
                 description='DP')

figpath = '/home/david/Documents/faraday-experiment/results/Faraday/Faraday_2026_05/test/'
figpath_server=figpath
data_path = '/mnt/compartido3'
startDate = '2026/05/25'
endDate = '2026/05/26'
startTime='11:00:00'
endTime='23:59:59'
Lag = '0'
db_range = ['15','25']
db_range = ['13','21']

read_unit = controller.addReadUnit(datatype='VoltageReader',
                                   path=data_path,
                                   startDate=startDate,
                                   endDate=endDate,
                                   startTime=startTime,
                                   endTime=endTime,
                                   online=0,
                                   getByBlock='True',
                                   nTries=20,#120,
                                   nTries_file =20,#9,
                                   walk=1,
                                   delay=30)

proc_voltage = controller.addProcUnit(datatype='VoltageProc',inputId=read_unit.getId())

op3 = proc_voltage.addOperation(name='CombineChannels')
op3.addParameter(name='sub_list', value='[[0,1]]')
op3.addParameter(name='sum_list', value='[[0,1]]')

#op42 = proc_voltage.addOperation(name='selectChannels')
#op42.addParameter(name='channelList', value='(0,1)')

op2 = proc_voltage.addOperation(name='ProfileSelector')
op2.addParameter(name='profileRangeList', value='0,139')

op4 = proc_voltage.addOperation(name='deFlip')
op4.addParameter(name='channelList', value='1')

op3 = proc_voltage.addOperation(name='filterByHeights')
op3.addParameter(name='window', value='10') #IPP 1500 km

op5 = proc_voltage.addOperation(name='LagsReshape')
#op5.addParameter(name='NSCAN', value='198')

proc_spectra = controller.addProcUnit(datatype='SpectraLagProc',inputId=proc_voltage.getId())
proc_spectra.addParameter(name='nFFTPoints', value='12')
proc_spectra.addParameter(name='nProfiles', value='12')
proc_spectra.addParameter(name='ByLags', value='True')
proc_spectra.addParameter(name='nLags', value='11')
proc_spectra.addParameter(name='LagPlot', value=Lag)

opObj11 = proc_spectra.addOperation(name='NoisePlot')
opObj11.addParameter(name='id', value='3')
opObj11.addParameter(name='wintitle', value='Noise')
opObj11.addParameter(name='xmin', value='0')
opObj11.addParameter(name='xmax', value='24')
#opObj11.addParameter(name='ymin', value='28')
#opObj11.addParameter(name='ymax', value='42')
# opObj11.addParameter(name='save', value=figpath)

op9 = proc_spectra.addOperation(name='RTIPlot')
op9.addParameter(name='id', value='20')
op9.addParameter(name='wintitle', value='RTI')
op9.addParameter(name='xmin', value='0')
op9.addParameter(name='xmax', value='24')
op9.addParameter(name='zmin', value=db_range[0])
op9.addParameter(name='zmax', value=db_range[1])
op9.addParameter(name='showprofile', value='1')
op9.addParameter(name='timerange', value=str(24))
# op9.addParameter(name='save', value=figpath)

op44 = proc_spectra.addOperation(name='selectHeights')
op44.addParameter(name='minIndex', value='0')
op44.addParameter(name='maxIndex', value='66')

op19 = proc_spectra.addOperation(name='SpectraPlot')
op19.addParameter(name='id', value='205')
op19.addParameter(name='wintitle', value='Spc DP Inst')
op19.addParameter(name='showprofile', value='1')
op19.addParameter(name='zmin', value=db_range[0])
op19.addParameter(name='zmax', value=db_range[1])

op08 = proc_spectra.addOperation(name='IntegrationFaradaySpectra')
op08.addParameter(name='n', value='214') #107

op6 = proc_spectra.addOperation(name='removeDCLagFlip')

op9 = proc_spectra.addOperation(name='SpectraPlot')
op9.addParameter(name='id', value='20')
op9.addParameter(name='wintitle', value='Spectra DP')
op9.addParameter(name='zmin', value=db_range[0])
op9.addParameter(name='zmax', value=db_range[1])
#op9.addParameter(name='xaxis', value='frequency')
op9.addParameter(name='showprofile', value='1')
op9.addParameter(name='show', value='1')
op9.addParameter(name='save', value=figpath)

op10 = proc_spectra.addOperation(name='SpectraDataToFaraday')

op11 = proc_spectra.addOperation(name='FlagBadHeightsSpectra')

op13 = proc_spectra.addOperation(name='DoublePulseACFs_PerLag')

op14 = proc_spectra.addOperation(name='FaradayAngleAndDPPower')

op15 = proc_spectra.addOperation(name='IGRFModel')

op16 = proc_spectra.addOperation(name='ElectronDensityFaraday')
op16.addParameter(name='NSHTS', value='60')
op16.addParameter(name='RATE', value='1.8978873e-6')

op17 = proc_spectra.addOperation(name='NormalizeDPPowerRoberto_V2')

op18 = proc_spectra.addOperation(name='DPTemperaturesEstimation')
op18.addParameter(name='IBITS', value='16')

op19 = proc_spectra.addOperation(name='ACFs')


op20 = proc_spectra.addOperation(name='ACFsPlot')
op20.addParameter(name='id', value='177')
op20.addParameter(name='wintitle', value='ACFs')
op20.addParameter(name='ymin', value='180')
op20.addParameter(name='ymax', value='600')
op20.addParameter(name='save', value=figpath_server)

op21 = proc_spectra.addOperation(name='EDensityPlot')
op21.addParameter(name='id', value='179')
op21.addParameter(name='wintitle', value='Electron Density')
op21.addParameter(name='ymin', value='180')
op21.addParameter(name='ymax', value='920')
op21.addParameter(name='xmin', value='1e3')
op21.addParameter(name='xmax', value='1e7')
op21.addParameter(name='save', value=figpath_server)

op22 = proc_spectra.addOperation(name='TempsDPPlot')
op22.addParameter(name='id', value='175')
op22.addParameter(name='wintitle', value='Temperatures')
op22.addParameter(name='ymin', value='180')
op22.addParameter(name='ymax', value='650')
op22.addParameter(name='save', value=figpath_server)

op013 = proc_spectra.addOperation(name='DataSaveCleaner')

op16 = proc_spectra.addOperation(name='DenRTIPlot')
op16.addParameter(name='id', value='174')
op16.addParameter(name='wintitle', value='Electron Density RTI')
op16.addParameter(name='xmin', value='0')
op16.addParameter(name='xmax', value='24')
op16.addParameter(name='ymin', value='180')
op16.addParameter(name='ymax', value='920')
op16.addParameter(name='zmin', value='1e4')
op16.addParameter(name='zmax', value='1e7')
op16.addParameter(name='xrange', value=str(24))
op16.addParameter(name='save', value=figpath_server)

op17 = proc_spectra.addOperation(name='ETempRTIPlot')
op17.addParameter(name='id', value='175')
op17.addParameter(name='wintitle', value='Electron Temperature RTI')
op17.addParameter(name='xmin', value='0')
op17.addParameter(name='xmax', value='24')
op17.addParameter(name='ymin', value='180')
op17.addParameter(name='ymax', value='600')
op17.addParameter(name='zmin', value='100')
op17.addParameter(name='zmax', value='4000')
op17.addParameter(name='xrange', value=str(24))
op17.addParameter(name='save', value=figpath_server)

op18 = proc_spectra.addOperation(name='ITempRTIPlot')
op18.addParameter(name='id', value='176')
op18.addParameter(name='wintitle', value='Ion Temperature RTI')
op18.addParameter(name='xmin', value='0')
op18.addParameter(name='xmax', value='24')
op18.addParameter(name='ymin', value='180')
op18.addParameter(name='ymax', value='600')
op18.addParameter(name='zmin', value='100')
op18.addParameter(name='zmax', value='4000')
op18.addParameter(name='xrange', value=str(24))
op18.addParameter(name='save', value=figpath_server)

one = {'gdlatr': 'lat', 'gdlonr': 'lon', 'inttms': 'paramInterval'} #reader gdlatr-->lat only 1D

two = {
    'gdalt': 'heightList',   #<----- nmonics
    'NE': ('DensityFinal', 0),
    'DNE': ('EDensityFinal', 0),
    'TE': ('ElecTempFinal', 0),
    'DTE': ('EElecTempFinal', 0),
    'TI': ('IonTempFinal', 0),
    'DTI': ('EIonTempFinal', 0),
    } #writer
#f=open('/home/cportilla/proyecto/Faraday/moder_test.txt','r')
#file_contents=f.read()
ind = ['gdalt']
meta = {
    'kinst': 10, #instrument code
    'kindat': 1800, #type of data
    'catalog': {
        'principleInvestigator': 'Danny Scipión',
        'expPurpose': 'Electron Density',
        #'sciRemarks': file_contents
        },
    'header': {
        'analyst': 'D. Hysell'
    }
}
#f.close()

op_writer = proc_spectra.addOperation(name='MADWriter')
op_writer.addParameter(name='path', value=figpath)
op_writer.addParameter(name='format', value='hdf5')
op_writer.addParameter(name='oneDDict', value=json.dumps(one))
op_writer.addParameter(name='twoDDict', value=json.dumps(two))
op_writer.addParameter(name='ind2DList', value=json.dumps(ind))
op_writer.addParameter(name='metadata', value=json.dumps(meta))


controller.start()
