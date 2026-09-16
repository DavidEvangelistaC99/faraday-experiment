'''
#!/usr/bin/python3.7.4
'''

from schainpy.controller import Project
import json

controller = Project()
controller.setup(id = '001',
                 name='Faraday',
                 description='DP')


figpath = '/home/david/Documents/DATA-3/Faraday/18_21_Aug_26/clean-test/20_Aug_26-4/'
# figpath = '/home/david/Documents/DATA-3/Faraday/18_21_Aug_26/20_Aug_26/'
figpath_server=figpath

# data_path = '/mnt/isr_all/Faraday'
data_path = '/mnt/share/2026_08/Faraday/main_radar/rawdata'

startDate = '2026/08/20'
endDate = '2026/08/20'
startTime='00:00:00'
endTime='23:59:59'

startTime_ = startTime[:2]
endTime_ = endTime[:2]

Lag = '0'
db_range = ['30','40']

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

'''opObj11 = proc_spectra.addOperation(name='NoisePlot')
opObj11.addParameter(name='id', value='3')
opObj11.addParameter(name='wintitle', value='Noise')
opObj11.addParameter(name='xmin', value='0')
opObj11.addParameter(name='xmax', value='24')
#opObj11.addParameter(name='ymin', value='28')
#opObj11.addParameter(name='ymax', value='42')
# opObj11.addParameter(name='save', value=figpath)'''



op9 = proc_spectra.addOperation(name='RTIPlot')
op9.addParameter(name='id', value='20')
op9.addParameter(name='wintitle', value='RTI')
op9.addParameter(name='xmin', value='0')
op9.addParameter(name='xmax', value='24')
op9.addParameter(name='zmin', value=db_range[0])
op9.addParameter(name='zmax', value=db_range[1])
op9.addParameter(name='showprofile', value='1')
op9.addParameter(name='timerange', value=str(24))
op9.addParameter(name='save_period', value=3600)
#op9.addParameter(name='show', value='1')
op9.addParameter(name='save', value=figpath)

controller.start()
