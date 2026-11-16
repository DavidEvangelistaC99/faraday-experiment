#!/usr/bin/env python
'''
Created on Nov 04, 2026
@author: David Evangelista
'''

import os, sys
from schainpy.controller import Project

desc = "Hybrid (LP+DP): High Altitude Experiment"
controllerObj = Project()
controllerObj.setup(id = '191', name='test01', description=desc)

#data path 
dpath = '/media/cportilla/HDD/DAVID-OP/data-sevag/7'

#results path
figpath = '/media/cportilla/HDD/DAVID-OP/DATA-3/HAE/28_Oct_2021'

t = ['10.','19.']
#change db range
#db_range = ['0','50']
#to RTI
#db_range = ['-5','25']

db_range = ['-2.5','15']

#date and time (LT)
startDate = '2021/10/28'
endDate = '2021/10/28'
startTime = '10:00:00'
endTime = '18:59:59'

readUnitConfObj = controllerObj.addReadUnit(datatype='VoltageReader',
                                            path=dpath,
                                            startDate=startDate,
                                            endDate=endDate,
                                            startTime=startTime,
                                            endTime=endTime,
                                            online=0,
                                            walk=1,
                                            expLabel='',
                                            getByBlock=0,
                                            delay=30)

procUnitConfObj1 = controllerObj.addProcUnit(datatype='VoltageProc', inputId=readUnitConfObj.getId())

opObj11 = procUnitConfObj1.addOperation(name='selectChannels')
opObj11.addParameter(name='channelList', value='0,1')

#we receive 128 profiles plus 22 noise profiles
opObj11 = procUnitConfObj1.addOperation(name='ProfileSelector')
opObj11.addParameter(name='profileRangeList', value='(0,127)')

opObj11 = procUnitConfObj1.addOperation(name='selectHeights', optype='other')
opObj11.addParameter(name='minIndex', value='0')
opObj11.addParameter(name='maxIndex', value='800')
#0, 800 indexes represent 0 to 3000 km

#reorder the data for spectrogram analysis
'''opObj11 =  procUnitConfObj1.addOperation(name='SSheightProfiles', optype='other')
opObj11.addParameter(name='step', value='5', format='int')
#nsamples is equal to pulse width samples
opObj11.addParameter(name='nsamples', value='100', format='int')'''

procUnitConfObj1SPC = controllerObj.addProcUnit(datatype='Spectra', inputId=procUnitConfObj1.getId())
procUnitConfObj1SPC.addParameter(name='nFFTPoints', value='128', format='int')
procUnitConfObj1SPC.addParameter(name='nProfiles', value='128', format='int')

opObj11 = procUnitConfObj1SPC.addOperation(name='IncohInt', optype='other')
opObj11.addParameter(name='n', value='10', format='int')
#nIncohInt 5 is equal to 30 sec

op9 = procUnitConfObj1SPC.addOperation(name='SpectraPlot')
op9.addParameter(name='id', value='20')
op9.addParameter(name='wintitle', value='Spectra DP')
op9.addParameter(name='zmin', value=db_range[0])
op9.addParameter(name='zmax', value=db_range[1])
# op9.addParameter(name='xaxis', value='frequency')
op9.addParameter(name='showprofile', value='1')
op9.addParameter(name='show', value='1')
op9.addParameter(name='save', value=figpath)

'''op9 = procUnitConfObj1SPC.addOperation(name='RTIPlot')
op9.addParameter(name='id', value='20')
op9.addParameter(name='wintitle', value='RTI')
op9.addParameter(name='xmin', value=t[0])
op9.addParameter(name='xmax', value=t[1])
op9.addParameter(name='zmin', value=db_range[0])
op9.addParameter(name='zmax', value=db_range[1])
op9.addParameter(name='showprofile', value='1')
op9.addParameter(name='timerange', value=str(24))
op9.addParameter(name='save_period', value=3600)
op9.addParameter(name='save', value=figpath)'''

'''opObj11 = procUnitConfObj1SPC.addOperation(name='CoherencePlot', optype='other')
opObj11.addParameter(name='id', value='101', format='int')
opObj11.addParameter(name='wintitle', value='High altitude', format='str')
opObj11.addParameter(name='coherence_cmap', value='jet', format='str')
opObj11.addParameter(name='xmin', value=t[0], format='float')
opObj11.addParameter(name='xmax', value=t[1], format='float')
#opObj11.addParameter(name='showprofile', value='1')
opObj11.addParameter(name='save', value=figpath, format='str')'''

'''
opObj11 = procUnitConfObj1SPC.addOperation(name='SpectrogramPlot', optype='other')
opObj11.addParameter(name='id', value='1', format='int')
opObj11.addParameter(name='wintitle', value='High altitude', format='str')
opObj11.addParameter(name='zmin', value=db_range[0], format='int')
opObj11.addParameter(name='zmax', value=db_range[1], format='int')
opObj11.addParameter(name='xmin', value=t[0], format='float')
opObj11.addParameter(name='xmax', value=t[1], format='float')
opObj11.addParameter(name='save', value=figpath)'''

controllerObj.start()

