'''
#!/usr/bin/python3.7.4
'''



from schainpy.controller import Project

import json

controller = Project()
controller.setup(id = '001',
                 name='Hybrid',
                 description='DP+LP')

read_unit = controller.addReadUnit(datatype='VoltageReader',
                                   #path='/home/roberto/Documents/data',
                                   #path='/data/DP_19_12_Campaing',
                                   #######path='/home/roberto/pum',
                                   path='/home/roberto/puma/2020_09/lp_faraday',
                                   ###path='/home/roberto/puma/2020_02/hybrid2',
                                   #path='/data',
                                   #path='/home/roberto/dataHybrid/hybrid',
                                   ##path='/home/roberto/dataHybrid',
                                   #######startDate='2019/10/14',
                                   startDate='2020/09/03',
                                   #######endDate='2019/10/20',
                                   endDate='2020/09/03',
                                   #startTime='05:58:26',
                                   startTime='15:00:00',
                                   endTime='23:59:59',
                                   online=0,
                                   getByBlock='True',
                                   walk=0,
                                   #nTries=2,
                                   #nFiles=1,
                                   delay=30)


#############################################
proc_unit = controller.addProcUnit(datatype='VoltageProc',inputId=read_unit.getId())
op1 = proc_unit.addOperation(name='selectChannels')
op1.addParameter(name='channelList', value='0,1,2,3')

op121 = proc_unit.addOperation(name='ProfileSelector')
op121.addParameter(name='profileRangeList', value='0,149')

#op02 = proc_unit.addOperation(name='DuplicateData', optype='self')

#op002 = proc_unit.addOperation(name='filterByHeights2', optype='self')
#op002.addParameter(name='window', value='2', format='int')

op2 = proc_unit.addOperation(name='filterByHeights')
op2.addParameter(name='window', value='4')

#opp2 = proc_unit.addOperation(name='test111', optype='other')


op3 = proc_unit.addOperation(name='CrossProdHybrid')
op3.addParameter(name='NLAG', value='16')
op3.addParameter(name='NRANGE', value='200')
op3.addParameter(name='NCAL', value='0')
op3.addParameter(name='DPL', value='11')
op3.addParameter(name='NDN', value='0')
op3.addParameter(name='NDT', value='67')
op3.addParameter(name='NDP', value='67')
op3.addParameter(name='NSCAN', value='128')
op3.addParameter(name='lagind', value='(0,1,2,3,4,5,6,7,0,3,4,5,6,8,9,10)')
op3.addParameter(name='lagfirst', value='(1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1)')
op3.addParameter(name='NAVG', value='16')
op3.addParameter(name='nkill', value='6')

#op_LP_3 = proc_unit.addOperation(name='RemoveDebris', optype='other')




figpath='/home/roberto/Pictures/Test_Hybrid_2019_09_Reader+Fitter'
#figpath='/home/roberto/Pictures'
figpath='/home/roberto/Pictures/Test_Hybrid_2020_02_Reader_V2'
figpath='/home/roberto/Pictures/Hybrid_develop_manual'


op082 = proc_unit.addOperation(name='RTILPPlot')
op082.addParameter(name='id', value='177')
op082.addParameter(name='wintitle', value='RTI LP')
op082.addParameter(name='xmin', value='0')
op082.addParameter(name='xmax', value='24')
#op082.addParameter(name='ymin', value='180')
#op082.addParameter(name='ymax', value='600')
op082.addParameter(name='zmin', value='32')
op082.addParameter(name='zmax', value='80')
#op16.addParameter(name='showprofile', value='1')
#op16.addParameter(name='timerange', value=str(24*60*60))
op082.addParameter(name='xrange', value=str(24*1))
#op082.addParameter(name='save', value=figpath)



op080 = proc_unit.addOperation(name='CrossProductsPlot')
op080.addParameter(name='id', value='177')
op080.addParameter(name='wintitle', value='CrossProducts')
#op7.addParameter(name='xmin', value='0')
#op7.addParameter(name='xmax', value='24')
#op080.addParameter(name='ymin', value='180')
#op080.addParameter(name='ymax', value='600')
#op16.addParameter(name='zmin', value='32')
#op16.addParameter(name='zmax', value='62')
#op16.addParameter(name='showprofile', value='1')

#op16.addParameter(name='timerange', value=str(24*60*60))
#op080.addParameter(name='save', value=figpath)
#op080.addParameter(name='save_period', value=10)



op082 = proc_unit.addOperation(name='RTIDPPlot')
op082.addParameter(name='id', value='177')
op082.addParameter(name='wintitle', value='RTI')
op082.addParameter(name='xmin', value='0')
op082.addParameter(name='xmax', value='24')
#op082.addParameter(name='ymin', value='180')
#op082.addParameter(name='ymax', value='600')
op082.addParameter(name='zmin', value='32')
op082.addParameter(name='zmax', value='62')
#op16.addParameter(name='showprofile', value='1')

#op16.addParameter(name='timerange', value=str(24*60*60))
op082.addParameter(name='xrange', value=str(24*1), format='int')
#op082.addParameter(name='save', value=figpath, format='str')



op081 = proc_unit.addOperation(name='NoiseDPPlot')
op081.addParameter(name='id', value='178')
op081.addParameter(name='wintitle', value='Noise')
op081.addParameter(name='xmin', value='0')
op081.addParameter(name='xmax', value='24')
#op080.addParameter(name='ymin', value='180')
#op080.addParameter(name='ymax', value='600')
op081.addParameter(name='zmin', value='20')
op081.addParameter(name='zmax', value='65')
#op16.addParameter(name='showprofile', value='1')
op081.addParameter(name='xrange', value=str(24*1))
#op16.addParameter(name='timerange', value=str(24*60*60))
#op081.addParameter(name='save', value=figpath)
'''


'''
op080 = proc_unit.addOperation(name='CrossProductsLPPlot')
op080.addParameter(name='id', value='1777')
op080.addParameter(name='wintitle', value='CrossProducts')
#op080.addParameter(name='ymin', value='180')
#op080.addParameter(name='ymax', value='600')
#op080.addParameter(name='save', value=figpath)




op77_LP = proc_unit.addOperation(name='XmitWaveformPlot') #only works as it should when matplotlib.version>=3.2
op77_LP.addParameter(name='id', value='65')
op77_LP.addParameter(name='wintitle', value='Xmit Waveform')
#op77_LP.addParameter(name='save', value=figpath)




op_LP_3 = proc_unit.addOperation(name='RemoveDebris')


'''
o=proc_unit.addOperation(name='printAttribute')
o.addParameter(name='attributes',value='flagNoData')
'''
'''
dataList=['kax','kay','kbx','kby','kax2','kay2','kbx2',
        'kby2','kaxbx','kaxby','kaybx','kayby','kaxay',
        'kbxby','FirstAVGDate','LastAVGDate',
        'FirstBlockDate',
        'utctime','TimeBlockSeconds','output_LP','noise_final']
metadataList=['NAVG','DH','H0','DPL','NDP','NR','FirstHeight',
            'MAXNRANGENDT','NLAG','heightList','timeZone','NRANGE','NSCAN']


op221 = proc_unit.addOperation(name='HDFWriter')
op221.addParameter(name='path', value='/home/roberto/hdf5_outputs/test')
op221.addParameter(name='dataList', value=dataList)
op221.addParameter(name='metadataList', value=metadataList)
op221.addParameter(name='blocksPerFile', value=1000)
'''

'''
op21 = proc_unit.addOperation(name='DatWriter_New_V1')
op21.addParameter(name='experiment', value='HP')
op21.addParameter(name='path', value='/home/roberto/DatOut_Sch/Test_Hybrid/2020_09/test_reader')
op21.addParameter(name='format', value='out')
'''



op4 = proc_unit.addOperation(name='IntegrationHP')
op4.addParameter(name='nint', value='2')


op5 = proc_unit.addOperation(name='SumFlipsHP')

op6 = proc_unit.addOperation(name='FlagBadHeights')

#op7 = proc_unit.addOperation(name='GetNoise')
#op7 = proc_unit.addOperation(name='NoisePower')

op8 = proc_unit.addOperation(name='DoublePulseACFs')

op9 = proc_unit.addOperation(name='FaradayAngleAndDPPower')

op10 = proc_unit.addOperation(name='IGRFModel')

op11 = proc_unit.addOperation(name='ElectronDensityFaraday')
op11.addParameter(name='NSHTS', value='41')
op11.addParameter(name='RATE', value='1.8978873e-6')

op12 = proc_unit.addOperation(name='NormalizeDPPower')
#op12.addParameter(name='cut0', value='170')
#op12.addParameter(name='cut1', value='650')

op022 = proc_unit.addOperation(name='ACFs')

op012 = proc_unit.addOperation(name='DPTemperaturesEstimation')
op012.addParameter(name='IBITS', value='16')

op0122 = proc_unit.addOperation(name='LongPulseAnalysis')
op0122.addParameter(name='NACF', value='90')



#figpath='/home/roberto/Pictures/Test_Hybrid_2020_02_V2'
#figpath='/home/roberto/Pictures/Test_Hybrid_2019_09_Reader+Fitter'
#'''
op026 = proc_unit.addOperation(name='ACFsPlot')
op026.addParameter(name='id', value='177')
op026.addParameter(name='wintitle', value='ACFs DP')
#op7.addParameter(name='xmin', value='0')
#op7.addParameter(name='xmax', value='24')
#op016.addParameter(name='ymin', value='180')
#op016.addParameter(name='ymax', value='600')
#op16.addParameter(name='timerange', value=str(24*60*60))
#op026.addParameter(name='save', value=figpath)


op026 = proc_unit.addOperation(name='ACFsLPPlot')
op026.addParameter(name='id', value='177')
op026.addParameter(name='wintitle', value='ACFs LP')
#op7.addParameter(name='xmin', value='0')
#op7.addParameter(name='xmax', value='24')
op026.addParameter(name='ymin', value='550')
#op026.addParameter(name='ymax', value='600')
#op16.addParameter(name='zmin', value='32')
#op16.addParameter(name='zmax', value='62')
#op16.addParameter(name='showprofile', value='1')

#op16.addParameter(name='timerange', value=str(24*60*60))
#op026.addParameter(name='save', value=figpath)


op026 = proc_unit.addOperation(name='EDensityHPPlot')
op026.addParameter(name='id', value='179')
op026.addParameter(name='wintitle', value='Electron Density')
#op7.addParameter(name='xmin', value='0')
#op7.addParameter(name='xmax', value='24')
op026.addParameter(name='ymin', value='150')
op026.addParameter(name='ymax', value='1350')
#op16.addParameter(name='zmin', value='1000')
#op16.addParameter(name='zmax', value='10000000')
#op16.addParameter(name='showprofile', value='1')

#op16.addParameter(name='timerange', value=str(24*60*60))
#op026.addParameter(name='save', value=figpath)


op0166 = proc_unit.addOperation(name='TempsHPPlot')
op0166.addParameter(name='id', value='1755')
op0166.addParameter(name='wintitle', value='Temperatures')
#op7.addParameter(name='xmin', value='0')
#op7.addParameter(name='xmax', value='24')
op0166.addParameter(name='ymin', value='160')
#op016.addParameter(name='ymax', value='600')
#op16.addParameter(name='zmin', value='32')
#op16.addParameter(name='zmax', value='62',)
#op16.addParameter(name='showprofile', value='1')

#op16.addParameter(name='timerange', value=str(24*60*60))
#op0166.addParameter(name='save', value=figpath)




op0177 = proc_unit.addOperation(name='FracsHPPlot')
op0177.addParameter(name='id', value='1700')
op0177.addParameter(name='wintitle', value='Composition')
#op7.addParameter(name='xmin', value='0')
#op7.addParameter(name='xmax', value='24')
#op016.addParameter(name='ymin', value='180')
#op016.addParameter(name='ymax', value='600')
#op16.addParameter(name='zmin', value='32')
#op16.addParameter(name='zmax', value='62')
#op16.addParameter(name='showprofile', value='1')

#op16.addParameter(name='timerange', value=str(24*60*60))
#op0177.addParameter(name='save', value=figpath)

#'''


op018 = proc_unit.addOperation(name='DataSaveCleanerHP')


op181 = proc_unit.addOperation(name='ETempRTIPlot')
op181.addParameter(name='id', value='176')
op181.addParameter(name='wintitle', value='ElectronTemp')
op181.addParameter(name='xmin', value='0')
op181.addParameter(name='xmax', value='24')
op181.addParameter(name='ymin', value='180')
op181.addParameter(name='ymax', value='1350')
op181.addParameter(name='zmin', value='500')
op181.addParameter(name='zmax', value='6500')
#op16.addParameter(name='showprofile', value='1')
op181.addParameter(name='xrange', value=str(24*1))
#op16.addParameter(name='timerange', value=str(24*60*60))
#op181.addParameter(name='save', value=figpath)


op182 = proc_unit.addOperation(name='ITempRTIPlot')
op182.addParameter(name='id', value='176')
op182.addParameter(name='wintitle', value='IonTemp')
op182.addParameter(name='xmin', value='0')
op182.addParameter(name='xmax', value='24')
op182.addParameter(name='ymin', value='180')
op182.addParameter(name='ymax', value='1350')
op182.addParameter(name='zmin', value='500')
op182.addParameter(name='zmax', value='6500')
#op16.addParameter(name='showprofile', value='1')
op182.addParameter(name='xrange', value=str(24*1))
#op16.addParameter(name='timerange', value=str(24*60*60))
#op182.addParameter(name='save', value=figpath)



op183 = proc_unit.addOperation(name='HFracRTIPlot')
op183.addParameter(name='id', value='176')
op183.addParameter(name='wintitle', value='H+ Frac')
op183.addParameter(name='xmin', value='0')
op183.addParameter(name='xmax', value='24')
op183.addParameter(name='ymin', value='180')
op183.addParameter(name='ymax', value='1350')
op183.addParameter(name='zmin', value='0')
op183.addParameter(name='zmax', value='1')
#op16.addParameter(name='showprofile', value='1')
op183.addParameter(name='xrange', value=str(24*1))
#op16.addParameter(name='timerange', value=str(24*60*60))
#op183.addParameter(name='save', value=figpath)


op184 = proc_unit.addOperation(name='HeFracRTIPlot')
op184.addParameter(name='id', value='176')
op184.addParameter(name='wintitle', value='He+ Frac')
op184.addParameter(name='xmin', value='0')
op184.addParameter(name='xmax', value='24')
op184.addParameter(name='ymin', value='180')
op184.addParameter(name='ymax', value='1350')
op184.addParameter(name='zmin', value='0')
op184.addParameter(name='zmax', value='0.5')
#op16.addParameter(name='showprofile', value='1')
op184.addParameter(name='xrange', value=str(24*1))
#op16.addParameter(name='timerange', value=str(24*60*60))
#op184.addParameter(name='save', value=figpath)



op16 = proc_unit.addOperation(name='DenRTIPlot')
op16.addParameter(name='id', value='174')
op16.addParameter(name='wintitle', value='Electron Density RTI')
op16.addParameter(name='xmin', value='0')
op16.addParameter(name='xmax', value='24')
op16.addParameter(name='ymin', value='180')
op16.addParameter(name='ymax', value='1350')
op16.addParameter(name='zmin', value='1000')
op16.addParameter(name='zmax', value='10000000')
#op16.addParameter(name='showprofile', value='1')

op16.addParameter(name='xrange', value=str(24*1))
#op16.addParameter(name='save', value=figpath)





'''


one = {'gdlatr': 'lat', 'gdlonr': 'lon'} #reader gdlatr-->lat only 1D





two_wr = {
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
    'DPHE+': ('EPheFinal', 0),
    'SRACF0': ('acfs_DP', 0),
    'DSRACF0': ('acfs_error_DP', 0),
    'RACF1': ('acfs_DP', 1),
    'DRACF1': ('acfs_error_DP', 1),
    'RACF2': ('acfs_DP', 2),
    'DRACF2': ('acfs_error_DP', 2),
    'RACF3': ('acfs_DP', 3),
    'DRACF3': ('acfs_error_DP', 3),
    'RACF4': ('acfs_DP', 4),
    'DRACF4': ('acfs_error_DP', 4),
    'RACF5': ('acfs_DP', 5),
    'DRACF5': ('acfs_error_DP', 5),
    'RACF6': ('acfs_DP', 6),
    'DRACF6': ('acfs_error_DP', 6),
    'RACF7': ('acfs_DP', 7),
    'DRACF7': ('acfs_error_DP', 7),
    'RACF8': ('acfs_DP', 8),
    'DRACF8': ('acfs_error_DP', 8),
    'RACF9': ('acfs_DP', 9),
    'DRACF9': ('acfs_error_DP', 9),
    'RACF10': ('acfs_DP', 10),
    'DRACF10': ('acfs_error_DP', 10),
    'RACF11': ('acfs_LP', 0),
    'DRACF11': ('acfs_error_LP', 0),
    'RACF12': ('acfs_LP', 1),
    'DRACF12': ('acfs_error_LP', 1),
    'RACF13': ('acfs_LP', 2),
    'DRACF13': ('acfs_error_LP', 2),
    'RACF14': ('acfs_LP', 3),
    'DRACF14': ('acfs_error_LP', 3),
    'RACF15': ('acfs_LP', 4),
    'DRACF15': ('acfs_error_LP', 4),
    'RACF16': ('acfs_LP', 5),
    'DRACF16': ('acfs_error_LP', 5),
    'RACF17': ('acfs_LP', 6),
    'DRACF17': ('acfs_error_LP', 6),
    'RACF18': ('acfs_LP', 7),
    'DRACF18': ('acfs_error_LP', 7),
    'RACF19': ('acfs_LP', 8),
    'DRACF19': ('acfs_error_LP', 8),
    'RACF20': ('acfs_LP', 9),
    'DRACF20': ('acfs_error_LP', 9),
    'RACF21': ('acfs_LP', 10),
    'DRACF21': ('acfs_error_LP', 10),
    'RACF22': ('acfs_LP', 11),
    'DRACF22': ('acfs_error_LP', 11),
    'RACF23': ('acfs_LP', 12),
    'DRACF23': ('acfs_error_LP', 12),
    'RACF24': ('acfs_LP', 13),
    'DRACF24': ('acfs_error_LP', 13),
    'RACF25': ('acfs_LP', 14),
    'DRACF25': ('acfs_error_LP', 14),
    'RACF26': ('acfs_LP', 15),
    'DRACF26': ('acfs_error_LP', 15)


    # 'DVIPN': ('data_output', 3),
    #'snl': 'data_SNR'
    } #writer

f=open('/home/roberto/moder_test.txt','r')
file_contents=f.read()


ind = ['gdalt']

meta = {
    'kinst': 10, #instrumnet code
    'kindat': 1802, #type of data  ###to be confirmed
    'catalog': {
        'principleInvestigator': 'Marco Milla',
        'expPurpose': 'Electron Density',
        'sciRemarks': file_contents
        },
    'header': {
        'analyst': 'David Hysell'
    }
}

f.close()




op_writer = proc_unit.addOperation(name='MADWriter')
#op_writer.addParameter(name='path', value='/home/roberto/hdf5_outputs')
op_writer.addParameter(name='path', value='/home/roberto/hdf5_outputs/Hybrid_2020_09_reader')
op_writer.addParameter(name='format', value='hdf5')
op_writer.addParameter(name='oneDDict', value=json.dumps(one))
op_writer.addParameter(name='twoDDict', value=json.dumps(two_wr))
op_writer.addParameter(name='ind2DList', value=json.dumps(ind))
op_writer.addParameter(name='metadata', value=json.dumps(meta))

#op_writer.addParameter(name='blocks', value=60)
'''




controller.start()
