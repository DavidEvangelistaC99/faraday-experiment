"""testMadrigalWebServices.py runs a test of the Madrigal Web Services interface
   for a given Madrigal server.

   usage:

   python testMadrigalWebServices.py [madrigal main url]

   If no madrigal main url given, tries to get it from madrigal.metadata if installed.
   Otherwise, prints usage with madrigal main url required.

   Prints "success" if no errors raised - no longer compares to std output file

   Assumes the basic Madrigal test files have been installed.

"""

# $Id: testMadrigalWebServices.py 6378 2017-12-11 20:42:54Z brideout $

success = True
try:

    import sys
    import traceback
    import tempfile
    import os, os.path
    
    import madrigalWeb.madrigalWeb

    if len(sys.argv) < 2:
        # see if we can get the url from madrigal.metadata
        try:
            import madrigal.metadata
            madDB = madrigal.metadata.MadrigalDB()
            madrigalUrl = madDB.getTopLevelUrl()
        except:
            print('usage: python testMadrigalWebServices.py <madrigal main url>')
            sys.exit(-1)
    else:
        madrigalUrl = sys.argv[1]
        
    tempDir = tempfile.gettempdir()

    # constants
    user_fullname = 'Bill Rideout - automated test'
    user_email = 'brideout@haystack.mit.edu'
    user_affiliation = 'MIT Haystack'

    testData = madrigalWeb.madrigalWeb.MadrigalData(madrigalUrl)

    instList = testData.getAllInstruments()

    # print out Millstone
    for inst in instList:
        if inst.code == 30:
            print(str(inst) + '\n')
            
    print('')

    print('Test of getExperiments:')
    expList = testData.getExperiments(30, 1998,1,19,0,0,0,1998,1,22,0,0,0)
    expList.sort()

    for exp in expList:
        # should be only one
        print(str(exp) + '\n')
        
    print('')

    print('Test of getExperimentFiles:')
    fileList = testData.getExperimentFiles(expList[0].id)

    thisFilename = None

    for file in fileList:
        if file.category == 1:
            print(str(file) + '\n')
            thisFilename = file.name
            break

    if thisFilename == None:
        thisFilename = fileList[0].name
        
    print('')

    print('Test of downloadFile - simple format:')
    result = testData.downloadFile(fileList[0].name, os.path.join(tempDir, "test.txt"), 
                                   user_fullname, user_email, user_affiliation, "simple")
    
    print('')

    print('Test of downloadFile - hdf5 format:')
    try:
        result = testData.downloadFile(fileList[0].name, os.path.join(tempDir, "test.hdf5"), 
                                       user_fullname, user_email, user_affiliation, "hdf5")
    except:
        traceback.print_exc()
        print('Test of downloadFile - hdf5 format failed')
        success = False
        
    print('')

    print('Test of isprint:')
    print(testData.isprint(thisFilename,
                                   'gdalt, ti',
                                   'filter=gdalt,500,600 filter=ti,1900,2.0E+3',
                                   user_fullname, user_email, user_affiliation))
    print('')

    print('Test 1 of madCalculator:')
    result = testData.madCalculator(1999,2,15,12,30,0,45,55,5,-170,-150,10,200,2.0E+2,0,
                                    'sdwht, kp', ['kinst'], [30])

    for line in result:
        print(line)

    print('')
    
    print('Test 2 of madTimeCalculator:')
    result = testData.madTimeCalculator(1999,2,15,12,30,0,1999,2,15,13,30,0,0.025E+1,'kp, ap3')
    
    for line in result:
        print(line)

    print('')
    
    print('Test of madCalculator2:')
    try:
        result = testData.madCalculator2(1999,2,15,12,30,0,[45,55],[-170,-150],[200,300],'bmag, pdcon',
                                         ['kp'],[1.0],['ti','te','ne'],
                                         [[1000,1000],[1100,1200],[1e+11,1.2e+11]])

        for line in result:
            print(line)
    except:
        traceback.print_exc()
        print('Test of madCalculator2 failed')
        success = False

    print('')
    
    # two tests of madCalculator3 - with and without 1 and 2D values
    print('Test 1 of madCalculator3:')
    try:
        result = testData.madCalculator3(yearList=[2001,2001], monthList=[3,3], dayList=[19,20],
                                         hourList=[12,12], minList=[30,40], secList=[20,0],
                                         latList=[[45,46,47,48.5],[46,47,48.2,49,50]],
                                         lonList=[[-70,-71,-72,-73],[-70,-71,-72,-73,-74]],
                                         altList=[[145,200,250,300.5],[2.0E+2,250,300,350,400]],
                                         parms='bmag,pdcon, ne_model',
                                         oneDParmList=[],
                                         oneDParmValues=[],
                                         twoDParmList=[],
                                         twoDParmValues=[])
    
        for line in result:
            print(line)
    except:
        traceback.print_exc()
        print('Test of madCalculator3 failed')
        success = False

    print('')
    
    print('Test 2 of madCalculator3:')
    try:
        result = testData.madCalculator3(yearList=[2001,2001], monthList=[3,3], dayList=[19,20],
                                         hourList=[12,12], minList=[30,40], secList=[20,0],
                                         latList=[[45,46,47,48.5],[46,47,48.2,49,50]],
                                         lonList=[[-70,-71,-72,-73],[-70,-71,-72,-73,-74]],
                                         altList=[[145,200,250,300.5],[200,250,300,350,400]],
                                         parms='bmag,pdcon, ne_model',
                                         oneDParmList=['kinst','elm'],
                                         oneDParmValues=[[31.0,31.0],[45.0,50.0]],
                                         twoDParmList=['ti','te','ne'],
                                         twoDParmValues=[[[1000,1000,1000,1000],[1000,1000,1000,1000,1000]],
                                                         [[1100,1200,1300,1400],[1500,1000,1100,1200,1300]],
                                                         [[1.0e10,1.0e+10,1.0e10,1.0e10],[1.0e10,1.0e10,1.0e10,1.0e10,1.0e10]]])
        for line in result:
            print(line)
    except:
        traceback.print_exc()
        print('Test 2 of madCalculator3 failed')
        success = False
        
        
    print('Test of geodeticToRadar')
    slatgd = 42.0
    slon = -70.0
    saltgd = 0.1
    try:
        result = testData.geodeticToRadar(slatgd, slon, saltgd, [50, 51,52], [-80.0, -70.0, -60.0], [200.0, 3.0E+2, 400.0])
        print(result)
    except:
        traceback.print_exc()
        print('Test of geodeticToRadar failed')
        success = False
        
        
    print('Test of radarToGeodetic')
    slatgd = 42.0
    slon = -70.0
    saltgd = 0.1
    # get azs, els, ranges from previous result
    azs = []
    els = []
    ranges = []
    for az, el, thisRange in result:
        azs.append(az)
        els.append(el)
        ranges.append(thisRange)
    try:
        result = testData.radarToGeodetic(slatgd, slon, saltgd, azs, els, ranges)
        print(result)
    except:
        traceback.print_exc()
        print('Test of radarToGeodetic failed')
        success = False
        
        
    print('Test of listFileTimes')
    expDir = 'experiments/1998/mlh/20jan98'
    try:
        result = testData.listFileTimes(expDir)
        print(result)
        print('')
    except:
        traceback.print_exc()
        print('Test of listFileTimes failed')
        success = False
        
    print('Test of getVersion')
    try:
        result = testData.getVersion()
        print(result)
        version = result
        print('')
    except:
        traceback.print_exc()
        print('Test of getVersion failed')
        success = False
        
    if testData.compareVersions(version, '2.7'):
        print('Test of downloadWebFile for this Madrigal site 3.0 or greater')
        for item in result:
            if item[0].find('plot') != -1:
                print('Downloading %s' % (item[0]))
                try:
                    testData.downloadWebFile(item[0], tempDir)
                except:
                    traceback.print_exc()
                    print('Test of downloadWebFile failed')
                    success = False
                break
                
        
    print('Test of traceMagneticField - igrf')
    year = 1998
    month = 1
    day = 20
    hour = 18
    minute = 30
    second = 0
    inputType = outputType = 0
    alts = [200, 300, 400]
    lats = [42,42,42]
    lons = [-70,-70,-70]
    model = 1 # igrf
    qualifier = 0 # comjugate
    stopAlt = 1000
    try:
        result = testData.traceMagneticField(year, month, day, hour, minute, second, 
                                             inputType, outputType, alts, lats, lons, model, qualifier, stopAlt)
        print(result)
        print('')
    except:
        traceback.print_exc()
        print('Test of traceMagneticField failed')
        success = False
        
        
    print('Test of traceMagneticField - Tsyganenko')
    year = 1998
    month = 1
    day = 20
    hour = 18
    minute = 30
    second = 0
    inputType = outputType = 0
    alts = [200, 300, 400]
    lats = [42,42,42]
    lons = [-70,-70,-70]
    model = 0 # Tsyganenko
    qualifier = 3 # Apex
    stopAlt = 1000
    try:
        result = testData.traceMagneticField(year, month, day, hour, minute, second, 
                                             inputType, outputType, alts, lats, lons, model, qualifier, stopAlt)
        print(result)
        print('')
    except:
        traceback.print_exc()
        print('Test of traceMagneticField failed')
        success = False
        
        
    

    print('')

    if success:
        print('\nMadrigalWebServices regression test result: SUCCESS')
    else:
        print('\nMadrigalWebServices regression test result: FAILURE')

except:
    traceback.print_exc()
    print('\nMadrigalWebServices regression test result: FAILURE')


