'''
views for docs app
 
@author: Bill Rideout
@contact: brideout@haystack.mit.edu

$Id: views.py 7140 2020-07-29 02:33:55Z brideout $
'''
# standard python imports
import os, os.path, sys
import random
import datetime
import re
import xml.dom.minidom
import hashlib

# django imports
import django, django.http
from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from . import forms

# madrigal imports
import madrigal.metadata
import madrigal.ui.web
import madrigal.openmadrigal
import subversion_operations

# helper methods
def get_exp_info(control_file, rev):
    """get_exp_info returns a tuple of (exp_desc, cyc_desc, control_file_text) for the input control_file and revision.
    
    For now grabs info from Subversion on apollo
    """
    exp_desc = ''
    cyc_desc = ''
    control_file_text = ''
    
    subObj = subversion_operations.reader()
    xmlFilename  = control_file.replace('.dat', '.xml')
    output = os.path.join('/tmp', control_file)
    sub_path = 'millstone_geospace_radar/control/ISRExperiments'
    
    try:
        subObj.getFile(os.path.join(sub_path, control_file), '/tmp', rev)
    except:
        try:
            subObj.getFile(os.path.join(sub_path, control_file), '/tmp')
        except:
            return((exp_desc, cyc_desc, control_file_text))
        
    f = open(output)
    control_file_text = f.read()
    f.close()
    
    try:
        subObj.getFile(os.path.join(sub_path, xmlFilename), '/tmp', rev)
    except:
        try:
            subObj.getFile(os.path.join(sub_path, xmlFilename), '/tmp')
        except:
            return((exp_desc, cyc_desc, control_file_text))
        
    f = open(output)
    xml_file_text = f.read()
    f.close()
    
    
    
    # now build complete xml file
    cvsText = '<?xml version="1.0"?>\n<ExperimentDescriptionFile>\n'
    cvsText += xml_file_text

    # closing tag
    cvsText += '</ExperimentDescriptionFile>'


    expDescDoc = xml.dom.minidom.parseString(cvsText)


    # get all experiments
    expElemList = expDescDoc.getElementsByTagName("Experiment")
    
    for exp in expElemList:
        
        # get one exp name from exp list
        nameEl = exp.getElementsByTagName("ExperimentName")[0]

        # get text node, convert from unicode to ascii
        expName = ''
        for elem in nameEl.childNodes:
            if elem.nodeType == elem.TEXT_NODE:
               expName += elem.data.encode('utf-8')

        # get exp descriptions from exp list
        expDescElList = exp.getElementsByTagName("ExperimentDescription")

        for expDescEl in expDescElList:

            # get text node, convert from unicode to ascii
            expDesc = ''
            for elem in expDescEl.childNodes:
                if elem.nodeType == elem.TEXT_NODE:
                   expDesc += elem.data.encode('utf-8')
            if len(exp_desc) == 0:
                exp_desc = expDesc # make sure we get something even if no id match

            # get one CVS_rev
            cvsEl = expDescEl.getElementsByTagName("CVS_rev")[0]

            cvs = ''
            for elem in cvsEl.childNodes:
                if elem.nodeType == elem.TEXT_NODE:
                   cvs += elem.data.encode('utf-8')

            if cvs == rev:
                # right rev - overwrite
                if len(expDesc) == 0:
                    exp_desc = expDesc
                

            

        # get cycle descriptions from exp list
        cycDescElList = exp.getElementsByTagName("CycleDescription")

        for cycDescEl in cycDescElList:

            # get text node, convert from unicode to ascii
            cycDesc = ''
            for elem in cycDescEl.childNodes:
                if elem.nodeType == elem.TEXT_NODE:
                   cycDesc += elem.data.encode('utf-8')
                   
            if len(cyc_desc) == 0:
                cyc_desc = cycDesc # make sure we get something even if no id match

            # get one CVS_rev
            cvsEl = cycDescEl.getElementsByTagName("CVS_rev")[0]

            cvs = ''
            for elem in cvsEl.childNodes:
                if elem.nodeType == elem.TEXT_NODE:
                   cvs += elem.data.encode('utf-8')

            if cvs == rev:
                # right rev - overwrite
                if len(cycDesc) == 0:
                    cyc_desc = cycDesc

        
    return((exp_desc, cyc_desc, control_file_text))
    
    

@csrf_exempt
def get_log_admin(request):
    """get_log_admin returns the admin access log page.  
    
    Inputs:
        request 
    """
    madDB = madrigal.metadata.MadrigalDB()
    bg_color = madDB.getBackgroundColor()
    
    if request.method == 'POST':
        form =forms.AdminLogForm(request.POST)
        if form.is_valid():
            
            cleaned_data = form.cleaned_data
            kinstList = [int(item) for item in eval(cleaned_data['kinst'])]
            startDT = datetime.datetime(cleaned_data['startDate'].year, cleaned_data['startDate'].month, cleaned_data['startDate'].day)
            endDT = datetime.datetime(cleaned_data['endDate'].year, cleaned_data['endDate'].month, cleaned_data['endDate'].day)
            madWebObj = madrigal.ui.web.MadrigalWeb(madDB)
            stageDir = os.path.join(madDB.getMadroot(), 'experiments/stage')
            tmpFile = os.path.join(stageDir, 'admin_%i.log' % (random.randint(0, 1000000)))
            madWebObj.filterLog(tmpFile, kinstList, startDT, endDT)
            f = open(tmpFile, 'r')
            thisFile = django.core.files.File(f)
            response = django.http.HttpResponse(thisFile, content_type='application/x-octet-stream')
            response['Content-Disposition'] = 'attachment; filename="' + os.path.basename(tmpFile) + '"'
            response['Content-Length'] = os.path.getsize(tmpFile)
            response['Set-Cookie'] = 'fileDownload=true; path=/'
            return(response)

    else:
        form =forms.AdminLogForm()
    
    responseDict = {'bg_color': bg_color}
    responseDict['form'] = form
    return render(request, 'get_log_admin.html', responseDict)


def get_latest_metadata_version(request):
    """get_latest_metadata_version returns the latest version of the specified metadata file
    
    Inputs:
        request - contains key fullPath - full path to the Madrigal file in Subversion relative to trunk
            (example: 'madroot/metadata/siteTab.txt')
    """
    madDB = madrigal.metadata.MadrigalDB()
    fullPath = request.GET['fullPath']
    
    subObj = subversion_operations.reader('OpenMadrigal')
    output = os.path.join('/tmp', os.path.basename(fullPath))
    path = os.path.join('trunk', fullPath)
    subObj.getFile(path, '/tmp')
    
    f = open(output)
    fileStr = f.read()
    f.close()
    try:
        os.remove(output)
    except:
        pass
    
    return(django.http.HttpResponse(fileStr))


def get_all_metadata_versions(request):
    """get_all_metadata_versions returns the a page listing all versions numbers available for a given metadata
    file, one line for each version string
    
    Inputs:
        request - contains key fullPath - full path to the Madrigal file in Subversion relative to trunk
            (example: 'madroot/metadata/siteTab.txt')
            
    """
    madDB = madrigal.metadata.MadrigalDB()
    fullPath = request.GET['fullPath']
    
    subObj = subversion_operations.reader('OpenMadrigal')
    path = os.path.join('trunk', fullPath)
    revDict = subObj.getFileRevsDates(path)
            
    retStr = ''
    for k in sorted(revDict):
        retStr += '%i\n' % (revDict[k])
    
    return(django.http.HttpResponse(retStr))


def get_metadata_version(request):
    """get_metadata_version returns the specifed version of the specified metadata file
    
    Inputs:
        request - contains key fullPath - full path to the Madrigal file in Subversion relative to trunk
            (example: 'madroot/metadata/siteTab.txt'), and key version
        
    """
    madDB = madrigal.metadata.MadrigalDB()
    fullPath = request.GET['fullPath']
    version = request.GET['version']
    
    subObj = subversion_operations.reader('OpenMadrigal')
    output = os.path.join('/tmp', os.path.basename(fullPath))
    path = os.path.join('trunk', fullPath)
    subObj.getFile(path, '/tmp', version)
    
    f = open(output, 'rb')
    text = f.read()
    f.close()
    try:
        os.remove(output)
    except:
        pass
    
    return(django.http.HttpResponse(text))


def get_open_madrigal_shared_files(request):
    """get_open_madrigal_shared_files returns the specifed version of the shared OpenMadrigal file 
    
    Inputs:
        request - contains key filename - in form siteName_siteID/expTab.txt or siteName_siteID/fileTab.txt
        
    For now hardcoded path to files
    """
    path = '/usr/local/apache2/htdocs/madrigal/distributionFiles/metadata'
    path3 = '/usr/local/apache2/htdocs/madrigal/distributionFiles/metadata3' # used because siteTab.txt for mad3 is incompatible
    filename = request.GET['filename']
    if filename == 'siteTab.txt':
        fullPath = os.path.join(path3, filename)
    else:
        fullPath = os.path.join(path, filename)
    f = open(fullPath, 'r')
    text = f.read()
    f.close()
    return(django.http.HttpResponse(text))
    
  
def get_madrigal_videos(request):
    """get_madrigal_videos returns the Madrigal video tutorials page 
    
    Inputs:
        request 
        
    """
    madDB = madrigal.metadata.MadrigalDB()
    bg_color = madDB.getBackgroundColor()
    responseDict = {'bg_color': bg_color}
    return render(request, 'get_video_tutorials.html', responseDict)

def get_exp_notes(request):
    """get_exp_notes returns the description of the control file - links from catalog in MLH files
    
    Inputs:
        request 
        
    """
    responseDict = {}
    responseDict['control_file'] = request.GET['exp']
    responseDict['rev'] = request.GET['rev']
    
    exp_desc, cyc_desc, control_file_text = get_exp_info(responseDict['control_file'],
                                                         responseDict['rev'])
    
    responseDict['exp_desc'] = exp_desc
    responseDict['cyc_desc'] = cyc_desc
    responseDict['control_file_text'] = control_file_text
    
    
    return render(request, 'exp_notes.html', responseDict)


def compare_to_archive(request):
    """compare_to_archive implements the old compareToArchive.py cgi script
    
    Now checks newest to oldest, up to 50 max
    """
    filePath = request.GET['filePath']
    fileTextMd5 = request.GET['fileTextMd5']
    strip = False
    if 'strip' in request.GET:
        strip = True
    madDB = madrigal.metadata.MadrigalDB()
    openMadObj = madrigal.openmadrigal.OpenMadrigal(madDB)
    revList = openMadObj.getAllRevisionNumbers(filePath)
    revList.reverse() # check newest first

    if len(revList) == 0:
        return(django.http.HttpResponse('None None'))
    
    # loop through all revisions, looking for a match
    for rev in revList[:50]:
        thisText = openMadObj.getCvsVersion(filePath, rev)
        if thisText is None:
            continue
        if strip:
            data = thisText.strip().encode()
            thisMd5 = hashlib.md5(data)
        else:
            thisMd5 = hashlib.md5(thisText.encode())
        if thisMd5.hexdigest() == fileTextMd5:
            return(django.http.HttpResponse('%s %s' % (revList[0], rev)))

    # no matches found
    return(django.http.HttpResponse('%s None' % (revList[0])))
 

def open_madrigal(request):
    """openmadrigal implements the openmadigal page
    """
    madDB = madrigal.metadata.MadrigalDB()
    madSiteObj = madrigal.metadata.MadrigalSite(madDB)
    responseDict = {}
    site_list = []
    for siteId, siteName in madSiteObj.getSiteList():
        if madSiteObj.getSiteServer(siteId).find('http') == -1:
            url = 'http://' + os.path.join(madSiteObj.getSiteServer(siteId),
                                           madSiteObj.getSiteDocRoot(siteId))
        else:
            url = os.path.join(madSiteObj.getSiteServer(siteId),
                               madSiteObj.getSiteDocRoot(siteId))
        site_list.append((url, siteName))
        
    responseDict['site_list'] = site_list
    return(render(request, 'openmadrigal.html', responseDict))

def madrigal_admin(request):
    """madrigal_admin implements the openmadigal admin page
    """
    responseDict = {}
    return(render(request, 'admin.html', responseDict))

def madrigal_download(request):
    """madrigal_download implements the openmadigal download page
    """
    responseDict = {}
    return(render(request, 'madDownload.html', responseDict))
    
        
       
    
    
    
