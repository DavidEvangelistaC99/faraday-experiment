'''
urls for docs app
 
@author: Bill Rideout
@contact: brideout@haystack.mit.edu

$Id: urls.py 7093 2020-01-23 19:46:59Z brideout $
'''


# django imports
from django.conf.urls import url
from . import views

urlpatterns = [ url(r'^getLogAdmin.py', 
                views.get_log_admin, name='getLogAdmin'),
               url(r'^getLatestMetadataVersion$', 
                views.get_latest_metadata_version, name='getLatestMetadataVersion'),
               url(r'^getAllMetadataVersions$', 
                views.get_all_metadata_versions, name='getAllMetadataVersions'),
               url(r'^getMetadataVersion$', 
                views.get_metadata_version, name='getMetadataVersion'),
               url(r'^getOpenMadrigalSharedFiles$', 
                views.get_open_madrigal_shared_files, name='getOpenMadrigalSharedFiles'),
               url(r'^getMadrigalVideos/$', 
                views.get_madrigal_videos, name='getMadrigalVideos'),
               url(r'^expNotes.py/$', 
                views.get_exp_notes, name='get_exp_notes'),
               url(r'^compareToArchive.py$', 
                views.compare_to_archive, name='compare_to_archive'),
               url(r'^openmadrigal/$', 
                views.open_madrigal, name='open_madrigal'),
               url(r'^madrigalAdmin/$', 
                views.madrigal_admin, name='madrigal_admin'),
               url(r'^madrigalDownload/$', 
                views.madrigal_download, name='madrigal_download'),
            ]