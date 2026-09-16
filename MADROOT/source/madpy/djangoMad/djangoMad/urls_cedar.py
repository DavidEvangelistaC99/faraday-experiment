from django.conf.urls import include, url
from django.contrib import admin
import madweb.views

urlpatterns = [
    url(r'^', include('madweb.urls')),
    url(r'^', include('apps.cedar.urls')),
    url(r'^$', madweb.views.index),
]
