"""set up file for the Python Madrigal Remote API

$Id: setup.py 7057 2019-10-24 21:46:41Z brideout $
"""
import os, os.path, sys

from distutils.core import setup
    
setup(name="madrigalWeb",
      version="3.1.12",
      description="Remote Madrigal Python API",
      author="Bill Rideout",
      author_email="brideout@haystack.mit.edu",
      url="http://cedar.openmadrigal.org",
      packages=["madrigalWeb"],
      keywords = ['Madrigal'],
      scripts=['madrigalWeb/globalIsprint.py', 'madrigalWeb/globalDownload.py',
               'madrigalWeb/examples/exampleMadrigalWebServices.py']
      )

    