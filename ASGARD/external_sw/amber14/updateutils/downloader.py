"""
This module contains the methods for downloading files from the internet. It
supports setting up a proxy to connect to the web and download the files.

It defines a number of limits that control how the process of downloading is
displayed to the user. These limits, all shown in bytes, are:

SHOW_PB_LIMIT     The smallest size for which a progress bar is shown

SHOW_BY_LIMIT     The largest size for which download size is shown sans units

SHOW_KB_LIMIT     The largest size for which download size is shown in KB

SHOW_MB_LIMIT     The largest size for which download size is shown in MB.
                  Anything larger is (scary!) shown in GB

DL_CHUNK_SIZE     Number of bytes we download at a time

Also defined are some constants to help with connecting to the internet. The
_WEB_ACCESS and _ACCESS_CODES variables are internal-use only. The DOMAIN_CHECK
is the website we should check for connection.  has_internet_access will check
for internet (unless it has been checked already).

needs_internet is a decorator that should adorn any function that requires
internet (except for those that won't be reachable until internet connection has
already been validated). This will execute the appropriate check and raise an
exception if no connection is detected before running the function as normal
"""
from __future__ import division
import re
import os
import sys
from updateutils.exceptions import NoInternetAccess
from updateutils.progressbar import AsciiProgressBar
try:
   from urllib2 import urlopen, addinfourl, HTTPError
except ImportError:
   from http.client import HTTPResponse as addinfourl
   from urllib.request import urlopen
   from urllib.error import HTTPError

# Limits
SHOW_PB_LIMIT = 102400
SHOW_BY_LIMIT = 1024
SHOW_KB_LIMIT = 1843200
SHOW_MB_LIMIT = 1073741824
DL_CHUNK_SIZE = 4096 # 4 KB chunks

_WEB_ACCESS = 'UNKNOWN'

_ACCESS_CODES = {'UNKNOWN' : None, 'YES' : True, 'NO' : False, 'DOWN' : False}

_NEEDS_PROXY = False
_PROXY_USER = None
_PROXY_ADDY = None

DOMAIN_CHECK = 'http://ambermd.org'

def set_proxy_info(address=None, user=None):
   """
   This sets the _NEEDS_PROXY variable with the given user name. If a proxy is
   necessary, then "has_internet_access" will request the proxy information as
   soon as internet access is required. This is triggered by the
   "needs_internet" decorator.
   """
   global _PROXY_USER, _NEEDS_PROXY, _PROXY_ADDY
   _NEEDS_PROXY = bool(address)
   _PROXY_ADDY = address
   _PROXY_USER = user

def has_internet_access(fatal=False):
   """ Checks to see if we have access to the web """
   global _WEB_ACCESS, _ACCESS_CODES, DOMAIN_CHECK, _NEEDS_PROXY
   if _ACCESS_CODES[_WEB_ACCESS] is None:
      # If we need to set up the proxy, do it now
      if _NEEDS_PROXY:
         setup_proxy(_PROXY_ADDY)
      try:
         host = urlopen(DOMAIN_CHECK)
         host.close()
         _WEB_ACCESS = 'YES'
      except IOError:
         err = sys.exc_info()[1]
         if 'timed out' in str(err):
            _WEB_ACCESS = 'DOWN'
         else:
            _WEB_ACCESS = 'NO'
   access = _ACCESS_CODES[_WEB_ACCESS]
   if fatal and not access:
      raise NoInternetAccess('Cannot connect to %s' % DOMAIN_CHECK)
      
   return _ACCESS_CODES[_WEB_ACCESS]

def needs_internet(fcn):
   """
   Apply this function as a decorator around all functions that require access
   to the internet.  It will check if we have internet access first, and if it
   does, then it will call the desired function. Otherwise, raise an exception
   complaining about the lack of connection
   """
   def new_fcn(*args, **kwargs):
      from updateutils.downloader import has_internet_access
      if has_internet_access(fatal=True):
         return fcn(*args, **kwargs)

   return new_fcn

class DownloadManager(object):
   """ Sets up downloading files """

   sizere = re.compile(r'Content-Length: (\d+)')
   output = sys.stdout

   def __init__(self, url, shortname=None):
      """
      Open up the URL, unless it is already an addinfourl instance, determine
      how big the download is, and from that determine how we will display the
      download
      """
      global SHOW_PB_LIMIT, SHOW_BY_LIMIT, SHOW_KB_LIMIT, SHOW_MB_LIMIT
      if isinstance(url, addinfourl):
         self.url = url
      else:
         self.url = urlopen(url)

      if shortname is None:
         self.shortname = self.url.url
      else:
         self.shortname = shortname

      self.bytesize = int(self.sizere.findall(str(self.url.info()))[0])
      self.show_progress = self.bytesize >= SHOW_PB_LIMIT

      # Determine how we show the user how big we are
      self.reported_size = '?? Unknown'
      if self.bytesize < SHOW_BY_LIMIT:
         self.reported_size = '%d Bytes' % self.bytesize
      elif self.bytesize < SHOW_KB_LIMIT:
         self.reported_size = '%.2f KB' % (self.bytesize / 1024)
      elif self.bytesize < SHOW_MB_LIMIT:
         self.reported_size = '%.2f MB' % (self.bytesize / 1024 / 1024)
      else:
         # YIKES!
         self.units = '%.2f GB(!)' % (self.bytesize / 1024 / 1024 / 1024)

   def download_to(self, destination):
      """ Downloads the file to a destination """
      global DL_CHUNK_SIZE
      # Set up the progress bar if we're reporting one
      self.output.write('Downloading %s (%s)\n' % 
                        (self.shortname, self.reported_size))
      self.output.flush()
      if self.show_progress:
         progress = AsciiProgressBar(self.bytesize)
      total_dled_bytes = 0
      dest = open(destination, 'wb')
      while total_dled_bytes < self.bytesize:
         dest.write(self.url.read(DL_CHUNK_SIZE))
         if self.show_progress:
            progress.update_add(DL_CHUNK_SIZE)
         total_dled_bytes += DL_CHUNK_SIZE
      dest.close()

def setup_proxy(address):
   """ Get the Proxy information from the user """
   global _PROXY_USER
   from getpass import getpass
   passwd = getpass("Proxy Password: ")
   # Strip out the "http:// from the proxy address
   address = address.replace('http://', '')
   # Add http proxy to environment variables
   if _PROXY_USER is None:
      os.environ['http_proxy'] = 'http://%s' % address
   else:
      if len(passwd) == 0:
         os.environ['http_proxy'] = 'http://%s@%s' % (_PROXY_USER, address)
      else:
         os.environ['http_proxy'] = 'http://%s:%s@%s' % (_PROXY_USER,
                                                         passwd, address)

def test():
   """ Test suite for the downloader """
   sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
   website = 'http://www.clas.ufl.edu/users/jswails1/bugfixes'
   dl1 = DownloadManager('%s/AmberTools/12.0/bugfix.19' % website, 'AT.19')
   dl2 = DownloadManager('%s/AmberTools/12.0/bugfix.26' % website)
   dl3 = DownloadManager('%s/12.0/bugfix.9.bz2' % website, 'Amber/bugfix.9')

   dl1.download_to('bugfix.18')
   dl2.download_to('bugfix.28')
   dl3.download_to('bugfix.9.bz2')
