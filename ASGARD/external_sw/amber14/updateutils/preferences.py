"""
This module contains the functionality required to write, parse, and add to
local preferences that allow users to customize various things. The currently
allowable things you may customize are:

   proxy -- This is the proxy address
   proxy user -- This is the user-name for the proxy (default is `whoami`)
   amber updates -- This is the location of the Amber bug fixes (if you want it
                    to be different than the Amber bugfix webpage)
   ambertools updates -- This is the location of the AmberTools bug fixes
"""

from copy import deepcopy
import re
from updateutils.exceptions import PreferenceError

# Some useful regular expressions for extracting the desired information
proxyre = re.compile(r'proxy = (\S+)', re.IGNORECASE)
proxyuserre = re.compile(r'proxy user = (\S+)', re.IGNORECASE)
amberre = re.compile(r'amber updates = (\S+)', re.IGNORECASE)
atoolsre = re.compile(r'ambertools updates = (\S+)', re.IGNORECASE)
domainre = re.compile(r'internet check = (\S+)', re.IGNORECASE)
      
def _value_or_none(val):
   """ Return None if val evaluates to False (e.g., an empty list) """
   if not val: return None
   if isinstance(val, list):
      return val[0]
   return val

def _write_ifnot_none(dest, key, indict):
   r"""
   Writes "key = val\n" to the destination if that key in the given dictionary
   is not None. Otherwise do nothing
   """
   val = indict[key]
   if val is not None:
      dest.write('%s = %s\n' % (key, val))

class Preferences(object):
   """ This is the class responsible for loading up the preferences """
   def __init__(self, filename='.updaterrc'):
      """
      This sets the file name and loads the preferences stored in filename
      """
      self.filename = filename
      self._preferences = {'proxy' : None, 'proxy user' : None,
            'amber updates' : None, 'ambertools updates' : None,
            'internet check' : None}
      self.load()
      self._originals = deepcopy(self._preferences)

   def load(self):
      """ Loads the preferences """
      global proxyre, proxyuserre, amberre, atoolsre
      try:
         text = open(self.filename, 'r').read()
      except IOError:
         # File does not exist, nothing to do
         return

      proxy = proxyre.findall(text)
      proxyuser = proxyuserre.findall(text)
      amber = amberre.findall(text)
      atools = atoolsre.findall(text)
      domain = domainre.findall(text)

      self._preferences['proxy'] = _value_or_none(proxy)
      self._preferences['proxy user'] = _value_or_none(proxyuser)
      self._preferences['amber updates'] = _value_or_none(amber)
      self._preferences['ambertools updates'] = _value_or_none(atools)
      self._preferences['internet check'] = _value_or_none(domain)

   def save(self):
      """ Saves the current preferences to the file """
      # Only save the preferences if they have been changed
      for key in self._preferences:
         if self._preferences[key] != self._originals[key]:
            break
      else:
         # If we reached here, everything is the same, so return without writing
         return
      # If we did not execute the 'else', write our file
      f = open(self.filename, 'w')
      write = lambda x: _write_ifnot_none(f, x, self._preferences)
      write('proxy')
      write('proxy user')
      write('amber updates')
      write('ambertools updates')
      write('internet check')
      f.close()

   def __getitem__(self, thing):
      """ Have attempts to access instances here pass through to preferences """
      return self._preferences[thing]
   
   def __setitem__(self, thing, val):
      """ Set preferences directly """
      if not thing in self._preferences.keys():
         raise PreferenceError('Unknown preference: %s' % thing)
      self._preferences[thing] = val

   def keys(self):
      return self._preferences.keys()

def test():
   """ Test suite for preferences """
   import os
   if os.path.exists('.updaterrc'): os.unlink('.updaterrc')
   prefs = Preferences('.updaterrc')
   prefs.load()
   failed = 0
   # Make sure everything is None
   passed = True
   for key in prefs.keys():
      if prefs[key] is not None:
         passed = False
   if passed:
      print('No preferences loaded yet: PASSED')
   else:
      print('Preferences loaded from empty file: FAILED')
      failed += 1
   
   prefs['proxy'] = 'http://this.that.com:800'
   prefs['proxy user'] = 'jswails'
   prefs['ambertools updates'] = \
     'http://www.clas.ufl.edu/users/jswails1/bugfixes/AmberTools/12.0/'
   prefs['amber updates'] = \
     'http://www.clas.ufl.edu/users/jswails1/bugfixes/12.0/'
   
   try:
      prefs['no option'] = 'nonsense'
      print('Allowed adding illegal preference. FAILED')
      failed += 1
   except PreferenceError:
      print('Blocked illegal preference: PASSED')

   prefs.save()
   if os.path.exists('.updaterrc'):
      print('Preference file saved: PASSED')
   else:
      print('Preference file not saved: FAILED')
      failed += 1

   # Try loading them up again
   prefs = Preferences('.updaterrc')
   prefs.load()
   passed = True
   if prefs['proxy'] != 'http://this.that.com:800':
      passed = False
      print('Proxy not parsed successfully: FAILED')
      print(prefs['proxy'])
      failed += 1
   if prefs['proxy user'] != 'jswails':
      passed = False
      print('Proxy user not parsed successfully: FAILED')
      print(prefs['proxy user'])
      failed += 1
   if prefs['ambertools updates'] != \
           'http://www.clas.ufl.edu/users/jswails1/bugfixes/AmberTools/12.0/':
      print(prefs['ambertools updates'])
      passed = False
      print('AmberTools location not parsed successfully: FAILED')
      failed += 1
   if prefs['amber updates'] != \
           'http://www.clas.ufl.edu/users/jswails1/bugfixes/12.0/':
      passed = False
      print('Amber location not parsed successfully: FAILED')
      failed += 1
      print(prefs['amber updates'])

   if passed:
      print('Preferences loaded successfully: PASSED')

   print('%d Preference tests failed' % failed)
