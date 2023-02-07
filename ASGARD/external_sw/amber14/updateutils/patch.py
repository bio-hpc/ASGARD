"""
This module contains the main patch classes that are used as bug fixes (and
occasionally upgrades) to Amber and AmberTools
"""

from os import path, makedirs, chmod, unlink, stat, rename
import re
import shutil
from subprocess import Popen, PIPE
import sys
from updateutils.exceptions import (PatchNotFound, NoPatchHeader, IllegalFile,
          BadPatchError, MissingProgram, PatchingError, UpdaterTypeError,
          PatchingWarning, PatchingUpdater)
from updateutils.utils import which
import warnings

class LocalPatch(object):
   """
   Basic Patch class. This type is specifically for ASCII-type patches which are
   stored locally on the machine, but should be sub-classed for different types
   of patches
   """

   def __init__(self, patch_file, prefix='bugfix'):
      """ Saves the file name and loads the header information """
      self.name = patch_file
      self.d2uconvert = False # Do we need to do a dos2unix conversion for cygwin?
      try:
         self.number = int(path.split(self.name)[1].replace(prefix + '.', ''))
      except (ValueError, TypeError):
         self.number = 0
      if not path.exists(patch_file):
         raise PatchNotFound('Could not open local patch file %s' % self.name)
      # Regular expression for the end of the header.
      end_line = re.compile(r'-+')
      header_ended = False
      # Load the header
      self.header = ''
      tmpfile = open(patch_file, 'r')
      for line in tmpfile:
         # Make sure the end_line matches and it's just a line of -'s
         if end_line.match(line.strip()) and not end_line.sub('', line.strip()):
            tmpfile.close()
            header_ended = True
            break
         self.header += line

      if not header_ended:
         warnings.warn('Patch %s does not have typical header information' %
                        self.name, NoPatchHeader)
         self.header = ''
      self.size = stat(self.name)

   def description(self):
      """ Glean the description out of the header """
      if not hasattr(self, 'header') or not self.header:
         return 'No header information in patch %s' % self.name
      if not 'Description' in self.header:
         return 'Could not find a valid description in %s' % self.name
      return self.header[self.header.index('Description'):]

   def author(self):
      """ Return the author from the header file """
      return self._get_attr(re.compile(r' *[Aa]uthor\(*s*\)*:*'))

   def date(self):
      """ Get the date the patch was made """
      return self._get_attr(re.compile(r' *[Dd]ate:* *'))

   def programs(self):
      """ Get the program edited from the header """
      progs = self._get_attr(re.compile(r' *[Pp]rogram\(*s*\)*:*'))
      if progs:
         return progs.strip().replace(',', ' ').split()
      return None

   def _get_attr(self, regexp):
      """ Generic method to retrieve an attribute from the header """
      if not hasattr(self, 'header') or not self.header:
         return None
      # Search for the given attribute in the header
      for line in self.header.split('\n'):
         if regexp.match(line):
            return regexp.sub('', line).strip()
      return None

   def matches_self(self):
      """ Sees if we match ourself """
      if not hasattr(self, '_matches_self'):
         self._files_edited()
      return self._matches_self

   def files_edited(self, key=None):
      """
      This is a wrapper around _files_edited to return a full list of files. It
      will take a key and optionally only return the files associated with the
      given key (or all files if key is None). Return is a single list
      """
      files = self._files_edited()
      if key is not None:
         return files[key]
      else:
         f = []
         # Skip the file modes
         for key in files:
            if key in ('new file modes', 'deleted file modes'): continue
            f += files[key]
         return f

   def _files_edited(self):
      """
      Find out which files will be edited. Returns a dictionary with the
      following keys referencing the following objects:
         
         'deleted files'  -- List of all files that are deleted in this patch
         'new files'      -- List of all files that are created in this patch
         'new file modes' -- List of permissions for all new files
         'modified files' -- List of all other files that are modified
      """
      # Prevent repeat parsing
      if hasattr(self, 'files'):
         return self.files
      # This regex matches lines that start like "+++ path/to/file.cpp"
      modfile = re.compile(r'\+\+\+ (\S+)')
      # See if we modify ourself in any way
      selffile = re.compile(r'\+\+\+ updater.py|\+\+\+ updateutils')
      # Track the old file names in case we delete a file, in which case the new
      # file name will be /dev/null (so we want the old file name instead)
      oldfile = re.compile(r'--- (\S+)')
      newfile = re.compile(r'new file mode (\d{3})+')
      delfile = re.compile('deleted file mode (\d{3})+')
      permbits = re.compile(r'(\d{3})')
      p = open(self.name, 'r')
      retdict = {'deleted files' : [], 'new files' : [],
                 'new file modes' : [], 'modified files' : [],
                 'deleted file modes' : []}
      oldfname = ''
      self._matches_self = False
      for line in p:
         if '\r' in line:
            self.d2uconvert = sys.platform == 'cygwin'
         if oldfile.match(line):
            oldfname = oldfile.match(line).groups()[0]
            continue
         if selffile.match(line):
            self._matches_self = True
            filename = line.lstrip('+++').strip().split()[0]
            if not filename in retdict['modified files']:
               retdict['modified files'].append(filename)
            continue
         elif modfile.match(line):
            filename = modfile.match(line).groups()[0]
            # Make sure our modifications do not have absolute paths and do not
            # travel outside AMBERHOME for security.
            if filename.startswith('../'):
               raise IllegalFile('Detected patched file outside AMBERHOME. Not '
                                 'allowed for security reasons [%s]' % filename)
            elif filename.startswith('/') and filename != '/dev/null':
               raise IllegalFile('Detected patched file with absolute path. '
                           'Not allowed for security reasons. [%s]' % filename)
            elif filename == '/dev/null':
               # Deleted files should only appear once!
               if (oldfname in retdict['new files'] or
                   oldfname in retdict['deleted files'] or
                   oldfname in retdict['modified files']):
                  raise BadPatchError('Deleted file appears multiple times.')
               retdict['deleted files'].append(oldfname)
            else:
               if oldfname == '/dev/null':
                  if (filename in retdict['new files'] or
                      filename in retdict['deleted files'] or
                      filename in retdict['modified files']):
                     raise BadPatchError('New file appears multiple times.')
                  retdict['new files'].append(filename)
                  # New files should appear only once!
               else:
                  # Make sure this is neither new or deleted
                  if (filename in retdict['new files'] or
                      filename in retdict['deleted files']):
                     raise BadPatchError('New/Deleted file appears multiple '
                                          'times.')
                  if not filename in retdict['modified files']:
                     retdict['modified files'].append(filename)
            oldfname = ''
            continue
         elif newfile.match(line):
            # oct code (necessary for chmod)
            permissions = int(permbits.findall(line).pop(), 8)
            retdict['new file modes'].append(permissions)
         elif delfile.match(line):
            # oct code (necessary for chmod)
            permissions = int(permbits.findall(line).pop(), 8)
            retdict['deleted file modes'].append(permissions)

      p.close()
      self.files = retdict
      return retdict

   def check_validity(self):
      """
      This checks that only AmberTools or only Amber are modified by this patch.
      This is primarily being added so developers can check their patch.
      """
      atools, amber = False, False
      for fname in self.files_edited():
         if fname.startswith('src/pmemd'):
            amber = True
         elif fname != '/dev/null':
            atools = True
      return atools != amber # this is identical to xor for booleans

   def restore_backups(self, remove_dir=True, reverse=False):
      """
      Restores all backup files
      """
      if not hasattr(self, 'backup_dest'):
         warnings.warn('No backup file directory to restore', PatchingWarning)
         return
      # If we're reversing the patch, restore deleted and modified. Otherwise,
      # restore new and modified
      if reverse:
         restores = self.files['modified files'] + self.files['new files']
      else:
         restores = self.files['modified files'] + self.files['deleted files']
      for fname in restores:
         shutil.copy(path.join(self.backup_dest, fname), fname)
      # Remove the backup directory if requested to
      if remove_dir: self.delete_backups()
   
   def delete_backups(self):
      """ Delete the backup directory """
      try:
         shutil.rmtree(self.backup_dest)
      except OSError:
         pass
      del self.backup_dest

   def check(self, reverse=False, ignore_prev_applied=True):
      """
      Runs patch with --dry-run so we don't do anything, but find out what it
      would do, anyway. Shortcut to apply with dryrun=True
      """
      return self.apply(reverse, ignore_prev_applied, dryrun=True)

   def apply(self, reverse=False, ignore_prev_applied=True, dryrun=False):
      """
      Apply the patch, in reverse if desired (to undo a patch). Return value is
      whether the patch succeeded or not.
      """
      # Create all new files. This is files['new files'] if reverse is False and
      # files['deleted files'] if reverse is True
      files = self._files_edited()
      if reverse:
         makefs = files['deleted files']
         rmfs = files['new files']
      else:
         makefs = files['new files']
         rmfs = files['deleted files']
      if not dryrun:
         for f in makefs:
            # Make sure the containing directory exists. If not, make it
            # (recursively)
            fpath, fname = path.split(f)
            if fpath and not path.isdir(fpath):
               makedirs(fpath)
            open(f, 'w').close()

      # Convert the file if need be
      if self.d2uconvert:
         p = open(self.name, 'r')
         p2 = open(self.name + '.tmp', 'w')
         for line in p:
            p2.write(line.replace('\r\n', '\n'))
         p.close(); p2.close()
         rename(self.name+'.tmp', self.name)
      p = open(self.name, 'r')
      patch = which('patch')
      if patch is None:
         raise MissingProgram('Cannot find patch. Install patch from your '
               'package manager (e.g., apt-get, yum, zypper, etc.) in order '
               'to apply updates')
      args = [patch, '-p0', '-N', '-s']
      if reverse:
         args += ['-R']
      if dryrun:
         args += ['--dry-run']

      process = Popen(args, stdin=p, stdout=PIPE, stderr=PIPE)

      (output, error) = process.communicate('')
      if process.wait() != 0:
         # See if we failed. If we failed, see if we failed badly (look for
         # FAILED), or see if we just ignored a previously applied hunk.
         output, error = output.decode('utf-8'), error.decode('utf-8')
         if 'FAILED' in output + error or not ignore_prev_applied:
            return False
         warnings.warn('Previously applied patch detected. Skipping.',
                       PatchingWarning)
      # Delete all deleted files. This is files['deleted files'] if reverse is
      # False and files['new files'] if reverse is True
      if not dryrun:
         # Delete any files that may be left over. But patch is _supposed_ to do
         # this for us, so squash any OSError's that are raised by trying to
         # remove a non-existent file
         for f in rmfs:
            try:
               unlink(f)
            except OSError:
               pass

      return True

   def _fix_modes(self, reverse=False):
      """
      Adjust the permissions of all added files if we were given that
      information in the patch file
      """
      if reverse:
         filekey = 'deleted files'
         permkey = 'deleted file modes'
      else:
         filekey = 'new files'
         permkey = 'new file modes'

      if len(self.files[filekey]) == 0:
         # Nothing to do
         return
      if len(self.files[permkey]) == 0:
         warnings.warn(('No permission information given in %s! All new files '
                       'will be read-only!') % self.name, PatchingWarning)
         return
      if len(self.files[filekey]) != len(self.files[permkey]):
         warnings.warn(('Permissions not found for every file in %s. All new '
                       'files will be read-only!') % self.name, PatchingWarning)
      # If we got this far, change all our permissions
      for i, fname in enumerate(self.files[filekey]):
         chmod(fname, self.files[permkey][i])

   def safe_apply(self, reverse=False, ignore_prev_applied=True):
      """
      Backs up the files, applies the patch, then either deletes the backups or
      restores them based on the result
      """
      
      succeeded = self.check(reverse, ignore_prev_applied)
      if not succeeded:
         raise PatchingError(('%s failed to apply. No changes made from this '
                              'patch') % self.name)
      # Disabled the locking and safe cleanup for now. This was originally
      # intended for the case where --dry-run did not work for patch, so a
      # failed patch could be effectively reversed. But since we use --dry-run,
      # this is only done to protect for interruptions.

      self.apply(reverse, ignore_prev_applied)
      self._fix_modes(reverse)

      # If we patch ourself, raise an exception to bail out of here so that
      # future patches can be applied with the 'updated' patcher (or
      # 'downgraded' if we're reversing)
      if self.matches_self():
         raise PatchingUpdater('Internal: Patching updater')

   # Compare patches on the basis of their numbers
   def __gt__(self, other): return self.number > other.number
   def __lt__(self, other): return self.number < other.number
   def __eq__(self, other): return self.number == other.number
   def __le__(self, other): return self.number <= other.number
   def __ge__(self, other): return self.number >= other.number
   def __ne__(self, other): return self.number != other.number
      
class CompressedLocalPatch(LocalPatch):
   """
   A Compressed version of an ASCII patch. Due to limitations in the compression
   modules in Python, the only way to reliably deal with compressed patches is
   to decompress them. This is a base class for all other compression
   """

   compressor = None
   extension = None
   args = []

   def __init__(self, patch_file, prefix='bugfix'):
      """ Unzips the file then calls the LocalPatch constructor """
      if not patch_file.endswith(self.extension):
         raise UpdaterTypeError('Bad extension for %s patch [%s]' %
                                (self.extension, patch_file))
      if not path.exists(patch_file):
         raise PatchNotFound('Could not open local patch file %s' % self.name)
      
      compressor = which(self.compressor)
      if compressor is None:
         raise MissingProgram(('Cannot find %s. Install compression software '
               'from your package manager (e.g., apt-get,  yum, zypper, etc.) '
               'in order to apply bzipped updates') % self.compressor)

      process = Popen([compressor] + self.args + [patch_file], 
                      stdout=PIPE, stderr=PIPE)

      output, error = process.communicate('')

      if process.wait():
         raise BadPatchError('Could not decompress %s with %s' %
                             (patch_file, self.compressor))
      LocalPatch.__init__(self, patch_file.replace(self.extension, ''), prefix)

class Bz2LocalPatch(CompressedLocalPatch):
   """ Compressed patch using bzip2 """
   compressor = 'bunzip2'
   extension = '.bz2'
   args = ['-f']

class GzipLocalPatch(CompressedLocalPatch):
   """ Compressed patch using gzip """
   compressor = 'gunzip'
   extension = '.gz'
   args = ['-f']

# List of all patch types with their required suffixes. For the time-being, add
# .bz2_ suffix since that was the temporary workaround for patch_amber.py. This
# can be removed when patch_amber.py is no longer relevant.
patch_suffixes = ['', '.bz2', '.gz', '.bz2_']
patch_types = {'' : LocalPatch, '.bz2' : Bz2LocalPatch, '.gz' : GzipLocalPatch,
               '.bz2_' : Bz2LocalPatch}

class RemotePatch(LocalPatch):
   """ Get the description from a remote ASCII patch """
   def __init__(self, open_url):
      """ Parses out the header, and that's it """
      # Regular expression for the end of the header.
      end_line = re.compile(r'-+')
      sizere = re.compile(r'Content-Length: (\d+)')
      self.size = int(sizere.findall(str(open_url.info()))[0])
      header_ended = False
      # Load the header
      self.header = ''
      for line in open_url:
         # Make sure the end_line matches and it's just a line of -'s
         if end_line.match(line.decode('utf-8').strip()) and \
            not end_line.sub('', line.decode('utf-8').strip()):
            header_ended = True
            break
         self.header += str(line.decode('utf-8'))

      if not header_ended:
         warnings.warn('Patch %s does not have typical header information' %
                        self.name, NoPatchHeader)
         self.header = ''

   # Kill the methods we don't want inherited
   def safe_apply(self, *args, **kwargs):
      raise PatchingError('Cannot apply an online patch.')

   def apply(self, *args, **kwargs):
      raise PatchingError('Cannot apply an online patch.')

   def files_edited(self, *args, **kwargs):
      raise PatchingError('Cannot parse a full online patch.')

   _files_edited = files_edited
   restore_backups = safe_apply
   delete_backups = safe_apply
   check = safe_apply
   _fix_modes = safe_apply

def test():
   """ Tests the patch implementation """
   from updateutils.downloader import urlopen
   print('Instantiating the local ASCII patch')
   patch1 = LocalPatch(path.join('test_patches', 'bugfix.1'))
   print('Instantiating the local Bzip2 patch')
   patch2 = Bz2LocalPatch(path.join('test_patches', 'bugfix.2.bz2'))
   print('Instantiating the local Gzip patch')
   patch3 = GzipLocalPatch(path.join('test_patches', 'bugfix.3.gz'))
   print('Done instantiating patches')

   print('\nPatch 1 statistics')
   print('Author: ' + patch1.author())
   print('Date:   ' + patch1.date())
   print('Programs: ' + ', '.join(patch1.programs()))
   print(patch1.description())
   print('New files: '+'\n           '.join(
                              patch1._files_edited()['new files']))
   print('Del files: '+'\n           '.join(
                              patch1._files_edited()['deleted files']))
   print('Others   : '+'\n           '.join(
                              patch1._files_edited()['modified files']))

   print('\nPatch 2 statistics')
   print('Author: ' + patch2.author())
   print('Date:   ' + patch2.date())
   print('Programs: ' + ', '.join(patch2.programs()))
   print(patch2.description())
   print('New files: '+'\n           '.join(
                              patch2._files_edited()['new files']))
   print('Del files: '+'\n           '.join(
                              patch2._files_edited()['deleted files']))
   print('Others   : '+'\n           '.join(
                              patch2._files_edited()['modified files']))

   print('\nPatch 3 statistics')
   print('Author: ' + patch3.author())
   print('Date:   ' + patch3.date())
   print('Programs: ' + ', '.join(patch3.programs()))
   print(patch3.description())
   print('New files: '+'\n           '.join(
                              patch3._files_edited()['new files']))
   print('Del files: '+'\n           '.join(
                              patch3._files_edited()['deleted files']))
   print('Others   : '+'\n           '.join(
                              patch3._files_edited()['modified files']))


   print('Applying patch 1')
   patch1.safe_apply()
   print('Applying patch 1 again strictly')
   try:
      patch1.safe_apply(ignore_prev_applied=False)
      print('FAILED TEST: unexpected success')
   except PatchingError:
      err = sys.exc_info()[1]
      print('Expected failure: %s: %s' % (type(err).__name__, err))

   print('Reversing patch 1')
   patch1.safe_apply(reverse=True)
   print('Reversing again, lax')
   try:
      patch1.safe_apply(reverse=True, ignore_prev_applied=True)
   except PatchingError:
      print('FAILED TEST: unexpected failure')

   print('Applying then reversing patch 2 and patch 3')
   patch2.safe_apply()
   patch3.safe_apply()
   patch2.safe_apply(reverse=True)
   patch3.safe_apply(reverse=True)

   print('\nTesting remote patch querying')
   remote = urlopen('http://ambermd.org/bugfixes/AmberTools/12.0/bugfix.1')
   patch4 = RemotePatch(remote)
   print('\nPatch 4 statistics')
   print('Author: ' + patch4.author())
   print('Date:   ' + patch4.date())
   print('Programs: ' + ', '.join(patch4.programs()))
   print(patch4.description())
   
   try:
      patch4._files_edited()
   except PatchingError:
      err = sys.exc_info()[1]
      print('Expected error: %s: %s' % (type(err).__name__, err))

   print('\nDone with patch tests.')
