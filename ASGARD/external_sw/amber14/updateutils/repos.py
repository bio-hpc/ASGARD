"""
This module contains classes for controlling bugfix management in a database or
repository, either remote or local.

This allows different repositories to have different standard prefixes
"""
import os
import shutil
import sys
from updateutils.downloader import (needs_internet, DownloadManager, urlopen,
                                    HTTPError)
from updateutils.exceptions import (RepoLocationError, APIError, FileNotFound,
                                    UpdaterTypeError, PatchingUpdater)
from updateutils.patchlist import (LocalPatchList, RemotePatchList, PatchList,
                                  patch_suffixes)

class PatchRepository(object):
   """ Base class for storing/accessing folders of patches online """
   
   patch_list_cls = PatchList
   output = sys.stdout
   empty_repo_msg = 'No Patches Available'

   def __init__(self, location, title, prefix='bugfix'):
      """
      Get some details about the location of the repository and instantiate the
      patch list
      """
      self.location = location
      # Make sure location ends with '/' for directories and remotes
      if not self.location.endswith('/'):
         self.location += '/'
      self.prefix = prefix
      self.patch_list = self.patch_list_cls(self.prefix)
      self.reponame = title

   def patch_details(self, first=1, verbose=None):
      """
      Write the details of all of the patches according to the given verbosity
      level
      """
      # Determine what patches we want to print here
      patches_to_print = []
      if self.prefix is None:
         # If we have no prefix, just print them all!
         patches_to_print = range(len(self.patch_list.names))
      else:
         for i, num in enumerate(self.patch_list.numbers):
            if first <= num:
               patches_to_print.append(i)
      # Set our verbosity level based on how many patches we have to print
      self._set_verbose(verbose, len(patches_to_print))
      
      prbuf = ''
      for j, i in enumerate(patches_to_print):
         fname = os.path.split(self.patch_list[i])[1]
         prbuf += fname
         if self._to_print['program']:
            progs = self.patch_list.programs(i)
            if isinstance(progs, list):
               prbuf += ' (modifies %s)\n' % (', '.join(progs))
            else:
               prbuf += ' (modifies %s)\n' % progs
         if self._to_print['date'] and self._to_print['author']:
            prbuf += 'Released on %s (written by %s)\n' % (
                        self.patch_list.date(i), self.patch_list.author(i))
         elif self._to_print['date']:
            prbuf += 'Released on %s\n' % self.patch_list.date(i)
         elif self._to_print['author']:
            prbuf += 'Written by %s\n' % self.patch_list.author(i)
         if self._to_print['description']:
            prbuf += 'Description:\n' + self.patch_list.description(i,fmt=True)
            prbuf += '\n'
         if self._to_print['files']:
            prbuf += 'Files edited: %s\n' % ('\n              '.join(
                                    self.patch_list.files_edited(i)))
         # If we only print our bugfix number, do so cleanly and only print 10
         # bug fixes per line, to keep it neat
         if self._only_print_number:
            if j < len(patches_to_print) - 1:
               prbuf += ','
            if j % 10 == 9:
               prbuf += '\n'
            else:
               prbuf += ' '

      if not prbuf:
         return self.empty_repo_msg

      return prbuf
            
   def _set_verbose(self, verbose, npatches):
      """
      Sets the verbosity level. The default will vary depending on how many
      patches are present. There are 5 levels of 'verbosity', corresponding to:

         0  -  Print only the truncated base name
         1  -  Also print the program they modify
         2  -  Also print the description
         3  -  Also print the author and posting date
         4  -  Also print all of the files that are modified
      """
      if verbose is None:
         if npatches < 2:
            verbose = 4
         elif npatches < 4:
            verbose = 3
         elif npatches < 6:
            verbose = 2
         elif npatches < 11:
            verbose = 1
         else:
            verbose = 0
      else:
         # Make sure verbose is between 0 and 4
         verbose = min(int(verbose), 4)
         verbose = max(verbose, 0)
      self._to_print = {'program' : True, 'author' : True, 'date' : True,
                        'files' : True, 'description' : True}
      self._only_print_number = False
      if verbose < 4:
         self._to_print['files'] = False
      if verbose < 3:
         self._to_print['author'] = False
         self._to_print['date'] = False
      if verbose < 2:
         self._to_print['description'] = False
      if verbose < 1:
         self._to_print['program'] = False
         self._only_print_number = True

   def fill_list(self):
      """ Poll the location to see what patches are available """
      pass

   def last(self):
      """ Last patch in the repo """
      self.fill_list()
      return self.patch_list.last()

   def first(self):
      """ First patch in the repo """
      return self.patch_list.first()

   def empty(self):
      """ Is this patch repository empty? """
      return self.patch_list.empty()

class RemotePatchRepository(PatchRepository):
   """ A repository of updates/patches online """

   patch_list_cls = RemotePatchList
   empty_repo_msg = 'No patches available'

   def __init__(self, location, title, prefix='bugfix'):
      PatchRepository.__init__(self, location, title, prefix)
      self.filled = False
      if self.prefix is None:
         raise UpdaterTypeError("Remote repository must have a prefix!")
      if not (self.location.startswith('http://') or 
              self.location.startswith('ftp://')):
         raise RepoLocationError('RemotePatchRepository cannot point to a local'
                                ' directory')
   
   def set_destination_repo(self, other):
      """
      If a patch is downloaded, it is transferred to a "downloaded"
      repository where it is held until it is applied. This will register that
      destination repository with this remote repository
      """
      if not isinstance(other, LocalUnappliedPatchRepository):
         raise APIError('Can only download to a LocalUnappliedPatchRepository')
      if other.prefix != self.prefix:
         raise APIError('Can only download to a LocalUnappliedPatchRepository '
                        'with the same file name prefix')
      self.destination = other

   def patch_details(self, first=1, verbose=None):
      """
      Override this, because in order to get these details we need to fill our
      list first. For remote repositories this is put off as long as possible to
      avoid having to go to the interwebs
      """
      self.fill_list()
      return PatchRepository.patch_details(self, first, verbose)

   @needs_internet
   def fill_list(self):
      """
      This goes through the website looking for bugfixes in numerical order and
      fills the list with them. It then returns a reference to this list
      """
      global patch_suffixes
      # Do not fill a second time. It could take too long
      if self.filled: return
      self.filled = True
      num = 1
      exists = True
      while exists:
         found = False
         for sfx in patch_suffixes:
            url = '%s%s.%s%s' % (self.location, self.prefix, num, sfx)
            try:
               p = urlopen(url)
               p.close()
               found = True
               break
            except HTTPError:
               continue
         # If we didn't find
         if not found:
            exists = False
            break
         num += 1
         # If we did, add it to the list
         self.patch_list.append(url)

      return self.patch_list
   
   @needs_internet
   def download_patches(self, first=1, last=None, 
                        apply_=False, download_first=False):
      """
      This deals with downloading the patches into a new (local) patch
      repository, beginning with the `first' patch.
      
      If apply_ is True, it will apply the patches immediately after
      downloading them.
      
      If download_first is True, it will download all patches before trying to
      apply any of them. It is True simply because no patches are applied

      If last is None, download all the rest. Otherwise, stop at last
      """
      # Fill our list of patches if we haven't already
      if not self.filled:
         self.fill_list()
      # If our list is empty, bail out.
      if self.patch_list.empty():
         return
      
      if last is None:
         last = self.patch_list.last()

      if not hasattr(self, 'destination'):
         raise APIError('You must register a local destination for remote '
                        'repositories to download patches')
      
      ndl = 0 # Number we've downloaded
      first_dl = 0
      # Now loop through our list of patches to find out which ones we want to
      # download. We download from `first' to the end
      for i, num in enumerate(self.patch_list.numbers):

         if num >= first and num <= last:
            # Do not download patches already in the destination
            if num in self.destination.patch_list:
               # Simply apply it if we already have it
               if apply_ and not download_first:
                  self.destination.apply(num)
               continue
            if ndl == 0:
               first_dl = num
               self.output.write('Downloading updates for %s\n' % self.reponame)
            ndl += 1
            # Replace .bz2_ with .bz2 to account for patch_amber.py fix.
            pname = os.path.split(self.patch_list[i])[1].replace('.bz2_','.bz2')
            manager = DownloadManager(self.patch_list[i], '%s/%s.%d' %
                                      (self.reponame, self.prefix, num))
            manager.download_to(os.path.join(self.destination.location, pname))
            # Refresh the bugfix repository now that we added a patch
            self.destination.fill_list()
            if apply_ and not download_first:
               self.destination.apply(num)

      if ndl > 0 and apply_ and download_first:
         self.destination.apply(first_dl, first_dl + ndl - 1)

   def _set_verbose(self, value, npatches):
      """
      Override this method to force file printing off -- not available for
      remote patches
      """
      PatchRepository._set_verbose(self, value, npatches)
      self._to_print['files'] = False

   def last(self):
      """ Last patch in the repo """
      self.fill_list()
      return PatchRepository.last(self)

   def first(self):
      """ First patch in the repo """
      self.fill_list()
      return PatchRepository.first(self)

   def empty(self):
      """ Is this patch repository empty? """
      self.fill_list()
      return PatchRepository.empty(self)

class _LocalPatchRepository(PatchRepository):
   """
   A repository of updates/patches on a local filesystem. This is actually a
   base class for those repos that manage applied and unapplied repositories.
   """

   patch_list_cls = LocalPatchList
   empty_repo_msg = 'No patches in folder'

   def __init__(self, location, title, prefix='bugfix'):
      PatchRepository.__init__(self, location, title, prefix)
      if (self.location.startswith('http://') or
            self.location.startswith('ftp://')):
         raise RepoLocationError('LocalPatchRepository cannot point '
                                 'to remote repos')
      self.absolute = (self.location.startswith(os.pathsep) or 
                       self.location.startswith('file:///'))
      # Strip off the file:// to return just the directory
      if self.location.startswith('file:///'):
         self.location = self.location.replace('file://', '')

      try:
         self.fill_list()
      except OSError:
         # The directory does not exist. Create it. No point in filling the list
         # afterwards since there obviously won't be anything in it. For
         # security, though, DON'T create the directory if the path is absolute.
         if self.absolute:
            raise FileNotFound('Patch directory with absolute path does not '
                            'exist. I will not create it for security reasons.')
         os.makedirs(self.location)

   def fill_list(self):
      """
      Finds out which bug fixes are present in this repository and populates the
      bugfix list. Ignore everything except files that start with prefix.
      """
      self.patch_list = LocalPatchList(self.prefix)
      for fname in os.listdir(self.location):
         if self.prefix is not None and fname.startswith(self.prefix):
            self.patch_list.append(os.path.join(self.location, fname))
         elif self.prefix is None:
            # Add all files if there is no prefix
            self.patch_list.append(os.path.join(self.location, fname))

   def apply(self, first, last=None):
      """
      Applies all of the patches between the `first' and last. If last is not
      specified, we just apply the first and be done with it

      self.reverse is used to determine if we apply in reverse or not
      """
      if not hasattr(self, 'destination'):
         raise APIError('Must have a place to move the applied patch')
      
      # See if we want to apply all of them
      if first == 'all':
         first, last = self.patch_list.first(), self.patch_list.last()
      # Otherwise first must be an integer, and see if we only want to apply 1
      if last is None:
         # If we reverse, we need to reverse _back_ to the first. So set last to
         # the last available patch
         if self.reverse:
            last = self.patch_list.last()
         else:
            last = first
      # Make sure first is lowest and last is highest to make reversing easier
      first, last = min(first, last), max(first, last)
      # If we reverse our patches, apply them in reverse order
      if self.reverse:
         iterseq = reversed(range(len(self.patch_list)))
      else:
         iterseq = range(len(self.patch_list))
      for i in iterseq:
         # If we do not have an ordered patch, treat 'first' and 'last' like
         # indexes into patch_list, rather than patch numbers
         if self.patch_list.disordered:
            pname = os.path.split(self.patch_list[i])[1]
            if self.reverse:
               self.output.write('Reversing %s\n' % pname)
            else:
               self.output.write('Applying %s\n' % pname)
            try:
               self.patch_list.apply(i, reverse=self.reverse)
               patches_self = False
            except PatchingUpdater:
               patches_self = True
            # Now move the patch
            if self.absolute:
               shutil.copy(self.patch_list[i],
                           os.path.join(self.destination.location, pname))
            else:
               os.rename(self.patch_list[i],
                         os.path.join(self.destination.location, pname))
            # If we patched ourself, bail out
            if patches_self:
               raise PatchingUpdater('DETECTED A PATCH TO THE UPDATER SCRIPT. '
                        'I am quitting now so future updates are applied '
                        'correctly. Make sure to run this script again (until '
                        'you stop receiving this message).')
            continue
         # Otherwise we have an ordered patch repository, so treat the numbers
         # like patch numbers rather than indexes
         num = self.patch_list.numbers[i]
         if num >= first and num <= last:
            pname = os.path.split(self.patch_list[i])[1]
            # Inform the user
            self.output.write(self.actionstr % (self.reponame, pname) + '\n')
            # Apply and move the patch
            try:
               self.patch_list.apply(i, reverse=self.reverse)
               patches_self = False
            except PatchingUpdater:
               patches_self = True
            # Re-generate our pname, since any compression suffix has been
            # removed.
            pname = os.path.split(self.patch_list[i])[1]
            # If we have an absolute path, just copy. Otherwise, move
            if self.absolute:
               shutil.copy(self.patch_list[i],
                     os.path.join(self.destination.location, pname))
            else:
               os.rename(self.patch_list[i],
                      os.path.join(self.destination.location, pname))
            # If we patched ourself, bail out
            if patches_self:
               raise PatchingUpdater('DETECTED A PATCH TO THE UPDATER SCRIPT. '
                        'I am quitting now so future updates are applied '
                        'correctly. Make sure to run this script again (until '
                        'you stop receiving this message).')
      # Now that we're done, refresh the list of patches in each repo
      self.fill_list()
      self.destination.fill_list()

class LocalUnappliedPatchRepository(_LocalPatchRepository):
   """
   This is where unapplied patches are stored. Patches from these repositories
   are NOT reversed when applied.
   """
   reverse = False
   actionstr = 'Applying %s/%s'
   empty_repo_msg = 'No unapplied patches found'

   def empty_repo(self):
      """
      This empties the repository of all unapplied patches. Refuse to do this
      for a repository with an absolute path, as this is unsafe
      """
      if self.empty():
         return
      if self.absolute:
         raise APIError('Cannot delete files from a repository with an absolute'
                       ' path. This is a security risk!')
      for fname in os.listdir(self.location):
         if self.prefix is None or fname.startswith(self.prefix):
            os.remove(os.path.join(self.location, fname))
      # We have emptied our repo, so get the list again
      self.fill_list()

   def download_patches(self, first=1, last=None, 
                        apply_=False, download_first=False):
      """
      This implements the remote repository version of the patch downloader so
      that a LocalUnappliedPatchRepository can act like a remote repo. Instead
      of downloading, however, we just copy
      """
      # If our list is empty, bail out.
      if self.patch_list.empty():
         return

      if last is None:
         last = self.patch_list.last()

      if not hasattr(self, 'destination'):
         raise APIError('You must register a local destination to download '
                        'patches')
      
      ndl = 0 # Number we've downloaded
      first_dl = 0
      # Now loop through our list of patches to find out which ones we want to
      # download. We download from `first' to the end
      for i, num in enumerate(self.patch_list.numbers):

         if num >= first and num <= last:
            # Do not download patches already in the destination
            if num in self.destination.patch_list:
               # Simply apply it if we already have it
               if apply_ and not download_first:
                  self.destination.apply(num)
               continue
            if ndl == 0:
               first_dl = num
               self.output.write('Copying updates for %s\n' % self.reponame)
            ndl += 1
            # Replace .bz2_ with .bz2 to account for patch_amber.py fix.
            pname = os.path.split(self.patch_list[i])[1].replace('.bz2_','.bz2')
            shutil.copy(self.patch_list[i],
                        os.path.join(self.destination.location, pname))
            # Refresh the bugfix repository now that we added a patch
            self.destination.fill_list()
            if apply_ and not download_first:
               self.destination.apply(num)

      if ndl > 0 and apply_ and download_first:
         self.destination.apply(first_dl, first_dl + ndl - 1)

   def set_destination_repo(self, other):
      """
      If a patch is downloaded, it is transferred to a "downloaded"
      repository where it is held until it is applied. This will register that
      destination repository with this remote repository
      """
      if not isinstance(other, _LocalPatchRepository):
         raise APIError('Can only move to a _LocalPatchRepository')
      if other.prefix != self.prefix:
         raise APIError('Can only move to a LocalAppliedPatchRepository '
                        'with the same file name prefix')
      self.destination = other

class LocalAppliedPatchRepository(_LocalPatchRepository):
   """
   This is where applied patches are stored. Patches from these repositories ARE
   reversed when applied
   """
   reverse = True
   actionstr = 'Reversing %s/%s'
   empty_repo_msg = 'No applied patches found'

   def set_destination_repo(self, other):
      """
      If a patch is downloaded, it is transferred to a "downloaded"
      repository where it is held until it is applied. This will register that
      destination repository with this remote repository
      """
      if not isinstance(other, LocalUnappliedPatchRepository):
         raise APIError('Can only unapply to a LocalUnappliedPatchRepository')
      if other.prefix != self.prefix:
         raise APIError('Can only unapply to a LocalUnappliedPatchRepository '
                        'with the same file name prefix')
      self.destination = other

def remote_patches(location, title, prefix='bugfix'):
   """
   This sets up a remote patch location. Since the remote repo _could_ be a
   local mirror, we first try to instantiate a remote patch repository, but fall
   back on a local unapplied patch repository if that doesn't work
   """
   try:
      return RemotePatchRepository(location, title, prefix=prefix)
   except RepoLocationError:
      return LocalUnappliedPatchRepository(location, title, prefix=prefix)

def test():
   """ Test suite for all of the repositories """
   sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
   print('Setting up the repositories...')
   AT_remote = RemotePatchRepository(
            'http://www.clas.ufl.edu/users/jswails1/bugfixes/AmberTools/12.0/',
            'AmberTools 12', prefix='bugfix')
   AT_local_unap = LocalUnappliedPatchRepository(os.path.join('.patches',
            'AmberTools_Unapplied_Patches'), 'AmberTools 12', prefix='bugfix')

   AT_local_ap = LocalAppliedPatchRepository(os.path.join('.patches',
            'AmberTools_Applied_Patches'), 'AmberTools 12', prefix='bugfix')

   print('Registering destination repos')
   AT_remote.set_destination_repo(AT_local_unap)
   AT_local_unap.set_destination_repo(AT_local_ap)
   AT_local_ap.set_destination_repo(AT_local_unap)
   
   # Testing out the remote repository
   print('Polling the available fixes in the remote repository')
   AT_remote.fill_list()
   print('Found %d remote patches' % len(AT_remote.patch_list))

   print('Printing details of remote repo')
   print(AT_remote.patch_details())
   print('\nPrinting details of patches 22 and up')
   print(AT_remote.patch_details(first=22))
   print('\nPrinting details of patches 27 and up')
   print(AT_remote.patch_details(first=27))
   print('\nPrinting details of patches 31 and up')
   print(AT_remote.patch_details(first=31))

   # Test out the interchange between local unapplied and remote
   print('Downloading patch files into the unapplied repository')
   AT_remote.download_patches(first=1, apply_=False)
   print('Flushing the local unapplied repository')
   AT_local_unap.empty_repo()
   print('Downloading the patch files again, and applying immediately')
   AT_remote.download_patches(first=1, apply_=True, download_first=False)
   print('List unapplied patch details:')
   print(AT_local_unap.patch_details())
   print('List applied patch details:')
   print(AT_local_ap.patch_details())
   print('Reversing the patches')
   AT_local_ap.apply('all')
