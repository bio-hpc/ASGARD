"""
This module contains the functionality required to upgrade to the next major
version of Amber/AmberTools. Since I don't know what form that upgrade will
take, I will enforce a very flexible one now that will be a little kludgey.

Here is how it will work:

An "upgrade" folder will be present on the Amber web with three files:

   message.txt
   upgrade.patch.bz2
   upgrade
 
message.txt -- a text file with the message that will be displayed to the user
     after applying the upgrade. It should also contain the last patch that
     should be applied before updater.patch is applied for both AmberTools and
     Amber. It will then apply up to that patch -or- reverse back to that patch
     for both Amber and AmberTools (or only one if the other is not listed).
     This file should have the format:

BEGIN DESCRIPTION
put a small description here. Keep it short and 1-liner
END DESCRIPTION
BEGIN MESSAGE
put your message here, as you want it formatted
END MESSAGE
AmberTools patch level = #
Amber patch level = #

upgrade.patch -- a patch applied to anything necessary that is required to bring
     it up-to-date with the newest version. (Can be gzipped or bzipped)

upgrade -- a script of any kind that will be downloaded to AMBERHOME, given
     execute permissions, and then run to complete the update. If it is
     ultimately unnecessary, it can simply be a no-op. The script will be tested
     for successful or unsuccessful application, so please make sure that it
     returns an appropriate error code. The upgrade file will be removed after
     it is run. It is run from AMBERHOME.
"""

import os
import re
from subprocess import Popen
import sys
from updateutils.downloader import (needs_internet, DownloadManager, urlopen,
                                    HTTPError)
from updateutils.exceptions import NotEnoughPatches
from updateutils.patch import Bz2LocalPatch

class Upgrader(object):
   output = sys.stdout
   """
   This is responsible for testing to see if an upgrade is available and
   providing the machinery to do it

   Carry a variable "_found_upgrade" that is None if it hasn't been searched
   yet, True if it exists, and False if it does not exist
   """
   def __init__(self, location):
      """ Set up where online we look for an upgrade """
      self.location = location
      # Make sure we end with a /
      if not self.location.endswith('/'):
         self.location += '/'
      self._found_upgrade = None

   @needs_internet
   def look_for_upgrade(self):
      """ Looks online for message.txt, upgrade.patch, and upgrade """
      if self._found_upgrade is not None:
         return self._found_upgrade
      try:
         urlopen(self.location + 'message.txt')
         urlopen(self.location + 'upgrade.patch.bz2')
         urlopen(self.location + 'upgrade')
         self._found_upgrade = True
      except HTTPError:
         self._found_upgrade = False

      return self._found_upgrade

   @needs_internet
   def load_upgrade_info(self):
      """ Parses the message.txt to figure out how we perform the upgrade """
      if not self.look_for_upgrade():
         # Nothing to do if we have no upgrade
         return
      msg = urlopen(self.location + 'message.txt')
      msgre = re.compile(r'BEGIN MESSAGE\s+(.*)\s+END MESSAGE', re.DOTALL)
      descre = re.compile(r'BEGIN DESCRIPTION\s+(.*)\s+END DESCRIPTION',
                          re.DOTALL)
      atoolsre = re.compile(r'AmberTools patch level = (\d+)')
      amberre = re.compile(r'Amber patch level = (\d+)')

      msgtext = msg.read()
      self.message = msgre.findall(msgtext)[0].strip()
      self.description = descre.findall(msgtext)[0].strip()
      try:
         self.atools_patch = int(atoolsre.findall(msgtext)[0])
      except IndexError:
         self.atools_patch = None
      try:
         self.amber_patch = int(amberre.findall(msgtext)[0])
      except IndexError:
         self.amber_patch = None

   def _set_to_patch(self, onlinerepo, unappliedrepo, appliedrepo, num):
      """
      Applies or reverses the given patches to a given point.
      
      If more patches have been applied than the number of patches we want, we
      reverse the extra patches and purge the unapplied repo.

      If fewer patches have been applied than the number of patches we want, we
      download and apply the necessary patches.
      """
      # If we are already at the correct number, nothing to do. Just return
      if appliedrepo.last() == num:
         return
      # If we are ahead, rewind and purge. Then return
      if appliedrepo.last() > num:
         appliedrepo.apply(num+1)
         unappliedrepo.empty_repo()
         return
      # If we are here, then we are behind. Download and apply all of the
      # necessary patches from the remote repository. First we have to check
      # that our remote repo has enough bug fixes, though.
      if onlinerepo.last() < num:
         raise NotEnoughPatches('Not enough online patches for the upgrade')
      onlinerepo.download_patches(first=appliedrepo.last()+1,
                                  last=num, apply_=True)

   def set_ambertools_patch(self, onlinerepo, unappliedrepo, appliedrepo):
      """ Sets the AmberTools branch to the appropriate patch level """
      # If we have nothing to do in this branch, do nothing
      if self.atools_patch is None:
         return
      self._set_to_patch(onlinerepo, unappliedrepo, appliedrepo,
                         self.atools_patch)

   def set_amber_patch(self, onlinerepo, unappliedrepo, appliedrepo):
      """ Sets the AmberTools branch to the appropriate patch level """
      # If we have nothing to do in this branch, do nothing
      if self.amber_patch is None:
         return
      self._set_to_patch(onlinerepo, unappliedrepo, appliedrepo,
                         self.amber_patch)

   @needs_internet
   def _apply_upgrade_patch(self):
      """
      This will apply the upgrade.patch file. No check is done to make sure that
      the patch levels are up-to-date, so this check must be done beforehand.

      This method downloads upgrade.patch to AMBERHOME, applies it, then deletes
      it
      """
      dl = DownloadManager(self.location + 'upgrade.patch.bz2',
                           'upgrade.patch.bz2')
      dl.download_to('upgrade.patch.bz2')
      patch = Bz2LocalPatch('upgrade.patch.bz2')
      # Check that the patch will work first. Then apply. Then delete
      self.output.write('Attempting the upgrade. Please wait...\n')
      if patch.apply(reverse=False, ignore_prev_applied=True, dryrun=True):
         patch.apply(reverse=False, ignore_prev_applied=True, dryrun=False)
         patch._fix_modes(reverse=False)
         os.unlink('upgrade.patch')
         return True
      else:
         self.output.write('Upgrade patch failed. Upgrade not performed.\n')
         return False

   def _show_message(self):
      """ Shows the message we want to display to users """
      self.output.write(self.message + 
                        '\n\nNOTE: This upgrade cannot be reversed!\n')

   @needs_internet
   def _post_op(self):
      """ This downloads the upgrade script and runs it """
      self.output.write('Finalizing the upgrade\n')
      dl = DownloadManager(self.location + 'upgrade', 'Upgrade Finalizer')
      dl.download_to('upgrade')
      os.chmod('upgrade', int('755', 8))
      process = Popen(['./upgrade'])
      if process.wait():
         self.output.write('Warning: Possible error detected in upgrade...\n')
      os.unlink('upgrade')

   def perform_upgrade(self):
      """
      This actually performs the upgrade. Must be called AFTER the patch levels
      are set appropriately
      """
      self._show_message()
      if self._apply_upgrade_patch():
         self._post_op()
