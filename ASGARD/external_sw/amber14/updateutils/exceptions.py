"""
List of exceptions and warnings used in the Amber updater script
"""

class AmberUpdaterError(Exception):
   """ Base error class for all Updater errors """

class AmberUpdaterWarning(Warning):
   """ Base warning class for all Updater warnings """

class ProgressBarError(AmberUpdaterError):
   """ A problem with one of the progress bars """

class LockWarning(AmberUpdaterWarning):
   """ A problem with instituting locks """

class RepoLocationError(AmberUpdaterError):
   """ If the location is invalid for the given PatchRepository class """

class PatchNotFound(AmberUpdaterError):
   """ If a given patch does not exist, but should """

class FileNotFound(AmberUpdaterError):
   """ Raised if a required file or folder does not exist """

class NoPatchHeader(AmberUpdaterWarning):
   """ If the patch file does not have a typical Amber header """

class IllegalFile(AmberUpdaterError):
   """ If we try to patch an illegal file (i.e., outside AMBERHOME) """

class BadPatchError(AmberUpdaterError):
   """ If the patch is no good for whatever reason """

class MissingProgram(AmberUpdaterError):
   """ If a necessary external program is missing """

class PatchingError(AmberUpdaterError):
   """ An error during the patch apply process """

class PatchingWarning(AmberUpdaterWarning):
   """ Nothing fatal, but worth reporting """

class UpdaterTypeError(AmberUpdaterError, TypeError):
   """ An error that occurs with using the wrong type variable here """

class APIError(AmberUpdaterError):
   """ Comes about from misusing the API """

class NoInternetAccess(AmberUpdaterError):
   """ If we cannot connect to the interwebs. """

class PatchingUpdater(AmberUpdaterError):
   """
   Raise this when we are patching this script so we can quit and allow future
   updates to occur with the 'fixed' version
   """

class NotEnoughPatches(AmberUpdaterError):
   """ Not enough patches have been found to perform some task """

class PreferenceError(AmberUpdaterError):
   """ Unknown preference """

class IllegalOptions(AmberUpdaterError):
   """ If users specify incompatible options on the command-line """

class UnknownArguments(AmberUpdaterWarning):
   """ If the user specified a command-line argument we don't recognize """

class NothingUpdated(AmberUpdaterWarning):
   """ If someone requested a specific update, but there was nothing to do """

class SetupError(AmberUpdaterError):
   """ If AMBERHOME is not set up correctly or updater.py has been moved """
