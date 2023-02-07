"""
Contains various utilities that don't really fit in anywhere else
"""
import os as _os
import sys as _sys

_asciiart = r"""       ,___,        ,___,   ,___,  ,_____,    ,_______,  ,_____,
      /     \       |    \./    |  |      \   |  ,____|  |      \
     /   ,   \      |           |  |  <>   |  |  |       |  <>   |
    /   / \   \     |  |\   /|  |  |      /   |  |____,  |      /
   /   /___\   \    |  | \./ |  |  |     |    |  ,____|  |     |
  /   ______,   \   |  |     |  |  |      \   |  |       |  ,,  \
 /   /       \   \  |  |     |  |  |  <>   |  |  |____,  |  | \  \
/___/         \___\ |__|     |__|  |______/   |_______|  |__|  \__\

                       <<< %s >>>

                       Version: %s
                       By: %s

"""

def which(program):
   """ Looks for a program in the PATH """
   def is_exe(prog):
      return _os.path.exists(prog) and _os.access(prog, _os.X_OK)

   # See if program has a full path. If it is, do not search PATH
   fpath, fname = _os.path.split(program)
   if fpath:
      return is_exe(program)

   for d in _os.environ['PATH'].split(_os.pathsep):
      if is_exe(_os.path.join(d, program)):
         return _os.path.join(d, program)

   return None

def print_updater_version(app, dest=None, exit=None):
   """ Prints the version information for the updater """
   global _asciiart
   from updateutils import __version__, __author__

   retstr = _asciiart % (app.full_parser.get_prog_name(), 
                         __version__, __author__)
   if app.preferences['ambertools updates'] is not None:
      retstr += ('Downloading AmberTools updates from mirror [%s]\n' %
                 app.preferences['ambertools updates'])
   if app.preferences['amber updates'] is not None:
      retstr += ('Downloading Amber updates from mirror [%s]\n' %
                 app.preferences['amber updates'])
   if hasattr(dest, 'write'):
      dest.write(retstr + '\n')
   if exit is not None:
      _sys.exit(exit)
   return retstr

def any(iterable):
   """ To support Python 2.X which doesn't have "any" """
   for it in iterable:
      if it: return True
   return False

def all(iterable):
   """ To support Python 2.X which doesn't have "all" """
   for it in iterable:
      if not it: return False
   return True
