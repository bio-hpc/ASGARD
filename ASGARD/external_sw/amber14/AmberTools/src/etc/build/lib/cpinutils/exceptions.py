"""
List of exceptions used in cpin creation
"""
import sys

def replace_excepthook(debug):
   " This function replaces sys.excepthook with one that suppresses tracebacks "
   def excepthook(exception_type, exception_value, tb):
      import traceback
      sys.stderr.write('%s: %s\n' % (exception_type.__name__, exception_value))
      sys.exit(1)
   if not debug:
      sys.excepthook = excepthook
