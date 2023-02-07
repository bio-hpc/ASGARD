"""
This module contains the classes necessary to display progress bars for various
downloads
"""
from __future__ import division
import sys
from updateutils.exceptions import ProgressBarError

class ProgressBar(object):
   """ Base progress bar class """

   def __init__(self, endvalue):
      """ Sets up the progress bar from a scale of 0 to endvalue """
      self.current_value = 0
      self.finished = False
      self.end_value = endvalue
      self._finalize_init()

   def update_add(self, value):
      """ Adds the passed value to the current value """
      if not isinstance(value, int) and not isinstance(value, float):
         raise ProgressBarError('update_add only takes a numeric value')
      if self.finished:
         raise ProgressBarError('ProgressBar already finished')
      self.current_value += value
      # Do not let us go over the max value
      self.current_value = max(self.current_value, 0)
      self.current_value = min(self.current_value, self.end_value)
      self._flush()

   def update_set(self, value):
      """ Sets the passed value to the current value """
      if not isinstance(value, int) and not isinstance(value, float):
         raise ProgressBarError('update_set only takes a numeric value')
      if self.finished:
         raise ProgressBarError('ProgressBar already finished')
      # Force the value to be between 0 and end value
      self.current_value = max(0, value)
      self.current_value = min(self.current_value, self.end_value)
      self._flush()

   def finish(self):
      """ Force finish """
      self.update_set(self.end_value)

   def reset(self):
      """ Resets the progress bar to be re-used """
      self.current_value = 0
      self.finished = False
      self._finalize_reset()

   def _flush(self):
      """
      Virtual method responsible for flushing the progress bar to the screen
      """
      raise ProgressBarError('Virtual method! Must implement _flush')

   def _finalize_reset(self):
      """ Virtual method responsible for finalizing the reset process """
      raise ProgressBarError('Virtual Method! Must implement _finalize_reset')

   def _finalize_init(self):
      """
      Does nothing in this class. Children can implement an 'extended'
      constructor this way
      """
      pass

class AsciiProgressBar(ProgressBar):
   output = sys.stdout
   BAR_SIZE = 50

   def _finalize_init(self):
      """ End of the Ascii Progress bar constructor """
      self.output.write('Downloading: [' + ' ' * self.BAR_SIZE + ']   0.0%')
      self.output.flush()
   
   def _flush(self):
      """ Flush the buffer to the screen """
      frac_done = self.current_value / self.end_value
      ndraw = int(frac_done * self.BAR_SIZE)
      # Back up to the beginning of the line and rewrite the whole progress bar
      self.output.write('\b' * (self.BAR_SIZE + 22))
      self.output.write('Downloading: [' + ':' * ndraw +
            ' ' * (self.BAR_SIZE - ndraw) + '] %5.1f%%' % (frac_done * 100))
      self.output.flush()

      # Print that we've finished
      if frac_done == 1:
         self.finished = True
         self.output.write(' Done.\n')
         self.output.flush()

   def _finalize_reset(self):
      """ Finalizes resetting the progress bar """
      self._finalize_init()


def test():
   """ Test suite """
   from time import sleep
   progress = AsciiProgressBar(1024)
   while not progress.finished:
      sleep(0.02)
      progress.update_add(8)
   progress.reset()
   while not progress.finished:
      sleep(0.02)
      progress.update_add(10)
   progress.reset()
   while not progress.finished:
      sleep(0.7)
      progress.update_add(128)
