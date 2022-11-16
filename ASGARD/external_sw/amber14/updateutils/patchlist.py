"""
This contains classes useful for storing lists of patches.
"""
from os import path
import re
from updateutils.downloader import needs_internet, urlopen
from updateutils.exceptions import UpdaterTypeError, PatchingError
from updateutils.patch import patch_suffixes, patch_types, RemotePatch

class PatchList(object):
   """ A list of patches """

   def __init__(self, prefix, order_matters=True):
      """
      Sets a list of bug fixes with the stored prefix. It also stores a list of
      classes that need to be used to instantiate each patch. Also keep track of
      all of our patch instances, so that we only instantiate them when we need
      to, but the instantiations are cached for quick access.

      Also see if we want to enforce that all bug fixes be added in order
      """
      global patch_suffixes
      self.prefix = prefix
      self.classes = []
      self.names = []
      self.numbers = []
      self.instances = []
      self.sizes = []
      self.order_matters = order_matters
      self.disordered = False
      self.patch_suffixes = patch_suffixes

   def __contains__(self, thing):
      """
      See if the bugfix we were sent is in our list already. Depending on what
      "thing" is, we will do something different:
         
         isinstance(thing, int) -- See if that bugfix # is in self.numbers

         thing.startswith(self.prefix) -- Extract the number and see if that is
                                          in self.numbers
         
         otherwise return thing in self.names
      """
      if isinstance(thing, int):
         return thing in self.numbers
      # If it is not an int, treat it like a string
      if self.prefix is not None and thing.startswith(self.prefix):
         patchnumre = re.compile(r'%s.(\d+)(?:\.bz2|\.gz)?' % self.prefix)
         if patchnumre.match(thing) and not patchnumre.sub('', thing):
            return int(patchnumre.match(thing).groups()[0]) in self.numbers
      elif self.prefix is None:
         # If prefix is not set, see if the file name is in the directory
         for fname in self.names:
            if path.split(fname)[1] == thing:
               return True
      # Otherwise, see if it's in the names. This will most likely be false,
      # unless they are the same file.
      return thing in self.names

   def delete(self, thing):
      """
      This will delete a given patch if self.__contains__(thing)
      """
      if not thing in self: return
      # Now we know the patch is present. Find it. And destroy it.
      if isinstance(thing, int):
         del self[self.numbers.index(thing)]
      elif self.prefix is not None and thing.startswith(self.prefix):
         patchnumre = re.compile(r'%s.(\d+)(?:\.bz2|\.gz)?' % self.prefix)
         if patchnumre.match(thing) and not patchnumre.sub('', thing):
            num = int(patchnumre.match(thing).groups()[0])
            del self[self.numbers.index(num)]
         else:
            del self[self.names.index(thing)]
      else:
         del self[self.names.index(thing)]

   def index(self, thing):
      """ Finds the index of the requested patch """
      if isinstance(thing, int):
         # Look for this patch number
         return self.numbers.index(thing)
      if isinstance(thing, str):
         if self.prefix is not None and thing.startswith(self.prefix):
            patchnumre = re.compile(r'%s.(\d+)(?:\.bz2|\.gz)?' % self.prefix)
            if patchnumre.match(thing) and not patchnumre.sub('', thing):
               num = int(patchnumre.match(thing).groups()[0])
               return self.numbers.index(num)
            else:
               return self.names.index(thing)
         elif self.prefix is None:
            # See if file name matches anything in this directory
            for i, fname in enumerate(self.names):
               if path.split(fname)[1] == thing:
                  return i
         return self.names.index(thing)

      raise ValueError('%s cannot be found in PatchList' % thing)

   def get_size(self, idx):
      """ Need a getter here because we have to instantiate the patch first """
      self._instantiate(self, idx)
      return self.sizes[idx]

   def __delitem__(self, i):
      """ Deletes the desired bugfix """
      del self.names[i], self.numbers[i], self.instances[i], self.classes[i]

   def __getitem__(self, i):
      """ PatchList[i] returns the name of the i'th patch """
      return self.names[i]

   def __iter__(self):
      """ Loop through all names """
      for name in self.names:
         yield name

   def __len__(self):
      """ Find out how many patches we have """
      return len(self.names)

   def empty(self):
      """ Is our list of patches empty? """
      return len(self.names) == 0

   def last(self):
      """ Returns the last patch in this list """
      try:
         return max(self.numbers)
      except ValueError:
         return 0

   def first(self):
      """ Returns the first patch in this list """
      try:
         return min(self.numbers)
      except ValueError:
         return 0

   def append(self, fname):
      """
      Add a new patch to the list. Insist that the prefix is consistent (or
      None), that bug fixes are added in order (if order_matters), and examine
      the suffix. See if it ends with any of the available patch_suffixes, then
      add the opening class to the classes list, file name to the names list,
      and None to the instances list (which will be filled with instances as
      they are created
      """
      number = 0
      pclass = None
      mysfx = ''
      fpath, basefname = path.split(fname)

      if self.prefix is not None:
         if not basefname.startswith(self.prefix):
            raise UpdaterTypeError('%s does not have prefix %s' % (fname,
                                   self.prefix))

      for sfx in self.patch_suffixes:
         # Skip over ASCII, since it has no suffix.
         if not sfx: continue
         if basefname.endswith(sfx):
            mysfx = sfx
            pclass = self.patch_types[sfx]
            break
      if pclass is None:
         # No suffix found? Must be ASCII
         pclass = self.patch_types['']

      # To support the remote patch list, pclass being -1 means there is no
      # online variant of that patch type and it must be downloaded. So have
      # pclass just be None as a placeholder
      if pclass == -1: pclass = None
      try:
         n = basefname.replace(self.prefix + '.', '').replace(mysfx, '').strip()
         number = int(n)
         # If we already loaded this bugfix, bail out
         if number in self: return
      except ValueError:
         # Just skip any patch that does not correspond to a numbered patch
         return
      except TypeError:
         # This means self.prefix is None
         self.disordered = True

      if self.disordered:
         # Do not add the same patch twice
         if basefname in self or fname in self: return
         self.names.append(fname)
         self.instances.append(None)
         self.classes.append(pclass)
         self.numbers.append(-1)
         self.sizes.append(0)
         return
      # Otherwise, keep us in order
      for i, nbr in enumerate(self.numbers):
         if number < nbr:
            self.names.insert(i, fname)
            self.instances.insert(i, None)
            self.classes.insert(i, pclass)
            self.numbers.insert(i, number)
            self.sizes.insert(i, 0)
            return
      # If we get here, we want to append to the end of the list
      self.names.append(fname)
      self.instances.append(None)
      self.classes.append(pclass)
      self.numbers.append(number)
      self.sizes.append(0)

   def author(self, idx):
      """ Get the author of the requested patch """
      return self._info_wrapper(idx, 'author')

   def date(self, idx):
      """ Get the date of the requested patch """
      return self._info_wrapper(idx, 'date')

   def description(self, idx, fmt=False):
      """ Get the description of the requested patch """
      desc = self._info_wrapper(idx, 'description')
      if fmt:
         desc = desc.replace('Description:', '').replace('description:', '')
         desclines = desc.strip().split('\n')
         desc = ''
         for line in desclines:
            desc += '\t' + ' '.join(line.strip().split()) + '\n'
      return desc

   def programs(self, idx):
      """ Get the description of the requested patch """
      return self._info_wrapper(idx, 'programs')

   def files_edited(self, idx):
      """ Get the list of edited files """
      return self._apply_wrapper(idx, 'files_edited')

   def matches_self(self, idx):
      return self._apply_wrapper(idx, 'matches_self')

   def _info_wrapper(self, idx, attr):
      """
      Wraps around one of the information-querying functions in the Patch
      classes. All this does is run the requested function.
      """
      return getattr(self._instantiate(idx), attr)()

   def _apply_wrapper(self, idx, attr):
      """
      Wraps around one of the patch-applying functions in the Patch classes. All
      this does is run the requested function, since all patch instances should
      implement these.
      """
      return getattr(self._instantiate(idx), attr)()


class LocalPatchList(PatchList):
   """ A list of patches on a local repository """

   def __init__(self, prefix, order_matters=True):
      """ Local patches have all the rights """
      global patch_types
      PatchList.__init__(self, prefix, order_matters)
      self.patch_types = patch_types

   def _instantiate(self, idx, force_new=False):
      """
      Instantiates a Patch object and returns it, and registers that instance in
      self.instances. If it has already been instantiated, just return that
      instance rather than create a new one
      """
      # Guard against an unavailable class type in this List
      if self.classes[idx] is None:
         return None

      if force_new or self.instances[idx] is None:
         p = self.classes[idx](self.names[idx], self.prefix)
         self.instances[idx] = p
         self.sizes[idx] = p.size

         # If our original version was gzipped or bzipped, then instantiating
         # unzips the file, so we need to change all of its attributes to ASCII
         # and strip the suffix
         if self.classes[idx] != self.patch_types['']:
            self.classes[idx] = self.patch_types['']
            for sfx in self.patch_suffixes:
               if sfx == '': continue
               self.names[idx] = self.names[idx].replace(sfx, '')

      return self.instances[idx]

   def apply(self, idx, reverse=False, ignore_prev_applied=True):
      """ Applies the requested patch (or all of them) """
      if idx == 'all':
         for idx in range(len(self.names)):
            if reverse:
               i = len(self.names) - idx - 1
            else:
               i = idx
            p = self._instantiate(i)
            p.safe_apply(reverse, ignore_prev_applied)
      else:
         p = self._instantiate(idx)
         return p.safe_apply(reverse, ignore_prev_applied)

class RemotePatchList(LocalPatchList):
   """ This is a list of patches on the remote """

   sizere = re.compile(r'Content-Length: (\d+)')

   def __init__(self, prefix, order_matters=True):
      """
      Settle for the LocalPatchList constructor, except overwrite the resulting
      patch_types dict for the remote types we can instantiate
      """
      PatchList.__init__(self, prefix, order_matters)
      # Now overwrite self.patch_types
      self.patch_types = {'' : RemotePatch, '.bz2' : -1, '.gz' : -1, 
                          '.bz2_' : -1}

   @needs_internet
   def _instantiate(self, idx, force_new=False):
      """
      Instantiates a Patch object and returns it, and registers that instance in
      self.instances. If it has already been instantiated, just return that
      instance rather than create a new one
      """
      # Guard against an unavailable class type in this List
      if self.classes[idx] is None:
         # First we need to get the size
         self.sizes[idx] = int(
                   self.sizere.findall(str(urlopen(self.names[idx]).info()))[0])
         return None

      if force_new or self.instances[idx] is None:
         urlobj = urlopen(self.names[idx])
         p = self.classes[idx](urlobj)
         self.sizes[idx] = p.size
         self.instances[idx] = p
         return p
      return self.instances[idx]

   def _apply_wrapper(self, idx, attr):
      """ No remote patch can do any applying/parsing """
      raise PatchingError('Patches must be downloaded before being applied')

   def _info_wrapper(self, idx, attr):
      """ Disable all of the methods decorated with this """
      try:
         return getattr(self._instantiate(idx), attr)()
      except AttributeError:
         return 'No information available'

def test():
   """ Test suite for the Patch lists """
   print('Initializing the local patch list and filling it')
   pl1 = LocalPatchList('bugfix')
   pl1.append('test_patches/bugfix.1')
   pl1.append('test_patches/bugfix.2.bz2')
   pl1.append('test_patches/bugfix.3.gz')

   print('Applying all three bug fixes')
   pl1.apply(0); pl1.apply(1); pl1.apply(2)
   print('Then reversing them')
   pl1.apply(0, True); pl1.apply(1, True); pl1.apply(2, True)

   print('Applying all three bug fixes again, using all')
   pl1.apply('all')
   print('Reversing all three bug fixes again, using all')
   pl1.apply('all', True)

   print('Initializing the remote patch list and filling it\n')
   pl2 = RemotePatchList('bugfix')
   basesite = 'http://www.clas.ufl.edu/users/jswails1/bugfixes/12.0'
   pl2.append(path.join(basesite, 'bugfix.1'))
   pl2.append(path.join(basesite, 'bugfix.2'))
   pl2.append(path.join(basesite, 'bugfix.3'))
   pl2.append(path.join(basesite, 'bugfix.4'))
   pl2.append(path.join(basesite, 'bugfix.5'))
   pl2.append(path.join(basesite, 'bugfix.6'))
   pl2.append(path.join(basesite, 'bugfix.7'))
   pl2.append(path.join(basesite, 'bugfix.8'))
   pl2.append(path.join(basesite, 'bugfix.9.bz2'))

   for i, p in enumerate(pl2):
      print('Patch %s:' % p)
      print('Author:   ' + pl2.author(i))
      print('Date:     ' + pl2.date(i))
      print('Programs: ' + pl2.programs(i))
      print(pl2.description(i))

   print('Done.')
