"""
Provides a class for reading CHARMM-style files. The key component to these
files is that the ! character is a comment character and everything after ! is
ignored.
"""
from chemistry.exceptions import CharmmFileError

class CharmmFile(object):
    """
    A CHARMM file that recognizes the "!" character as a 'comment' token. It
    can be iterated over and generally treated like a file object, but only
    spits out strings that have been truncated at its first comment character.
    
    There is currently no way to recognize a ! as a _non_ comment character,
    since allowing an escape character does not seem to be common practice and
    would likely introduce negative performance implications.
    """

    def __init__(self, fname, mode='r'):
        if mode not in ('r', 'w'):
            raise ValueError('Cannot open CharmmFile with mode "%s"' % mode)
        if mode == 'r':
            self.status = 'OLD'
        else:
            self.status = 'NEW'
        try:
            self._handle = open(fname, mode)
        except IOError, e:
            raise CharmmFileError(str(e))
        self.closed = False
        self.line_number = 0

    def write(self, *args, **kwargs):
        return self._handle.write(*args, **kwargs)

    def __iter__(self):
        # Iterate over the file
        for line in self._handle:
            try:
                idx = line.index('!')
                end = '\n'
            except ValueError:
                # There is no comment...
                idx = None
                end = ''
            yield line[:idx] + end

    def readline(self):
        self.line_number += 1
        line = self._handle.readline()
        try:
            idx = line.index('!')
            end = '\n'
        except ValueError:
            idx = None
            end = ''
        return line[:idx] + end

    def readlines(self):
        return [line for line in self]

    def read(self):
        return ''.join(self.readlines())

    def close(self):
        self._handle.close()
        self.closed = True

    def rewind(self):
        """ Return to the beginning of the file """
        self._handle.seek(0)

    def __del__(self):
        try:
            self.closed or self._handle.close()
        except AttributeError:
            # It didn't make it out of the constructor
            pass

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class CharmmStreamFile(object):
    """
    The stream file is broken down into sections of commands delimited by the
    strings:
        read <section> <options>
        ....
        ....
        end
    This object provides iterators over those sections and a file-like API for
    dealing with the text.

    """
    def __init__(self, fname):
        self.lines = CharmmFile(fname, 'r').readlines()
        self.line_number = 0

    def __iter__(self):
        return iter(self.lines)

    def rewind(self):
        """ Return to the beginning of the file """
        self.line_number = 0

    def next_section(self):
        """
        Fast-forwards the file to the next CHARMM command section

        Returns: (str, list)
            - The first string is the line defining the section that's being
              returned
            - The list is a list of all lines contained in the section
              excluding the "read <blah>" and "end" lines.

        Notes:
            The line pointer will be set to the line defining the 
        """
        lines = []
        while self.line_number < len(self.lines):
            line = self.lines[self.line_number].strip()
            if line[:4].lower() == 'read':
                title = line.strip()
                self.line_number += 1
                line = self.lines[self.line_number]
                while line and not line.strip().lower().startswith('end'):
                    lines.append(line)
                    self.line_number += 1
                    line = self.lines[self.line_number]
                return title, lines
            self.line_number += 1
        # No sections left
        return None, None


    def __del__(self): 

        pass
