#! PYTHONEXE

" This program will deconvolute pH REM cpout files into pH-based cpout files "

import re, os, sys
from optparse import OptionParser

class CpoutError(Exception):
    """ Raised for an error reading/parsing CPOUT file """

#-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

class RemCpout(object):
    """ cpout file from REM """

    #======================================================

    def __init__(self, fil):
        if type(fil).__name__ != 'file':
            raise TypeError('RemCpout: expected file type!')
        self.file = fil

    #======================================================

    def get_next_record(self):
        """ Gets the next record in this cpout file """
        lines = []
        line = self.file.readline()
        # Hit EOF
        if not line: return
        while line.strip():
            lines.append(line.strip())
            line = self.file.readline()
        return TitrationRecord(lines)

#-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

class TitrationRecord(object):
    """ A Titration Record from Cpout file """

    full_re = re.compile(r'Solvent pH: *([+-]?\d+\.\d+)')
    rec_re = re.compile(r'Residue *(\d+) State: *(\d+) pH: *([+-]?\d+\.\d+)')
    ph_rec = re.compile(r' *pH: *[+-]?\d+\.\d+')

    #======================================================

    def __init__(self, lines):
        self.lines = lines
        fullmatch = self.full_re.match(lines[0])
        if fullmatch:
            self.full = True
            self.pH = float(fullmatch.groups()[0])
        else:
            self.full = False
            recmatch = self.rec_re.match(lines[0])
            if not recmatch:
                raise CpoutError('Did not find expected line in cpout record!')
            self.pH = float(recmatch.groups()[2])

    #======================================================

    def write_to(self, outfile):
        """ Writes this record somewhere """
        if self.full:
            outfile.write(os.linesep.join(self.lines[:4]))
            outfile.write(os.linesep)
            outfile.write(os.linesep.join([self.ph_rec.sub('',l)
                                           for l in self.lines[4:]]))
            outfile.write(os.linesep)
            outfile.write(os.linesep)
            return
        outfile.write(os.linesep.join([self.ph_rec.sub('',l)
                                       for l in self.lines]))
        outfile.write(os.linesep)
        outfile.write(os.linesep)

#-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

def main():
    usage = '%prog [Options] cpout1 cpout2 cpout3 ... cpoutN'
    parser = OptionParser(usage=usage)
    parser.add_option('-p', '--prefix', dest='prefix', default=None,
                      help='Prefix of output cpout files. Output cpout files '
                      'are named <prefix>.pH_<pH>')
    opt, arg = parser.parse_args()

    try:
        cpout_list = [RemCpout(open(f, 'r')) for f in arg]
    except IOError:
        print 'Could not find file %s' % f
        sys.exit(1)

    if not opt.prefix:
        print 'Need a file prefix'
        parser.print_help()
        sys.exit(1)

    records = [cpout.get_next_record() for cpout in cpout_list]

    pH_list = [rec.pH for rec in records]

    files = {}
    for pH in pH_list:
        files[pH] = open(opt.prefix + '.pH_%.2f' % pH, 'w')

    while records[0]:
        for rec in records:
            if rec: rec.write_to(files[rec.pH])
        records = [cpout.get_next_record() for cpout in cpout_list]

    print 'Done!'


if __name__ == '__main__':
    import time
    start = time.time()
    main()
    print 'Took: %.1f sec' % ((time.time() - start))
