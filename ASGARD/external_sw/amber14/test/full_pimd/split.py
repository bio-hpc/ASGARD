from sys import argv
from string import atoi
from string import split

input = open( argv[1], "r" )
ncopy = atoi( argv[2] )
output = argv[3]

title = input.readline()

head  = input.readline()

natom = atoi( split( head )[0] )

natomCL = natom/ncopy

pos = []
nline = natom/2
for i in xrange(nline) :
    line = input.readline()
    items = split( line )
    for item in items: pos.append( item )

vel = []
nline = natom/2
for i in xrange(nline) :
    line = input.readline()
    items = split( line )
    for item in items: vel.append( item )

tail = input.readline()

for icopy in xrange(ncopy) :
    outfile = "%s.%d" % (output,icopy+1)
    print "outfile:", outfile
    out = open( outfile, "w" )
    out.write( title )
    out.write( "%6d\n" % natomCL )

    for iatomCL in xrange(natomCL) :
        iatom = iatomCL * ncopy + icopy
        out.write( "%12s" % pos[ 3*iatom ] )
        out.write( "%12s" % pos[ 3*iatom+1 ] )
        out.write( "%12s" % pos[ 3*iatom+2 ] )
        if( iatomCL%2 == 1 ) : out.write( "\n" )

    for iatomCL in xrange(natomCL) :
        iatom = iatomCL * ncopy + icopy
        out.write( "%12s" % vel[ 3*iatom ] )
        out.write( "%12s" % vel[ 3*iatom+1 ] )
        out.write( "%12s" % vel[ 3*iatom+2 ] )
        if( iatomCL%2 == 1 ) : out.write( "\n" )
    out.write( tail )

