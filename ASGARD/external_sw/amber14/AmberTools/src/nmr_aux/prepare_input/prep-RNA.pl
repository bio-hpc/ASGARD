#!/usr/bin/env perl -w

#===========================================
#  USER DEFINED VARIABLES

$cyana_dir="../../../cyana/pbs/batch1";   #location of the CYANA pdb files
$number_of_strands=1;                     #number of RNA strands (1 or 2)
$start_strand1=125;                       #start of first strand (with 5'-phosphate)

$end_strand1=0;                           #end of first strand (=0 if molecule contains 1 strand)
$start_strand2=0;                         #start of second strand (=0 if molecule contains 1 strand)

#============================================

@pdbfiles = glob("${cyana_dir}/*.pdb");    #get all the CYANA pdb files
$tmpoutfile="for-leap.pdb";                #temp pdb file for tleap 

foreach $file (@pdbfiles) {
    @fileroot = split(/\//,$file);   #  split filename by '\'
    $n = @fileroot;                  #  get the index of the last element of the fileroot array
    $n--;			     #  index starts from 0
    $filename = $fileroot[$n];	     #  get the filename only without dir path
    $filename =~ /(.*)\.(.*)/;       #  store the two parts of the filename seperated by the '.' char
    $filename_no_ext=$1;             #  get the filename without the suffix
   
### Open cyana pdb files for read, and temp file "for-leap.pdb" for write 
    open(IN, $file) or die "Unable to open \'$file\'\n";
    open(OUT, "> $tmpoutfile") or die "Unable to open \'$tmpoutfile\' for writing \n";

    while($line = <IN>) {
         next unless $line =~ /^ATOM|^HETATM/;   # read in and process each line starts with the words
         ( $atno, $atname, $alt, $resname, $chain, $resno, $x, $y, $z ) = unpack("x7 a4 x a4 a a4 a1 a4 x4 a8 a8 a8",$line);

         #
         #  do modifications necessary here:
         #
         
         #   Generic substitutions, should work for all RNA's:
         
         $alt = " ";
         
         $atname = "HO2'" if $atname eq " HO2";
         $atname = "H5''" if $atname eq " H5\"";
         
         $resname = "  G" if $resname eq "RGUA";
         $resname = "  C" if $resname eq "RCYT";
         $resname = "  A" if $resname eq "RADE";
         $resname = "  U" if $resname eq "URA ";
         
         #=========================================================================
         #   System-specific changes: edit for your particular problem:
         #
         #       (In this example, we convert the cyana residue numbers to Amber
         #       residue numbers, skip the cyana linker residues, and remove 
         #       the 5'-terminal phosphates from each of the Amber strands.)

   if ($number_of_strands == 1) { 
         next if ($resno ==  ${start_strand1} || $resno == ${start_strand2}) && ($atname eq " P  " || $atname eq " OP1" || $atname eq " OP2");
   }

   elsif ($number_of_strands == 2) {
         $start_skip1=$end_strand1+1;
         $end_skip1=$start_strand2-1;
         $chain = "B" if $resno > ${end_strand1};  #end of the first RNA strand
         next if $resno >= ${start_skip1} && $resno <= ${end_skip1};  # skip cyana linker residues
         next if ($resno ==  ${start_strand1} || $resno == ${start_strand2}) && ($atname eq " P  " || $atname eq " OP1" || $atname eq " OP2");
   }
         
printf OUT "ATOM   %4s %-4s %3s %1s%4d    %8.3f%8.3f%8.3f\n", $atno,$atname,$resname,$chain, $resno,$x,$y,$z;
         }
	 close(IN); 
         close(OUT);

#
#   generate the leap.in file for tleap, similar to the cat <<EOF in csh script
#   There CAN NOT BE any space in the beginning of the lines for the "print <<EOF .... EOF" part
#
         open(OUT, ">leap.in") or die "Unable to open file leap.in for writing\n";
print OUT <<EOF;
source leaprc.ff12SB
set default nocenter on
x = loadpdb for-leap.pdb
saveamberparm x amber.top amber.rst
savepdb x amber.pdb
quit
EOF
         close(OUT);
   
    	 system("tleap -f leap.in > tleap.errs");
	 system("mv amber.pdb ${filename}");
         system("mv amber.rst ${filename_no_ext}.rst");
}

