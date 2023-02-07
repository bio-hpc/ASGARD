# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'
# vim: set filetype=perl :

######################### We start with some black magic to print on failure.
BEGIN { $| = 1; }
END {print "1..1\nnot ok 1   Premature end of Script\n" unless $loaded;}
######################### End of black magic.

$testnum = 1;

eval { require Math::VectorReal; };
if ( $@ ) {
  # The module is NOT available so fake a successfull test
  print STDERR "  Failed to locate module  Math::VectorReal  under test.\n";
  print STDERR "  Aborting testing process...\n";
  exit 0
}

# Matrix tests are useless if  Math::MatrixReal has not been installed
# However it is NOT an error if that module is not installed so abort
# if that is the case.
eval { require Math::MatrixReal; };
if ( $@ ) {
  # The module is NOT available so fake a successfull test
  print STDERR "  WARNING:  Module  Math::MatrixReal  is not available.\n";
  print STDERR "  As this module can use it but does not require it, I will\n";
  print STDERR "  skip the testing of the Math::MatrixReal interface.....\n";
  print "1..16\n";
  $skip_matrix_tests = 1;
} else {
  print "1..26\n";
}
$loaded = 1;
print "ok ", $testnum++, "\tModule Math::VectorReal Located\n";

# ---------------------------------------
# Compare a script and its previous output

sub check_script_output {
 my $script = shift;

  print "not " unless -f $script.'.pl';  # script found?
  print "ok ", $testnum++, "\t--- $script script ----\n";

  open(OUT, "$script.out") || die;
  open(TEST, "- |") or exec('perl', $script.'.pl') or exit 0;

  $/='';  # read sections by paragraph
  while( $t = <TEST> ) {
    $o = <OUT>;
    ($testname) = split(/\n/,$t);

    # ignore any white space funny busness
    $t =~ s/\s+/ /; $t =~ s/^\s//; $t =~ s/\s$//;
    $o =~ s/\s+/ /; $o =~ s/^\s//; $o =~ s/\s$//;
    print "not "  if $t ne $o;
    print "ok ", $testnum++, "\t$testname\n"
  }
  close(OUT);
  close(TEST);

  print "not "  if $?;
  print "ok ", $testnum++, "\t--- exit of $script script ----\n";

}

# ======================================

&check_script_output( "vector_test" );  # tests 2..16

exit 0 if $skip_matrix_tests;   # proceed with matrix interface tests?
print "ok ", $testnum++, "\tModule Math::MatrixReal Located\n";

&check_script_output( "matrix_test" );  # tests 18..26

&check_script_output( "synopsis" );     # tests 27..39

# ======================================
