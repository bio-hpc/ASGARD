#!/usr/bin/perl
BEGIN {                    # Add executable dir to INC to find CPin.pm
  $0 =~ /^(.*?)\/([^\/]+)$/;	# $1 Path $2 filename
  push @INC, $1;
}
use strict;
use IO::File;
use CPin;

main();

sub main {
  if (@ARGV == 0 or not (-f $ARGV[0])) {
    usageMessage();
    exit;
  }
  my $cpin = new CPinNamelist(new IO::File($ARGV[0]));
  shift @ARGV;
  $ARGV[0] = '-' unless defined $ARGV[0]; # Input from STDIN or filename args


  
  my ($pH,$firstRecord) = (7.0,1);
  my (@States, @StateCounts, @AttemptTrans, @ActualTrans);
  my $protCounts = getProtCounts($cpin);
  foreach my $i (0..@$protCounts-1) {
    $StateCounts[$i] = [ map {0} 0..@{$protCounts->[$i]}-1 ];
    $AttemptTrans[$i] = $ActualTrans[$i] = 0;
  }

  while (@ARGV) {
    
    my $filename = shift @ARGV;
    (-r $filename) or $filename eq '-' or die "Can't open $filename for reading.\n";
    my $sfh = new IO::File($filename);
  
    while (<$sfh>) {            # Loop to parse cpout files
      if (/^Residue\s+(\d+)\s+State:\s+(\d+)/) {
	my ($residue, $state) = ($1, $2);
        if (not $firstRecord) {
          $AttemptTrans[$residue]++; # Minor bug: this increments every residue on a "full" record
          $ActualTrans[$residue]++
            if abs($protCounts->[$residue][$state] - $protCounts->[$residue][$States[$residue]]) != 0;
        }
	$States[$residue] = $state;
      } elsif (/^Solvent pH:\s+(\S+)/) {
        $pH = $1;
      } elsif ("\n" eq $_) {    # Newline means end of record for this MC step
        if ($firstRecord) {
          $firstRecord = 0;
        } else {
          foreach my $i (0..@States-1) {
            $StateCounts[$i][$States[$i]]++;
          }
        }
      } else {
        # Ignore unneeded lines (e.g. Time, Time step, etc.)
      }
    }
  }
  
  my $totprot = 0;
  for (my $resIndex = 0; $resIndex <@States; $resIndex++) {
    my $res = $cpin->getRes($resIndex);
    printf "%3s %4s: ", $res->ResName(), $res->ResNum();
    my $resArr = $StateCounts[$resIndex];
    # Final arg of calcOffset: 1 -> Allow +/-"Inf" predictions; 0 -> Only numeric predictions
    my $offset = calcOffset($resArr, $protCounts->[$resIndex], 1); 
    my $frac = calcFraction($resArr, $protCounts->[$resIndex]);
    $totprot += $frac;
    my ($predict, $exp, $err);
    if (substr($offset, -1, 1) eq 'f') { # -Inf or Inf string
      $predict = $offset;
      $err = $offset;
      $exp = sprintf("%2.3f",$res->pKa()) if defined $res->pKa();
    } else {
      $predict = sprintf("% 2.3f",$pH + $offset);
      if (defined $res->pKa()) {
        $exp = sprintf("%2.3f",$res->pKa());
        $err = sprintf("% 2.3f",$predict - $exp);
      }
      $offset = sprintf("% 2.3f",$offset);
    }
    print "Offset $offset  Pred $predict";
    print "  Exp $exp  Err $err" if defined $exp;
    printf("  Frac Prot %2.3f  Transitions %9d",$frac, $ActualTrans[$resIndex]);
    print "\n";
#      for (my $state = 0; $state < @{$StateCounts[$resIndex]}; $state++) {
#        print "      State $state:",($StateCounts[$resIndex][$state] or 0),"\n";
#      }
  }
  printf "\n\nAverage total molecular protonation: %3.3f\n", $totprot;

}

sub calcOffset {
  my ($resArr, $relProt, $infOK) = @_;
  my ($prot,$deprot) = getProtPopulations($resArr, $relProt);
  if ($infOK) {
    if ($deprot == 0) {
      return "Inf";
    } elsif ($prot == 0) {
      return "-Inf";
    }
  } else {
    $deprot += 0.00001;
    $prot += 0.00001;
  }
  return log($deprot/$prot)/(-log(10));
}

sub calcFraction {
  my ($resArr,$relProt) = @_;
  my ($prot, $deprot) = getProtPopulations($resArr, $relProt);
  return $prot/($deprot+$prot);
}

sub getProtPopulations {
  my ($resArr, $relProt) = @_;
  my @sortedProt = sort {$a <=> $b} @$relProt;
  my ($min, $max) = ($sortedProt[0],$sortedProt[$#sortedProt]);
  $max - $min == 1 or warn "Can't handle more than 2 protonation levels; results are probably wrong\n";
  my $prot = sum(map {$relProt->[$_] == $max ? $resArr->[$_] : 0} 0..$#$resArr);
  my $deprot =  sum(map {$relProt->[$_] == $min ? $resArr->[$_] : 0} 0..$#$resArr);
  return ($prot,$deprot);
}

sub sum {
  my $total=0;
  foreach my $element (@_) {
    $total+=$element if defined $element;
  }
  return $total;
}

# Return array of residude number by state number, where value is relative protonation
sub getProtCounts {
  my ($cpin) = @_;
  my $protCounts = [];
  foreach my $resIndex (0..$cpin->ResCount()-1) {
    my $stateData = $CPinData::CHRG{$cpin->getRes($resIndex)->ResName()};
    $protCounts->[$resIndex] = [ map {$stateData->[$_][1]} 0..@$stateData-1 ];
  }
  return $protCounts;
}
sub usageMessage {
  print <<USEMSG;
CalcpKa -- Analyze constant pH SANDER output to calculate pKa values

caclpka.pl <cpinfile> <cpoutfile> [ <cpoutfile> [ <cpoutfile> ... ] ]

USEMSG
  return;
}
