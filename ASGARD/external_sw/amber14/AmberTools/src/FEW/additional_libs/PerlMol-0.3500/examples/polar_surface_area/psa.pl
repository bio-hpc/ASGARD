#!/usr/bin/perl
#
# http://www.daylight.com/support/contrib/ertl/tpsa.pl
# This script is a conversion of Peter Ertl's c program tpsa.c into Perl.
# I added a few additional features, like the ability to use files as input
# and output.  -- Eric Allen
# 
# +++++++++++++++++++++++++
# 
#  Calculation of polar surface area based on fragment contributions (TPSA)
#  Peter Ertl, Novartis Pharma AG, Basel, August 2000
#  peter.ertl@pharma.novartis.com || peter.ertl@altavista.net
#  The program uses the SMILES Daylight Toolkit module.
#
#  For description of the methodology see :
#  P. Ertl, B. Rohde, P. Selzer
#  Fast Calculation of Molecular Polar Surface Area as a Sum of Fragment-based
#  Contributions and Its Application to the Prediction of Drug Transport 
#  Properties, J.Med.Chem. 43, 3714-3717, 2000
#
#  The program reads a list of SMILES strings from stdin, and writes SMILES
#  and the calculated TPSA on stdout.
#  The program expects "problematic" nitrogens in pentavalent form, i.e. nitro 
#  groups as -N(=O)=O, N-oxides as c1ccccn1=O, and azides as -N=N#N.  
#
#  This contrib program comes with absolutely no warranty of any kind.
#
# +++++++++++++++++++++++++
#
# Revision History
#
# Date         Who                 What
# -----------  ------------------  --------------------------------------------
# 08-MAR-2001  Eric Allen          Converted c program to Perl.
# 13-MAY-2005  Konrad Koehler      Replaced Daylight Dependencies with PerlMol modules.
# 

use strict;

# KFK: Commented out the following line and replaced by the following 9 lines:
#use DayPerl;
use Chemistry::Mol;
use Chemistry::Atom;
use Chemistry::Pattern;
use Chemistry::File::SMILES;
use Chemistry::File::SMARTS;
use Chemistry::Bond::Find;
use Chemistry::Ring;
use Chemistry::Ring::Find ':all';
use Chemistry::Ring 'aromatize_mol';


my %inputArguments = @ARGV;
my $calculatedOne = 0;
my $usageString = "usage:\n\tcalculatepsa.pl 'smiles string'\nor\n\tcalculatepsa.pl -f filename [-o outputfilename]\n\n";

my $line;
my $smiles;
my $id;
my $psa;
my $message;
my $filename;
my $outfilename;

my $writeToFile = 0;

if($#ARGV > 0)
{
	# means we have multiple arguements to deal with
	
	if(defined($inputArguments{-f}))
	{
		open INFILE,"<$inputArguments{-f}" or die "Could not open input file $inputArguments{-f}.\n";
	}
	else
	{
		die "$usageString";
	}
	
	if(defined($inputArguments{-o}))
	{
		open OUTFILE,">$inputArguments{-o}" or die "Could not open output file $inputArguments{-o}.\n";
		$writeToFile = 1;
	}
	
	while($smiles = <INFILE>)
	{
		chomp $smiles;
		$message = '';
		($psa,$message) = &CalculatePSA($smiles);
		if($writeToFile)
		{
			if(length($message) > 0)
			{
				print OUTFILE "$smiles $message\n";
			}
			else
			{
				print OUTFILE "$smiles\t$psa\n";
			}
		}
		else
		{
			if(length($message) > 0)
			{
				print "$smiles $message\n";
			}
			else
			{
				print "$smiles\t$psa\n";
			}
		}
	}
	close INFILE;
	if($writeToFile)
	{
		close OUTFILE;
	}	
}
elsif($#ARGV == 0)
{
	# means we have a single smiles string to use
	
	$smiles = $ARGV[0];
	chomp $smiles;
	($psa,$message) = &CalculatePSA($smiles);
	if(length($message) > 0)
	{
		print "$smiles $message\n";
	}
	else
	{
		print "$smiles\t$psa\n";
	}
}
else
{
	$smiles = <STDIN>;
	chomp $smiles;
	while(length($smiles) > 0)
	{
		$calculatedOne = 1;
		($psa,$message) = &CalculatePSA($smiles);
		if(length($message) > 0)
		{
			print "$smiles $message\n";
		}
		else
		{
			print "$smiles\t$psa\n";
		}
		$smiles = <STDIN>;
		chomp $smiles;
	}
	unless($calculatedOne)
	{
		print "$usageString";
	}
}


#################
                 # CalculatePSA
#################

sub CalculatePSA
{
	my ($smiles) = @_;
		
	my $atomic_number_nitrogen = 7;
	my $atomic_number_oxygen = 8;
	my $max_ring_count = 3;
	
	my $psa = 0;
	my $message;
	my $molecule;
	my $atoms;
	my $rings;
	my $bonds;
	my $oneAtom;
	my $oneBond;
	my $atomicNumber;
	my $isIn3Ring;
	my $hydrogenCount;
	my $formalCharge;
	my $isAromatic;
	my $oneRing;
	my $ringAtoms;
	my $oneRingAtom;
	my $resetRings;
	my $numberOfNeighbors;
	my $singleBondCount;
	my $doubleBondCount;
	my $tripleBondCount;
	my $aromaticBondCount;
	my $bondType;
	my $polarity;
	my $deAllocated;
	
	my $NTriple = 23.79;
	my $NHDouble = 23.85;
	my $NH2Single = 26.02;
	my $NH2PlusDouble = 25.59;
	my $NH3PlusSingle = 27.64;
	my $NSingleDouble = 12.36;
	my $NDoubleTriple = 13.60;
	my $NHSingleSingle = 12.03;
	my $NHSingleSingle3Ring = 21.94;
	my $NPlusSingleTriple = 4.36;
	my $NHPlusSingleDouble = 13.97;
	my $NH2PlusSingleSingle = 16.61;
	my $nAromaticAromatic = 12.89;
	my $nHAromaticAromatic = 15.79;
	my $nHPlusAromaticAromatic = 14.14;
	my $NSingleSingleSingle = 3.24;
	my $NSingleSingleSingle3Ring = 3.01;
	my $NSingleDoubleDouble = 11.68;
	my $NPlusSingleSingleDouble = 3.01;
	my $NPlusSingleSingleSingle = 4.44;
	my $nAromaticAromaticAromatic = 4.41;
	my $nSingleAromaticAromatic = 4.93;
	my $nDoubleAromaticAromatic = 8.39;
	my $nPlusAromaticAromaticAromatic = 4.10;
	my $nPlusSingleAromaticAromatic = 3.88;
	my $NPlusSingleSingleSingleSingle = 0.00;
	my $nonStandardNBasePolarity = 30.5;
	my $neighborNBasePolarity = 8.2;
	my $hydrogenBasePolarity = 1.5;
	my $ODouble = 17.07;
	my $OHSingle = 20.23;
	my $OMinusSingle = 23.06;
	my $OSingleSingle = 9.23;
	my $OSingleSingle3Ring = 12.53;
	my $oAromaticAromatic = 13.14;
	my $nonStandardOBasePolarity = 28.5;
	my $neighborOBasePolarity = 8.6;
	
	
    # KFK: replaced the following line with the next two lines:
	# $molecule = dt_smilin($smiles);
    $molecule = Chemistry::Mol->parse($smiles, format => 'smiles');
    aromatize_mol($molecule); # <--- AROMATIZE!!!

	unless(defined($molecule))
	{
		$message = "SMILES string invalid";
		return ($psa,$message);
	}

    # KFK: replaced the next two lines with the following line:
	#$atoms = dt_stream($molecule,TYP_ATOM);
	#$rings = dt_stream($molecule,TYP_CYCLE);
	

    # KFK: added the next block of code to recplace the following block:
    my @atoms = $molecule->atoms;

    foreach (@atoms) {

      $oneAtom       = $_;
      $oneAtom->collapse_hydrogens;
      $atomicNumber  = $oneAtom->Z;
      
      if(($atomicNumber == $atomic_number_nitrogen) || ($atomicNumber == $atomic_number_oxygen)) {
      
      $hydrogenCount = $oneAtom->implicit_hydrogens;
      $formalCharge  = $oneAtom->formal_charge;
      if ($formalCharge == "") { $formalCharge = 0; }
      $isAromatic    = $oneAtom->aromatic;
      if (find_ring($oneAtom, size => '3') == "") {
        $isIn3Ring = 0;
      } else {
        $isIn3Ring = 1;
      }
      $numberOfNeighbors = $oneAtom->neighbors($oneAtom);
      $singleBondCount   = 0;
      $doubleBondCount   = 0;
      $tripleBondCount   = 0;
      $aromaticBondCount = 0;
      my @bonds = $oneAtom->bonds($oneAtom);
      foreach (@bonds) {
        $oneBond = $_;
        my $bo = $oneBond->order;
        if ($bo == 1) {  $singleBondCount++; }
        if ($bo == 2) {  $doubleBondCount++; }
        if ($bo == 3) {  $tripleBondCount++; }
        if ($oneBond->aromatic) {  $aromaticBondCount++; }
      }
	
#	while($oneAtom = dt_next($atoms))
#	{
#		$atomicNumber = dt_number($oneAtom);
#		if(($atomicNumber == $atomic_number_nitrogen) || ($atomicNumber == $atomic_number_oxygen))
#		{
#			$hydrogenCount = dt_hcount($oneAtom);
#			$formalCharge = dt_charge($oneAtom);
#			$isAromatic = dt_aromatic($oneAtom);
#			$isIn3Ring = 0;
#			RINGCHECK:
#			while($oneRing = dt_next($rings))
#			{
#				if(dt_count($oneRing,TYP_ATOM) > $max_ring_count)
#				{
#					next RINGCHECK;
#				}
#				$ringAtoms = dt_stream($oneRing,TYP_ATOM);
#				while($oneRingAtom = dt_next($ringAtoms))
#				{
#					if($oneAtom == $oneRingAtom)
#					{
#						$isIn3Ring = 1;
#						last RINGCHECK;
#					}
#				}
#			}
#			$resetRings = dt_reset($rings);
#			
#			$numberOfNeighbors = 0;
#			$singleBondCount = 0;
#			$doubleBondCount = 0;
#			$tripleBondCount = 0;
#			$aromaticBondCount = 0;
#			$bonds = dt_stream($oneAtom,TYP_BOND);
#			while($oneBond = dt_next($bonds))
#			{
#				$numberOfNeighbors++;
#				$bondType = dt_bondtype($oneBond);
#				if($bondType == DX_BTY_SINGLE)
#				{
#					$singleBondCount++;
#				}
#				elsif($bondType == DX_BTY_DOUBLE)
#				{
#					$doubleBondCount++;
#				}
#				elsif($bondType == DX_BTY_TRIPLE)
#				{
#					$tripleBondCount++;
#				}
#				elsif($bondType == DX_BTY_AROMAT)
#				{
#					$aromaticBondCount++;
#				}
#			}
#			
			$polarity = -1;
			if($atomicNumber == $atomic_number_nitrogen)
			{
				if($numberOfNeighbors == 1)
				{
					if(($hydrogenCount == 0) && ($tripleBondCount == 1) && ($formalCharge == 0))
					{
						# N#
						$polarity = $NTriple;
					}
					elsif(($hydrogenCount == 1) && ($doubleBondCount == 1) && ($formalCharge == 0))
					{
						# [NH]=
						$polarity = $NHDouble;
					}
					elsif(($hydrogenCount == 2) && ($singleBondCount == 1) && ($formalCharge == 0))
					{
						# [NH2]-
						$polarity = $NH2Single;
					}
					elsif(($hydrogenCount == 2) && ($doubleBondCount == 1) && ($formalCharge == 1))
					{
						# [NH2+]=
						$polarity = $NH2PlusDouble;
					}
					elsif(($hydrogenCount == 3) && ($singleBondCount == 1) && ($formalCharge == 1))
					{
						# [NH3+]-
						$polarity = $NH3PlusSingle;
					}
				}
				elsif($numberOfNeighbors == 2)
				{
					if(($hydrogenCount == 0) && ($singleBondCount == 1) && ($doubleBondCount == 1) && ($formalCharge == 0))
					{
						$polarity = $NSingleDouble;
					}
					elsif(($hydrogenCount == 0) && ($tripleBondCount == 1) && ($doubleBondCount == 1) && ($formalCharge == 0))
					{
						$polarity = $NDoubleTriple;
					}
					elsif(($hydrogenCount == 1) && ($singleBondCount == 2) && ($formalCharge == 0) && !($isIn3Ring))
					{
						$polarity = $NHSingleSingle;
					}
					elsif(($hydrogenCount == 1) && ($singleBondCount == 2) && ($formalCharge == 0) && ($isIn3Ring))
					{
						$polarity = $NHSingleSingle3Ring;
					}
					elsif(($hydrogenCount == 0) && ($tripleBondCount == 1) && ($singleBondCount == 1) && ($formalCharge == 1))
					{
						$polarity = $NPlusSingleTriple;
					}
					elsif(($hydrogenCount == 1) && ($doubleBondCount == 1) && ($singleBondCount == 1) && ($formalCharge == 1))
					{
						$polarity = $NHPlusSingleDouble;
					}
					elsif(($hydrogenCount == 2) && ($singleBondCount == 2) && ($formalCharge == 1))
					{
						$polarity = $NH2PlusSingleSingle;
					}
					elsif(($hydrogenCount == 0) && ($aromaticBondCount == 2) && ($formalCharge == 0))
					{
						$polarity = $nAromaticAromatic;
					}
					elsif(($hydrogenCount == 1) && ($aromaticBondCount == 2) && ($formalCharge == 0))
					{
						$polarity = $nHAromaticAromatic;
					}
					elsif(($hydrogenCount == 1) && ($aromaticBondCount == 2) && ($formalCharge == 1))
					{
						$polarity = $nHPlusAromaticAromatic;
					}
				}
				elsif($numberOfNeighbors == 3)
				{
					if(($hydrogenCount == 0) && ($singleBondCount == 3) && ($formalCharge == 0) && !($isIn3Ring))
					{
						$polarity = $NSingleSingleSingle;
					}
					elsif(($hydrogenCount == 0) && ($singleBondCount == 3) && ($formalCharge == 0) && ($isIn3Ring))
					{
						$polarity = $NSingleSingleSingle3Ring;
					}
					elsif(($hydrogenCount == 0) && ($singleBondCount == 1) && ($doubleBondCount == 2) && ($formalCharge == 0))
					{
						$polarity = $NSingleDoubleDouble;
					}
					elsif(($hydrogenCount == 0) && ($singleBondCount == 2) && ($doubleBondCount == 1) && ($formalCharge == 1))
					{
						$polarity = $NPlusSingleSingleDouble;
					}
					elsif(($hydrogenCount == 1) && ($singleBondCount == 3) && ($formalCharge == 1))
					{
						$polarity = $NPlusSingleSingleSingle;
					}
					elsif(($hydrogenCount == 0) && ($aromaticBondCount == 3) && ($formalCharge == 0))
					{
						$polarity = $nAromaticAromaticAromatic;
					}
					elsif(($hydrogenCount == 0) && ($singleBondCount == 1) && ($aromaticBondCount == 2) && ($formalCharge == 0))
					{
						$polarity = $nSingleAromaticAromatic;
					}
					elsif(($hydrogenCount == 0) && ($aromaticBondCount == 3) && ($formalCharge == 1))
					{
						$polarity = $nPlusAromaticAromaticAromatic;
					}
					elsif(($hydrogenCount == 0) && ($singleBondCount == 1) && ($aromaticBondCount == 2) && ($formalCharge == 1))
					{
						$polarity = $nPlusSingleAromaticAromatic;
					}
				}
				elsif($numberOfNeighbors == 4)
				{
					if(($hydrogenCount == 0) && ($singleBondCount == 4) && ($formalCharge == 1))
					{
						$polarity = $NPlusSingleSingleSingleSingle;
					}
				}
				
				if($polarity < 0) # for an N with non-standard valence
				{
					$polarity = $nonStandardNBasePolarity - ($numberOfNeighbors * $neighborNBasePolarity) + ($hydrogenCount * $hydrogenBasePolarity);
				}
				
				if($polarity < 0)
				{
					$polarity = 0;
				}
			}
			elsif($atomicNumber == $atomic_number_oxygen)
			{
				if($numberOfNeighbors == 1)
				{
					if(($hydrogenCount == 0) && ($doubleBondCount == 1) && ($formalCharge == 0))
					{
						$polarity = $ODouble;
					}
					elsif(($hydrogenCount == 1) && ($singleBondCount == 1) && ($formalCharge == 0))
					{
						$polarity = $OHSingle;
					}
					elsif(($hydrogenCount == 0) && ($singleBondCount == 1) && ($formalCharge == -1))
					{
						$polarity = $OMinusSingle;
					}
				}
				elsif($numberOfNeighbors == 2)
				{
					if(($hydrogenCount == 0) && ($singleBondCount == 2) && ($formalCharge == 0) && !($isIn3Ring))
					{
						$polarity = $OSingleSingle;
					}
					elsif(($hydrogenCount == 0) && ($singleBondCount == 2) && ($formalCharge == 0) && ($isIn3Ring))
					{
						$polarity = $OSingleSingle3Ring;
					}
					elsif(($hydrogenCount == 0) && ($aromaticBondCount == 2) && ($formalCharge == 0))
					{
						$polarity = $oAromaticAromatic;
					}
				}
				
				if($polarity < 0) # O with non standard polarity
				{
					$polarity = $nonStandardOBasePolarity - ($numberOfNeighbors * $neighborOBasePolarity) + ($hydrogenCount * $hydrogenBasePolarity);
				}
				
				if($polarity < 0)
				{
					$polarity = 0;
				}
			}
			
			$psa += $polarity;
		}
	}
	

# KFK: commented out the following deallocation calls which are no-longer needed:
#	$deAllocated = dt_dealloc($atoms);
#	$deAllocated = dt_dealloc($bonds);
#	$deAllocated = dt_dealloc($rings);
#	$deAllocated = dt_dealloc($oneAtom);
#	$deAllocated = dt_dealloc($oneBond);
#	$deAllocated = dt_dealloc($oneRing);
#	$deAllocated = dt_dealloc($ringAtoms);
#	$deAllocated = dt_dealloc($oneRingAtom);
#	$deAllocated = dt_dealloc($molecule);

	return $psa,$message;
}

#################
                 #
#################


