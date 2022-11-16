# Package for reading/writing/manipulating AMBER constant pH
# cpin files.
# Also includes a Namelist parser and data, to reduce file count
# Programmed by John Mongan <jmongan@mccammon.ucsd.edu>, 2003
#
package CPin;
use IO::File;
use Data::Dumper;
use Carp;

sub new {
  my $type = shift;
  my $this = {ResCount => 0, StepPeriod => undef, pH => undef, SysName => undef};
  bless $this, $type;

  return $this;
}

sub output {
  my $this = shift;
  my $out = "";
  $this->{ResCount} = $this->ResCount();
  $out .= join("\n",map {$this->{$_}} qw(ResCount StepPeriod pH))."\n";
  foreach my $res (@{$this->{Residue}}) {
    $res->{AtomCount} = $res->{State}[0]->AtomCount();
    $res->{StateCount} = $res->StateCount();
    $out .= join("\n", map {$res->{$_}} qw(StateCount InitState FirstAtom AtomCount))."\n";
    foreach my $state (@{$res->{State}}) {
      $out .= join("\n", map {$state->{$_}} qw(Energy Protonation))."\n";
      $out .= join("\n", @{$state->{Charge}})."\n";
    }
  }
  $out .= "# System: ".$this->SysName()."\n";
  for (my $i = 0; $i < $this->ResCount(); $i++) {
    $out .= "# Residue: ".$this->getRes($i)->ResName()." ".$this->getRes($i)->ResNum();
    $out .= " ".$this->getRes($i)->pKa() if defined($this->getRes($i)->pKa());
    $out .= "\n";
  }
  return $out;
}

sub outputNml {
  my $this = shift;
  my $nml = Namelist->new();
  $nml->Name("CNSTPH");
  my $ndata = $nml->{data};

#   my @RES_FIELDS = qw(StateCount InitState FirstAtom AtomCount);
#   my @NML_RES_FIELDS = qw(NUM_STATES FIRST_STATE FIRST_ATOM NUM_ATOMS);

  my %ResMap;
  
  my ($atomCnt,$stateCnt) = (0,0);
  
  $ndata->{TRESCNT} = $this->{ResCount} = $this->ResCount();
  $ndata->{RESNAME}[0] = "'System: ".$this->SysName()."'" if defined $this->SysName();
  for (my $itres = 0; $itres < $this->{ResCount}; $itres++) {
    my $res = $this->{Residue}[$itres];
    $ndata->{STATEINF}[$itres]{NUM_STATES} = $res->StateCount();
    $ndata->{STATEINF}[$itres]{FIRST_ATOM} = $res->{FirstAtom};
    $ndata->{STATEINF}[$itres]{NUM_ATOMS} = $res->{State}[0]->AtomCount();
    $ndata->{RESSTATE}[$itres] = $res->{InitState};
    my $resname = join("", "'Residue: ", $res->ResName()," ",$res->ResNum());
    $resname = join("", $resname, " ", $res->pKa()) if defined $res->pKa();
    $ndata->{RESNAME}[$itres+1] = $resname."'";
    if (defined $res->ResName() and exists $ResMap{ $res->ResName() }) {
      ($ndata->{STATEINF}[$itres]{FIRST_STATE}, $ndata->{STATEINF}[$itres]{FIRST_CHARGE}) =
        @{$ResMap{ $res->ResName()}};
    } else {
      $ndata->{STATEINF}[$itres]{FIRST_STATE} = $stateCnt;
      $ndata->{STATEINF}[$itres]{FIRST_CHARGE} = $atomCnt;
      $ResMap{ $res->ResName() } = [$stateCnt, $atomCnt] if defined $res->ResName();
      foreach my $state (@{$res->{State}}) {
        $ndata->{PROTCNT}[$stateCnt] = $state->{Protonation};
        $ndata->{STATENE}[$stateCnt++] = $state->{Energy};
        foreach my $chrg (@{$state->{Charge}}) {
          $ndata->{CHRGDAT}[$atomCnt++] = $chrg;
        }
      }
    }
  }
  return $nml->namelistString();
}

sub setStatesFromOutputFile {
  my $this = shift;
  my ($filename) = @_;
  my $fh = new IO::File($filename) or die "Can't open $filename";
  $fh->seek(-4096,2);		# 4096 bytes from file end
  my @state;
  while (<$fh>) {
    if (/TRES\s+(\d+)\s+STATE\s+(\d+)/) {
      my ($residue, $state) = ($1, $2);
      $state[$residue] = $state;
    }
  }
  $this->setStates(@state);
}

sub setStates {
  my $this = shift;
  my @state = @_;
  for (my $i = 0; $i < @{$this->{Residue}}; $i++) {
    $this->{Residue}[$i]{InitState} = $state[$i] if defined $state[$i];
  }
}

sub StepPeriod {
  my $this = shift;
  $this->{StepPeriod} = $_[0] if defined $_[0];
  carp "Deprecated -- period is now an mdin parameter";
  return $this->{StepPeriod};
}
sub pH {
  my $this = shift;
  $this->{pH} = $_[0] if defined $_[0];
  carp "Deprecated -- pH is now an mdin parameter";
  return $this->{pH};
}
sub SysName {
  my $this = shift;
  $this->{SysName} = $_[0] if defined $_[0];
  return $this->{SysName};
}
sub getRes {
  my $this = shift;
  return $this->{Residue}[$_[0]];
}
sub ResCount {
  my $this = shift;
  return scalar(@{$this->{Residue}});
}
sub AddResidue {
  my $this = shift;
  my $newRes = Residue->new(@_);
  defined($newRes) or die "Couldn't create new residue";
  push @{$this->{Residue}}, $newRes;
}
sub DeleteResidueNum {
  my $this = shift;
  my ($resNum) = @_;
  splice @{$this->{Residue}}, $resNum, 1;
}
package Residue;
sub new {
  my $type = shift;
  my ($resName, $resNum, $firstAtom) = @_;
  @_ == 3 or die "Need resName, Number and First Atom to construct residue";
  my $this = {StateCount => 0, InitState => 0, FirstAtom =>$firstAtom,
	      ResName => $resName, ResNum => $resNum, pKa => undef};
  bless $this, $type;
  my $resData = $CPinData::CHRG{$resName};
  $this->FirstAtom($firstAtom);
  foreach my $state (@$resData) {
    push @{$this->{State}}, new State(@$state);
  }
  return $this;
}
sub InitState {
  my $this = shift;
  $this->{InitState} = $_[0] if defined $_[0];
  return $this->{InitState};
}
sub pKa {
  my $this = shift;
  $this->{pKa} = $_[0] if defined $_[0];
  return ($this->{pKa});
}
sub ResNum {
  my $this = shift;
  $this->{ResNum} = $_[0] if defined $_[0];
  return $this->{ResNum};
}
sub ResName {
  my $this = shift;
  $this->{ResName} = $_[0] if defined $_[0];
  return $this->{ResName};
}
sub FirstAtom {
  my $this = shift;
  $this->{FirstAtom} = $_[0] if defined $_[0];
  return $this->{FirstAtom};
}
sub StateCount {
  my $this = shift;
  return scalar(@{$this->{State}});
}


package State;
sub new {
  my $type = shift;
  my $this = {Energy => shift, Protonation => shift};
  $this->{Charge} = \ @_;
  bless $this, $type;
  return $this;
}

sub AtomCount {
  my $this = shift;
  return scalar(@{$this->{Charge}});
}

# Parse a standard namelist based cpin file
package CPinNamelist;
BEGIN {@CPinNamelist::ISA = qw (CPin); }

sub new {
  my $class = shift;
  my $this = $class->SUPER::new();
  my ($nmlFH) = @_;
  my $nml = new Namelist();
  $nml->readFile($nmlFH);
  $this->_loadData($nml);
  return $this;
}

sub _loadData {
  my $this = shift;
  my ($nml) = @_;
  my $data = $nml->{data};
  if ($data->{RESNAME}[0] =~ /^'System:\s+(\S+)\s*'/ ) {
    $this->SysName($1);
  } 
  for (my $i = 0; $i < $data->{TRESCNT};$i++) {
    my $residue = {};
    bless $residue, 'Residue';
    $residue->{StateCount} = $data->{STATEINF}[$i]{NUM_STATES};
    $residue->{InitState} =     # scalar if only one titrating residue
      (ref $data->{RESSTATE} ? $data->{RESSTATE}[$i] : $data->{RESSTATE});
    $residue->{FirstAtom} = $data->{STATEINF}[$i]{FIRST_ATOM};
    $residue->{AtomCount} = $data->{STATEINF}[$i]{NUM_ATOMS};
    my $firstState = $data->{STATEINF}[$i]{FIRST_STATE};
    if ($data->{RESNAME}[$i+1] =~ /^'Residue:\s+(\S+)\s+(\d+)\s*([0-9.]+)?/) {
      $residue->ResName($1); $residue->ResNum($2);
      $residue->pKa($3) if defined $3;
    }
    for (my $j = 0; $j < $residue->{StateCount}; $j++) {
      my $state = {};
      bless $state, 'State';
      $state->{Energy} = $data->{STATENE}[$firstState+$j];
      $state->{Protonation} = $data->{PROTCNT}[$firstState+$j];
      my $chrgBegin = $data->{STATEINF}[$i]{FIRST_CHARGE}+$j*$residue->{AtomCount};
      $state->{Charge} =[ @{$data->{CHRGDAT}}[$chrgBegin .. 
                                            $chrgBegin + $residue->{AtomCount} - 1 ] ];
      $residue->{State}[$j] = $state;
    }
    $this->{Residue}[$i] = $residue;
  }
}

# Parse an old-style (deprecated) "charges" format file to produce cpin object
package CPinFile;
BEGIN { @CPinFile::ISA = qw(CPin); }
  
sub new {
  my $class = shift;
  my $this = $class->SUPER::new();
  $this->_readCPinFile(@_);
  return $this;
}

sub _readCPinFile {
  my $this = shift;
  my ($cf) = @_;
  $this->{ResCount} = nextVal($cf);
  $this->{StepPeriod} = nextVal($cf);
  $this->{pH} = nextVal($cf);
  for (my $i = 0; $i < $this->{ResCount};$i++) {
    my $residue = {};
    bless $residue, 'Residue';
    $residue->{StateCount} = nextVal($cf);
    $residue->{InitState} = nextVal($cf);
    $residue->{FirstAtom} = nextVal($cf);
    $residue->{AtomCount} = nextVal($cf);
    for (my $j = 0; $j < $residue->{StateCount}; $j++) {
      my $state = {};
      bless $state, 'State';
      $state->{Energy} = nextVal($cf);
      $state->{Protonation} = nextVal($cf);
      for (my $k = 0; $k < $residue->{AtomCount}; $k++) {
	$state->{Charge}[$k] = nextVal($cf);
      }
      $residue->{State}[$j] = $state;
    }
    $this->{Residue}[$i] = $residue;
  }
  my $i = 0;
  while (<$cf>) {
    if (/^# System:\s+(\S+)/) {
      $this->SysName($1);
    } elsif (/^# Residue:\s+(\S+)\s+(\d+)\s*([0-9.]+)?/) {
      $this->getRes($i)->ResName($1);
      $this->getRes($i)->ResNum($2);
      $this->getRes($i)->pKa($3) if defined $3;
      $i++;
    }
  }
}

sub nextVal {
  my ($fh) = @_;
  my $line = $fh->getline();
  defined $line or die "Unexpected end of file.";
  chomp $line;
  return $line;
}


# Parse a PDB file to produce the cpin object
package CPinFromPDBFile;
BEGIN { @CPinFromPDBFile::ISA = qw(CPin); }
  
sub new {
  my $class = shift;
  my $this = $class->SUPER::new();
  $this->_readPDBFile(@_);
  return $this;
}

sub _readPDBFile {
  my $this = shift;
  my ($pdb) = @_;
  my $oldResNum = 0;
  while (<$pdb>) {
    my ($garbage, $atomNum, $atomName, $resName, $resNum) = split;
    if (defined $resName and exists $CPinData::CHRG{$resName} and $resNum != $oldResNum) {
      $oldResNum = $resNum;
      push @tResName, $resName;
      push @tResNum, $resNum;
      push @tResAtomNum, $atomNum;
    }
  }

  for (my $i = 0; $i < @tResName; $i++) {
    $this->AddResidue($tResName[$i], $tResNum[$i], $tResAtomNum[$i]);
  }
}

sub setpKas {
  my $this = shift;
  my $pkas = $CPinData::PKA{$this->SysName()};
  return unless defined $pkas;
  for (my $i = 0; $i < $this->ResCount(); $i++) {
    my $res = $this->getRes($i);
    if (defined $pkas->{$res->ResNum()}) {
      $res->pKa( $pkas->{$res->ResNum()} );
    }
  }
}

# Package to parse and write F90 style namelists. The namelist object
# represents the namelist data as a hash in its data member. Namelist arrays
# are represented with Perl arrays and namelist derived types with Perl
# hashes. All hash keys are upper cased. For instance, the namelist entry
# stateinf(4)%first_atom would be located at
# $nml->{data}{STATEINF}[4]{FIRST_ATOM} namelistString() returns the data
# structure as a fortran namelist string. parseNamelist performs the inverse
# operation, populating the namelist object from a string representation of
# an F90 namelist. readFile is a convenience method to feed parseNamelist
# from an open filehandle. Both of these method expect a single, complete
# namelist as input.
#
# Parsing namelists is somewhat more difficult than it might at first appear,
# and the parser has some limitations, particularly regarding character
# strings. Only single quoted strings are allowed. Single quotes cannot
# appear within strings (no escaping is recognized). Since at least some
# fortran implementations insert newlines into strings at column 80 of the
# namelist output, newlines are not considered significant. This means
# there's no way to represent a newline in the string if you need to.
#
# Programmed by John Mongan, Autumn 2003, UCSD Bioinformatics
#
package Namelist;
use Carp;
use strict;
use IO::File;

sub new {
  my $class = shift;
  my $this = {};
  $this->{name} = undef;
  $this->{data} = {};
  bless $this, $class;
}

# Accessor method for namelist's name
#
sub Name {
  my $this = shift;
  my ($arg) = @_;
  my $oldval = $this->{name};
  $this->{name} = $arg if defined $arg;
  return $oldval;
}
# Populate data structure from filehandle (assumed to contain only a
# single namelist)
#
sub readFile {
  my $this = shift;
  my ($fh) = @_;
  local $/;
  undef $/;
  my $nmlString = <$fh>;
  $this->parseNamelist($nmlString);
  $fh->close();
}

# populate namelist object from $str argument. $str contains a single,
# complete F90 style namelist.
sub parseNamelist {
  my $this = shift;
  my ($str) = @_;
  
  if ( $str =~ s/^\s*\&(\w+)\s*// ) { 
    $this->{name} = $1;
  } else {
    carp "Can't find namelist start delimiter, aborting";
    return;
  }
  # Tokens are delimited by equals sign, comma and or
  # whitespace. Equals sign itself is a token; other delimiters are
  # removed. Delimiters within a single quoted character string are
  # not significant and must be ignored. The following statement attempts
  # to implement these rules. 
  my @tokens =  grep((defined($_) and $_ ne ""),
                     split(/
                             (?:(?:[\s,]+|(?<==)) # Comma\space delim or lookbehind =
                                ((?:\d+\*)?       # optional repeat specifier 
                                '[^']*'))         # char string (only ' delim significant) 
                             | [\s,]+             # comma\space delimiter
                             | (=)                # = delim, with parens to capture =
                           /x, $str));
  tr/\n//d foreach (@tokens);   # Newlines are not significant
  for (my $i = 0; $i < @tokens - 1; ) { # Increment intentionally omitted
    if ($tokens[$i] !~ /^[a-zA-Z_]/) {
      carp "Expecting variable name, found ".$tokens[$i];
      $i++;
    }
    if ($tokens[$i+1] eq '=' and ($i+4 >= @tokens or $tokens[$i+4] eq '=')) { # Single value
      if ($tokens[$i+2] =~ /^[^']+\*/) { # Array, all same value
        $this->set($tokens[$i], [ $tokens[$i+2] ] );
      } else {                  # Scalar
        $this->set($tokens[$i],$tokens[$i+2]);
      }
      $i += 3;
    } elsif ($tokens[$i + 3] =~ /^[-'"0-9]/) { # Array
      my $j = $i+2;
      $j++ while ( $j < @tokens and $tokens[$j] ne '='); # Find end of array
      $j -= 2;
      $this->set($tokens[$i], [@tokens[$i+2..$j]]);
      $i = $j + 1;
    } else {
      carp "Can't process token: ".$tokens[$i];
      $i++;
    }
  }
  if ($tokens[@tokens - 1] ne '/') {
    carp "No end namelist token ('/') found";
  }
  return 1;
}

# Set a value in the perl data structure using a namelist specifier
# for the location
sub set {
  my $this = shift;
  my ($name, $value) = @_;
  croak "Need a name arg in set" unless $name ne "";
  my @tokens = split(/%/, uc $name);
  my $location = \$this->{data};
  foreach my $tok (@tokens) {
    if ($tok =~ /(\w+)\((\d+)\)/) { # if array element
      $location = \ $$location->{$1}[$2];
    } else {
      $location = \ $$location->{$tok};
    }
  }
  if (ref $value) {             # If array, expand repeats
    for (my $i = 0; $i < @$value; $i++) {
      if ($value->[$i] =~/^([^']+)\*(.+)$/) { # if entry has repeat
        my ($repeat, $val) = ($1, $2);
        splice @$value, $i, 1, ($val) x $repeat;
        $i += $repeat - 1;
      }
    }
  }
  $$location = $value;
}

# Return string representation (F90 style) of name list object
# (Wrapper for recursive _genEntryString)
#
sub namelistString {
  my $this = shift;
  my $str = join "", "&",$this->{name},"\n ";
  my $col = 1;
  $this->_genEntryString($this->{data},\$str,"", \$col, 1);
  $str .= "\n/\n";
  return $str;
}

# Private function that does real work of generating the namelist text
# from the perl data structure. Generates lists to intialize arrays
# where possible
#
sub _genEntryString {
  my $this = shift;
  my ($entry, $str, $prefix, $col, $toplevel) = @_;
  if (not ref $entry) {         # Simple scalar value, write output,
                                # inserting newline as necessary
    my $entryLen = length($prefix) + length($entry) + 3;
    if ($$col + $entryLen > 80) {
      $$str .= "\n ";
      $$col = 1;
    }
    $$col += $entryLen;
    $$str .= "$prefix=$entry, ";
  } elsif (ref $entry eq "HASH") { # Recurse to handle each hash element
                                   # separately, special case top level since
                                   # it doesn't represent a derived type
    foreach my $key (sort keys %$entry) {
      if ($toplevel) {
        $this->_genEntryString($entry->{$key}, $str, $prefix."$key", $col,0);
      } else {
        $this->_genEntryString($entry->{$key}, $str, $prefix."\%$key", $col,0);
        
      }
    }
  } elsif (ref $entry eq "ARRAY") {
    if ($prefix !~ /%/ and not ref $entry->[0]) { # Simple array, init with list
      $$str .= "\n " unless $$col == 1; 
      $$str .= "$prefix=";
      $$col = 2 + length($prefix);
      foreach my $el (@$entry) {
        my $entryLen = length($el) + 1;
        if ($$col + $entryLen > 80) {
          $$str .= "\n ";
          $$col = 1;
        }
        $$col += $entryLen;
        $$str .= "$el,";
      }
    } else {                    # Complex (derived type) array, init elements separately
      for (my $i = 0; $i < @$entry; $i++) {
        $this->_genEntryString($entry->[$i], $str, $prefix."($i)", $col,0);
      }
    }
  }
}

package CPinData;
no strict;

%CHRG = ("AS2" => [		# ASH backbone charges, extra q on CB
		  [0,		# Relative energy
		   0,		# Rel protonation
		   -0.415700,	# ASP
		   0.271900,
		   0.034100,
		   0.086400,
		   -0.031600-0.147,
		   0.048800-0.061,
		   0.048800-0.061,
		   0.646200+0.1532,
		   -0.555400-0.246,
		   -0.637600-0.1638,
		   0.474700-0.4747,
		   0.597300,
		   -0.567900,
		    0.474700-0.4747
		  ],
		  [26.4491 + 1.3818*4.0,	# Relative energy
		   1,		# Rel protonation
		   -0.415700,	# ASH syn
		   0.271900,
		   0.034100,
		   0.086400,
		   -0.031600,
		   0.048800,
		   0.048800,
		   0.646200,
		   -0.555400,
		   -0.637600,
		   0.474700,
		   0.597300,
		   -0.567900,
		   0.0
		  ],
		  [26.4491 + 1.3818*4.0,	# Relative energy
		   1,		# Rel protonation
		   -0.415700,	# ASH anti
		   0.271900,
		   0.034100,
		   0.086400,
		   -0.031600,
		   0.048800,
		   0.048800,
		   0.646200,
		   -0.555400,
		   -0.637600,
		   0.0,
		   0.597300,
		   -0.567900,
		   0.474700
		  ]
		 ],
	"AS4" => [		# ASH backbone charges, extra q on CB
		  [0,		# Relative energy
		   0,		# Rel protonation
		   -0.415700,	# ASP
		   0.271900,
		   0.034100,
		   0.086400,
		   -0.031600-0.147,
		   0.048800-0.061,
		   0.048800-0.061,
		   0.646200+0.1532,
		   -0.555400-0.246,
		   -0.637600-0.1638,
		   0.474700-0.4747,
		   0.597300,
		   -0.567900,
		    0.474700-0.4747,
		   0.0,
		   0.0
		  ],
		  [27.1769 + 1.3818*4.0-0.6*log(2),	# Relative energy
		   1,		# Rel protonation
		   -0.415700,	# ASH syn on O 2
		   0.271900,
		   0.034100,
		   0.086400,
		   -0.031600,
		   0.048800,
		   0.048800,
		   0.646200,
		   -0.555400,
		   -0.637600,
		   0.474700,
		   0.597300,
		   -0.567900,
		   0.0,
		   0.0,
		   0.0
		  ],
		  [27.1769 + 1.3818*4.0-0.6*log(2),	# Relative energy
		   1,		# Rel protonation
		   -0.415700,	# ASH anti on O 2
		   0.271900,
		   0.034100,
		   0.086400,
		   -0.031600,
		   0.048800,
		   0.048800,
		   0.646200,
		   -0.555400,
		   -0.637600,
		   0.0,
		   0.597300,
		   -0.567900,
		   0.474700,
		   0.0,
		   0.0
		  ],
		  [27.1769 + 1.3818*4.0-0.6*log(2),	# Relative energy
		   1,		# Rel protonation
		   -0.415700,	# ASH syn on O 1
		   0.271900,
		   0.034100,
		   0.086400,
		   -0.031600,
		   0.048800,
		   0.048800,
		   0.646200,
		   -0.637600,
		   -0.555400,
		   0.0,
		   0.597300,
		   -0.567900,
		   0.0,
		   0.4747,
		   0.0
		  ],
		  [27.1769 + 1.3818*4.0-0.6*log(2),	# Relative energy
		   1,		# Rel protonation
		   -0.415700,	# ASH anti on O 1
		   0.271900,
		   0.034100,
		   0.086400,
		   -0.031600,
		   0.048800,
		   0.048800,
		   0.646200,
		   -0.637600,
		   -0.555400,
		   0.0,
		   0.597300,
		   -0.567900,
		   0.0,
		   0.0,
		   0.4747
		  ]
		 ],
	"GL4" => [		# GLH backbone charges, extra q on CB
		  [0,		# Relative energy
		   0,		# Rel protonation
		   -0.4157,	# GLU (deprotonated)
		   0.2719,
		   0.0145,
		   0.0779,
		   -0.0398,
		   -0.0173,
		   -0.0173,	
		   0.0136,	
		   -0.0425,	
		   -0.0425,
		   0.8054,	
		   -0.8188,	
		   -0.8188,	
		   0.0,	
		   0.5973,	
		   -0.5679,	
		   0.0,
		   0.0,
		   0.0
		  ],
		  [8.93081-0.08144 + 1.3818*4.4-0.6*log(2+2*exp(-1.673/0.6)),	# Relative energy 
		   1,		# Rel protonation
		   -0.4157,	# GLH syn O 2
		   0.2719,	
		   0.0145,	
		   0.0779,	
		   -0.0071,	
		   0.0256,	
		   0.0256,	
		   -0.0174,	
		   0.043,	
		   0.043,	
		   0.6801,	
		   -0.5838,	
		   -0.6511,	
		   0.4641,	
		   0.5973,	
		   -0.5679,
		   0.0,
		   0.0,
		   0.0
		  ],
		  [8.93081-0.08144 + 1.3818*4.4-0.6*log(2+2*exp(-1.673/0.6)),	# Relative energy 
		   1,		# Rel protonation
		   -0.4157,	# GLH anti O 2
		   0.2719,	
		   0.0145,	
		   0.0779,	
		   -0.0071,	
		   0.0256,	
		   0.0256,	
		   -0.0174,	
		   0.043,	
		   0.043,	
		   0.6801,	
		   -0.5838,	
		   -0.6511,	
		   0.0,	
		   0.5973,	
		   -0.5679,
		   0.4641,
		   0.0,
		   0.0
		  ],
		  [8.93081-0.08144 + 1.3818*4.4-0.6*log(2+2*exp(-1.673/0.6)),	# Relative energy 
		   1,		# Rel protonation
		   -0.4157,	# GLH syn O 1
		   0.2719,	
		   0.0145,	
		   0.0779,	
		   -0.0071,	
		   0.0256,	
		   0.0256,	
		   -0.0174,	
		   0.043,	
		   0.043,	
		   0.6801,	
		   -0.6511,	
		   -0.5838,	
		   0.0,	
		   0.5973,	
		   -0.5679,
		   0.0,
		   0.4641,
		   0.0
		  ],
		  [8.93081-0.08144 + 1.3818*4.4-0.6*log(2+2*exp(-1.673/0.6)),	# Relative energy 
		   1,		# Rel protonation
		   -0.4157,	# GLH anti O 1
		   0.2719,	
		   0.0145,	
		   0.0779,	
		   -0.0071,	
		   0.0256,	
		   0.0256,	
		   -0.0174,	
		   0.043,	
		   0.043,	
		   0.6801,	
		   -0.6511,	
		   -0.5838,	
		   0.0,	
		   0.5973,	
		   -0.5679,
		   0.0,
		   0.0,
		   0.4641
		  ]
		 ],
	"HIP" => [
		  [0, # Energy correction from MD titration
		   2,		# Protonation
		   -0.347900,	# HIP
		   0.274700,
		   -0.135400,
		   0.121200,
		   -0.041400,
		   0.081000,
		   0.081000,
		   -0.001200,
		   -0.151300,
		   0.386600,
		   -0.017000,
		   0.268100,
		   -0.171800,
		   0.391100,
		   -0.114100,
		   0.231700,
		   0.734100,
		   -0.589400,
		  ],
		  [-2.74 - 1.3818*6.5 - 0.06,	# Energy (from TI)
		   1,		# protonation
		   -0.347900,	# HID
		   0.274700,
		   -0.135400,
		   0.121200,
		   -0.111,
		   0.040200,
		   0.040200,
		   -0.026600,
		   -0.381100,
		   0.364900,
		   0.205700,
		   0.139200,
		   -0.572700,
		   0.00000,
		   0.129200,
		   0.114700,
		   0.734100,
		   -0.589400,
		  ],
		  [-6.22 - 1.3818*7.1 - 0.24,	# Energy (from TI)
		   1,		# prot
		   -0.347900,	# HIE
		   0.274700,
		   -0.135400,
		   0.121200,
		   -0.1012,
		   0.036700,
		   0.036700,
		   0.186800,
		   -0.543200,
		   0.00000,
		   0.163500,
		   0.143500,
		   -0.279500,
		   0.333900,
		   -0.220700,
		   0.186200,
		   0.734100,
		   -0.589400,
		  ]
		 ],
	"LYS" => [
		  [-15.26 +1.3818*10.4 + 0.05, 		# Relative energy
		   3,		# Protonation
		   -0.3479,	# LYS
		   0.2747,
		   -0.24,
		   0.1426,
		   -0.0094,
		   0.0362,
		   0.0362,
		   0.0187,
		   0.0103,
		   0.0103,
		   -0.0479,
		   0.0621,
		   0.0621,
		   -0.0143,
		   0.1135,
		   0.1135,
		   -0.3854,
		   0.34,
		   0.34,
		   0.34,
		   0.7341,
		   -0.5894
		  ],
		  [0,		# Energy correction from MD titration
		   2,		# Protonation
		   -0.3479,	# LYN
		   0.2747,
		   -0.24,
		   0.1426,
		   -0.10961,
		   0.034,
		   0.034,
		   0.06612,
		   0.01041,
		   0.01041,
		   -0.03768,
		   0.01155,
		   0.01155,
		   0.32604,
		   -0.03358,
		   -0.03358,
		   -1.03581,
		   0,
		   0.38604,
		   0.38604,
		   0.7341,
		   -0.5894
		  ]
		 ],
	 "TYR" => [
		   [0,		# Relative energy
		    1,		# Protonation
		    -0.4157,	# TYR
		    0.2719,
		    -0.0014,
		    0.0876,
		    -0.0152,
		    0.0295,
		    0.0295,
		    -0.0011,
		    -0.1906,
		    0.1699,
		    -0.2341,
		    0.1656,
		    0.3226,
		    -0.5579,
		    0.3992,
		    -0.2341,
		    0.1656,
		    -0.1906,
		    0.1699,
		    0.5973,
		    -0.5679,
		   ],
		   [-65.0791 - 1.3818*9.6,		# Energy
		    0,		# Protonation
		    -0.4157,	# TYM
		    0.2719,
		    -0.0014,
		    0.0876,
		    -0.0858,
		    0.019,
		    0.019,
		    -0.213,
		    -0.103,
		    0.132,
		    -0.498,
		    0.132,
		    0.777,
		    -0.814,
		    0,
		    -0.498,
		    0.132,
		    -0.103,
		    0.132,
		    0.5973,
		    -0.5679,
		   ]
		  ]
       );

# HEWL pKas from Improving the Continuum Dielectric Approach to Calculating pKas of Ionizable Groups in Proteins
# Demchuk E, Wade RC
# J. Phys. Chem. 1996; 100: 17373-17387
%PKA = ( "HEWL" => {15 => 5.71,	# HIS ##
		    7 => 2.85,	# GLU
		    35 => 6.2,
		    18 => 2.66,	# ASP
		    48 => 2.5, ##
		    52 => 3.68, ##
		    66 => 2.0, ##
		    87 => 2.07,
		    101 => 4.09,
		    119 => 3.2,
		    20 => 10.3,	# TYR
		    23 => 9.8,
		    53 => 12.1,
		    1 => 10.8,	# LYS
		    13 => 10.5,
		    33 => 10.6,
		    96 => 10.8,
		    97 => 10.3,
		    116 => 10.4
		   }
       );

1;
