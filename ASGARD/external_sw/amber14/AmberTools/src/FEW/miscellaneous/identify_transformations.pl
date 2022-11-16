#/usr/bin/perl

# Script for the identification of the optimal transformations for a set of compounds
# based on Kruskal's algorithm. This program uses the perl package Graph::Kruskal available
# at CPAN (http://www.cpan.org). As input it requires a file containing a matrix of similarity 
# scores for the compounds that shall be analyzed in Tab-separated format. It is assumed that
# the smaller the score, the larger the similarity between the compounds. The first line and
# the first column of the matrix should contain the names of the structures to be analyzed.

use FindBin;
my $current_dir = $FindBin::Bin;
my @split_current_dir = split(/\//, $current_dir);
pop(@split_current_dir);
our $assumed_FEW_dir .= join("/", @split_current_dir);

# Check presence of Graph::Kruskal module
my @add_modules = ("Kruskal");
my @names_add_modules = ("Graph::Kruskal");

eval "use Graph::Kruskal qw(:all); 1" or die "\nThe Perl-module '$names_add_modules[$module_n]' required for running identify_transformations.pl
was not found. This module can be obtained under the terms of the same licence 
as Perl itself from http://www.cpan.org
Please install the module '$names_add_modules[$module_n]' before executing this script.\n\n";


if(@ARGV < 2){
	print "Usage: kruskal.pl <number of structures to regard> <score matrix file>\n";
	exit;
}

# Define vortices
my $vortices_no = $ARGV[0];

if($vortices_no == 0){
	print "Please specify the number of vertices, i.e. structures, to regard.\n";
	exit;
}

my @v_str;
foreach my $v (1...($vortices_no)){
	push(@v_str, $v);
}
&define_vortices(@v_str);

# Read the score matrix entries
my $matrix_file = $ARGV[1];

if($matrix_file eq ""){
	print "Please specify the name of the file containing the score matrix.\n";
	exit;
}

my $line_count=0;
my $col_count=0;
my @edges;
my %struct_names;

open(SCORE_MAT, "$matrix_file")||die "Cannot open the score file $matrix_file for reading.\n";
my $header = <SCORE_MAT>;

while(my $s_line = <SCORE_MAT>){
	$line_count++;
	chomp($s_line);
	
	my @scores = split(/\t/, $s_line);
	$struct_names{$line_count} = $scores[0];
	
	foreach my $s (1...$#scores){
		$col_count = $s;
		if($col_count > $line_count){
			push(@edges, $line_count);
			push(@edges, $col_count);
			push(@edges, $scores[$s]);
		}
	}
	
	$col_count = 0;
}

&define_edges(@edges);
my @kruskal = &kruskal();

open(KRUSK, ">kruskal.out") || die "Cannot open file kruskal.out\n";
print KRUSK "# List of transformations identified as optimal by Kruskal's algorithm: #\n";
print KRUSK "#########################################################################\n\n";
print KRUSK "Tansformation\t\tSimilarity score\n";
foreach my $k (@kruskal){
	if($k){
		my %krusk = %{$k};	
		print KRUSK $struct_names{$krusk{'from'}}."<->".$struct_names{$krusk{'to'}}."\t\t".$krusk{'cost'}."\n";
	} 
}
close KRUSK;
