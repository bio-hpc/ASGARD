#/usr/bin/perl

# Script for the determination of the STOP criterion for TI production.
############################################################################################
# Case 1 (method=1):
# The average dV/dL standard error for a specific lambda and the
# difference between the current standard error and the standard error calculated
# in the previous step is determined.
# Depending on the specified error limit "GO_ON" or "STOP" are output to the standard error
# file. Furthermore if the standard error significantly increases from one step to
# another "CAUTION" is printed to the output file to ensure that a warning will be issued.
#
# The autocorrelation time is determined based on all dV/dL values calculated during the
# first production run. If the number of dV/dL values is smaller than 10 x Tau, STOP is 
# printed to standard error file and the analyses will be terminated after the first round
# of production simulations.
##########################################################################################
# Case 2 (method=2):
# The precision of average dV/dL is determined employing the student's distribution.
# A significance level of 95% is used.
# If the determined precision is below the specified limit, "STOP" is written to the
# standard error file.
##########################################################################################
use FindBin qw($Bin);
use lib "$Bin/../additional_libs";
use PointEstimation;
use Distributions;
use Descriptive;


my $path = $ARGV[0];
my $alias = $ARGV[1];
my $lambda = $ARGV[2];
my $no_of_files = $ARGV[3];
my $method = $ARGV[4];
my $limit = $ARGV[5];

if(@ARGV < 5){
	print "Usage: perl convergenceCheck.pl <path> <alias of V0 structure> <lambda value> <no of files that shall be regarded> <method> <error limit>\n";
	exit;
}

# Reading of last standard error
my $stderr_last=0;
my $autocorr=0; 
my $t_file = "stderr_".$lambda;
if(-e $t_file){
	open(STDERR_LAST, $t_file) || die "Cannot open file $t_file for reading.\n";
	my $status = <STDERR_LAST>;
	$stderr_last = <STDERR_LAST>;
	chomp($stderr_last);
	$stderr_last =~ s/^\s+//;
	$stderr_last =~ s/\s+$//;
	$autocorr = <STDERR_LAST>;
	if($autocorr =~ m/Autocorrelation: (\d+)/){
		$autocorr = $1;
	}
}

# Read dV/dL values and calculate average, standard deviation and standard error
my $sum_dVdL = 0;
my $avg_dVdL = 0;
my @dVdL_data = ();
my $value_no = 0;
my $begin_read = 0;
my $time_interval_detected = 0;
my $save_time = 0;
my $time_inter = 0;

chdir($path);

# Read data and calculate average dV/dL
for(my $i=1; $i<=$no_of_files; $i++){
	
	my $p;
	if($i < 10){
		$p = "0$i";
	}
	else{
		$p = $i;
	}
	
	my $file_i = $path."/".$alias."_prod".$p."_v0_l".$lambda.".out";
	open(FILE, $file_i) || die "Cannot open $file_i for reading.\n";

	while(my $f_line = <FILE>){
		chomp($f_line);
			
		if($f_line =~ m/4\.  RESULTS/){
			$begin_read = 1;
		}
			
		if($f_line =~ m/A V E R A G E S/){
			$begin_read = 0;
			close FILE;
			last;
		}
		
		if(($begin_read == 1)&&($time_interval_detected == 0)&&($f_line =~ m/NSTEP =\s+(\d+)   TIME\(PS\) =\s+(\d+)\.(\d+)/)){
			my $time = $2.".".$3;
			
			if($save_time != 0){
				$time_inter = $time - $save_time;
				$time_interval_detected = 1;
			}
			$save_time = $time;
		}
			
		if(($begin_read == 1)&&($f_line =~ m/^ DV\/DL/)){
			my @value = split(/\s+/, $f_line);
			$sum_dVdL = $sum_dVdL + $value[3];
			push(@dVdL_data, $value[3]);
			$value_no++;
		}
	}
}

$avg_dVdL = $sum_dVdL / $value_no;
$avg_dVdL = sprintf("%.4f", $avg_dVdL);
$sum_dVdL = 0;


# Determine standard deviation
my $std_sum = 0;
my $value_no_std = 0;
my $std_dVdL;
my $std_err;

for(my $i=0; $i<@dVdL_data; $i++){
	$std_sum = ($std_sum + (($dVdL_data[$i] - $avg_dVdL)**2));
	$value_no_std++;
}

my $factor = (1 / ($value_no - 1));
$std_dVdL = $factor * $std_sum;
$std_dVdL = sqrt($std_dVdL);

# Determine total simulation time
my $total_time = $time_inter * $value_no;

# Check autocorrelation time in first run
my @res;
my $ratio;
if(! -e $t_file){
	my @data1 = @dVdL_data;
	my @data2 = @dVdL_data;
	$autocorr = determine_autocorrelation(\@data1, \@data2, \@res);
}

# Check whether the total number of dV/dL values is smaller than 10 x Tau
$ratio = $value_no / $autocorr;

# Determine standard error
my $factor_autocorr = $total_time / (2*$autocorr*$time_inter);
$std_err = 1/sqrt($factor_autocorr);
$std_err = $std_err*$std_dVdL;
$std_err = sprintf("%.3f", $std_err);

# Calculate difference
my $stderr_diff = $stderr_last - $std_err;
my $abs_diff = abs($stderr_diff);

# Determine precision according to student's distribution using a significance level of 95%
if($autocorr > 1){
	my @Y;
	my $count=0;
	my $sum=0;
	foreach my $d (@dVdL_data){
		$count++;
		$sum = $sum+$d;
		if($count == $autocorr){
			my $avg = $sum / $count;
			$count = 0;
			$sum = 0;
			push(@Y, $avg);
		}
	}
	@dVdL_data = @Y;
}

my $mue = determine_precision(\@dVdL_data);
$mue = sprintf("%.3f", $mue);
	
open(OUT, ">$t_file") || die "Cannot open file $t_file for writing.\n";
if(($method == 1)&&($ratio < 10)){
	print OUT "RATIO_TOO_SMALL\n";
}
elsif((($method == 1)&&($stderr_diff > 0)&&($abs_diff <= $limit))||
	  (($method == 2)&&($mue < $limit))){
	print OUT "STOP\n";
}
elsif(($method == 1)&&($stderr_last!=0)&&($stderr_diff < 0)&&($abs_diff > 0.05)){
	print OUT "CAUTION\n";
}
else{
	print OUT "GO_ON\n";
}
print OUT $std_err."\n";
print OUT "Autocorrelation: ".$autocorr."\n";
print OUT "Convergence level reached: ".$mue." kcal/mol";
close OUT;



# Subroutine for the determination of the autocorrelation time
# Original time correlation function written by H. Gohlke
sub determine_autocorrelation{
	my $autocorrelation_time;
	my $r_data1 = shift;
	my $r_data2 = shift;
	my $r_res   = shift;
	@$r_res = ();

	my $len1 = scalar(@$r_data1);
	my $len2 = scalar(@$r_data2);
	my $len;
	if ($len1 > $len2){
		$len = $len1;
	}
	else{
		$len = $len2;
	}

	# data1 lags data2 (i.e. is shifted to the right of it)
	# -> peak in positive lags (i.e. 0 .. $len-1 in @res)
	# data2 lags data1 (i.e. is shifted to the right of it)
	# -> peak in negative lags (i.e. $len .. 2*$len-1 in @res)
	#
	# * data1/2 are centered
	# * correlation values are finally normalized
	my $i;
	for($i = 0;$i < 2*$len; $i++){
		push @$r_res,0.0;
	}
	my $mean1 = 0.0;
	for($i=0; $i < $len1; $i++){
		$mean1 += $r_data1->[$i];
	}
	$mean1 /= $len1;

	my $mean2 = 0.0;
	for($i=0; $i < $len2; $i++){
		$mean2 += $r_data2->[$i];
	}
	$mean2 /= $len2;

	my $lag;
	my $norm1 = 0;
	my $norm2 = 0;
	my $norm;
	for($lag = 0; $lag < $len; $lag++){
		for($i = 0; $i < $len; $i++){
			$r_res->[$lag] +=	($i+$lag >= $len1 ? 0.0 : $r_data1->[$i+$lag]-$mean1) *
								($i      >= $len2 ? 0.0 : $r_data2->[$i]-$mean2);
			$r_res->[$len+$lag] +=	($i      >= $len1 ? 0.0 : $r_data1->[$i]-$mean1) *
									($i+$lag >= $len2 ? 0.0 : $r_data2->[$i+$lag]-$mean2);
			if($lag == 0){
				$norm1 += ($r_data1->[$i]-$mean1) * ($r_data1->[$i]-$mean1);
				$norm2 += ($r_data2->[$i]-$mean2) * ($r_data2->[$i]-$mean2);
			}
		}
		if($lag == 0){
			$norm = sqrt($norm1 * $norm2);
		}
		$r_res->[$lag] /= $norm;
		$r_res->[$len+$lag] /= $norm;
	}
	
	my $e_factor = 1/exp(1);
	for($lag = 0; $lag < $len; $lag++){
		if($r_res->[$lag] < $e_factor){
			$autocorrelation_time = $lag;
			last;
		}
	}
	
	return $autocorrelation_time;
}


# Determination of convergence level, i.e. +/-µ in kcal/mol
sub determine_precision{
	my $dVdL_data_ref = shift;
	my @dVdL_data = @{$dVdL_data_ref};
	
	my $stat = new Statistics::PointEstimation;
	$stat->set_significance(95); #set the significance(confidence) level to 95%
	$stat->add_data(@dVdL_data);
	my $mue = $stat->delta();
	return $mue;
}
