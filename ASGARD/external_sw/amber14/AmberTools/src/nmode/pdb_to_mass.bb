#! /bin/sh
#
#  read in pdb file and output "mass" file for use in
#    mdovrly
#
awk '
BEGIN {
	m["H"] = 1.008; m["C"] = 12.0; m["N"] = 14.0; m["O"] = 16.0;
	m["S"] = 32.0; m["Z"] = 65.38; m["F"] = 56.0; zero=0.0
}
$1=="ATOM" && ($3=="N"||$3=="C"||$3=="CA") { code = substr($3,1,1); 
	if (code in m) { n++; printf "%7.2f ", m[code]
		if (n==10) { n=0; printf "\n" }
	} else { print "error in mass parsing",code }
	next
}
{  n++; printf "%7.2f ", zero; if (n==10) { n=0; printf "\n" }
} ' $1 > $2
	
