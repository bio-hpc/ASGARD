BEGIN{i=1;s=""}
{while ( i < NF ){s = s" "$i;if( $i == "-np" || $i == "-n"){i++;s=s" 4";}; i++;}}
END{print s}
