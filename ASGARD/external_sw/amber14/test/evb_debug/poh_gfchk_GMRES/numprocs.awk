BEGIN{i=0}
{while ( i < NF ){if( $i == "-np" || $i == "-n"){i++;print $i; exit}; i++;}}
