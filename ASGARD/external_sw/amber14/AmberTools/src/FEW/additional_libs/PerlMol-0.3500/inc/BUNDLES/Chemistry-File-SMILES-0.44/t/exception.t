use strict;
use warnings;

use Test::More;
use Chemistry::File::SMILES;

BEGIN {
    if (eval 'use Test::Exception; 1') {
        #plan tests => 9;
        plan 'no_plan';
    } else {
        plan skip_all => "You don't have Test::Exception installed";
    }
}

# make sure that we throw some expected exceptions

###### MOL ######

throws_ok { Chemistry::Mol->parse('CC=(O)', format => 'smiles') } 
    qr/SMILES ERROR/, "bond before branch";

throws_ok { Chemistry::Mol->parse('CxC(=O)', format => 'smiles') } 
    qr/SMILES ERROR/, "invalid symbol";



