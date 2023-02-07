use Test::More;

use Chemistry::File::PDB;

my $mol = Chemistry::Mol->read("test.pdb");
my $macromol = Chemistry::MacroMol->read("test.pdb");

plan tests => 9;

is(scalar($mol->atoms),             139,        '$mol->atoms');
is(scalar($macromol->atoms),        139,        '$macromol->atoms');
is(scalar($macromol->domains),      10,         '$macromol->domains');
is($macromol->domains(4)->type,     'VAL',      '$macromol->domains(4)->name');
is($macromol->domains(4)->name,     'VAL4',     '$macromol->domains(4)->name');
is($mol->atoms(48)->name,           'CG1',      '$mol->atoms(44)->name');
is($mol->atoms(48)->symbol,         'C',        '$mol->atoms(44)->name');
is($mol->atoms(48)->attr("pdb/residue_name"),   'VAL4',  
    '$mol->atoms(44)->attr("pdb/residue_name")');
is($macromol->domains(4)->atoms_by_name('CG1')->attr("pdb/serial_number"),
    48,  
    q{$macromol->domains(4)->atoms_by_name('CG1')->attr("pdb/serial_number")});


