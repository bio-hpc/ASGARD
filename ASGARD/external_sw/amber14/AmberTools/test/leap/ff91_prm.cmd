logfile ff91_prm.log

set default OldPrmtopFormat on
source oldff/leaprc.ff91

logfile ff91_prm.log

x = loadpdb ff91/all_aminoan.p
saveamberparm x ./all_aminoan91.top ./all_aminoan91.crd

ncres = { NALA CALA NPRO CPRO }
x = loadpdbusingseq ff91/all_aminonc.p ncres
saveamberparm x ./all_aminonc91.top ./all_aminonc91.crd

strand = { HB RADE RPOM RURA RPOM RCYT RPOM RGUA HE }
x = loadpdbusingseq ff91/all_rna91.p strand
# adjust charges of terminal residues
set x.2.O5' charge -0.425
set x.8.O3' charge -0.514

#set x.2.O5' type OH
#set x.8.O3' type OH

saveamberparm x ./all_rna91.top ./all_rna91.crd

strand = { HB DCYT DPOM DGUA DPOM DADE DPOM DTHY HE }
x = loadpdbusingseq ff91/all_dna91.p strand
set x.2.O5' charge -0.425
set x.8.O3' charge -0.514
saveamberparm x ./all_dna91.top ./all_dna91.crd

quit
