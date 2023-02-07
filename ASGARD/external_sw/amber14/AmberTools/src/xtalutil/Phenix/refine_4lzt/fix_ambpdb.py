import iotbx.pdb

#======================================================================#
# When tleap creates an Amber topology and ambpdb recreates a pdb with #
# atom order that matches the topology, information about chain name,  #
# b-factors and occupancy is lost. This is a quick jiffy script (will  #
# most likely fail on more complex cases like multiple residues with   #
# same seq id) to replace the missing information using the original   #
# file.																   #
########################################################################

pre_amber_pdb_file="4lzti.pdb"
post_amber_pdb_file="4lzt.pdb"
output_pdb_file="out.pdb"


pdb_pre=iotbx.pdb.input(file_name=pre_amber_pdb_file)
pdb_h_pre = pdb_pre.construct_hierarchy()
pdb_post=iotbx.pdb.input(file_name=post_amber_pdb_file)
pdb_h_post = pdb_post.construct_hierarchy()

for chain_post in pdb_h_post.chains():
  for resi_post in chain_post.conformers()[0].residues():
      for atom_post in resi_post.atoms():
        #import code; code.interact(local=locals()); assert 0  
        #print resi_post.resname
        #if resi_post.resname=='WAT': print chain_post.id, resi_post.resseq, resi_post.resname
        for chain_pre in pdb_h_pre.chains():
          for resi_pre in chain_pre.conformers()[0].residues():
            if resi_post.resname=='NO3' and resi_pre.resname=='NO3':
				print resi_pre.resseq, resi_post.resseq
			
            if resi_pre.resseq==resi_post.resseq and resi_pre.resname==resi_post.resname:
              
              for atom_pre in resi_pre.atoms():
                if atom_pre.name == atom_post.name and atom_pre.i_seq==atom_post.i_seq:
					atom_post.b=atom_pre.b
					atom_post.occ=atom_pre.occ
					chain_post.id=chain_pre.id	
pdb_h_post.write_pdb_file(file_name=output_pdb_file, append_end=True,crystal_symmetry=pdb_pre.crystal_symmetry())
