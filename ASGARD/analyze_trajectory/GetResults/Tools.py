#
#	Class with several functions about distance, checking data, split ligand-protein...
#
import re, math, os
import shutil

from .Molecule import Molecule
from .Execute import Execute

class Tools(object): 
    def __init__(self, cfg):
        self.cfg = cfg
        self.lstTemplate = []
        self.execute = Execute(self.cfg)

    def distancia(self, array_lig, array_prot):
        """
            Distance between 2 points
        """
        x = math.pow(array_lig[0] - array_prot[0], 2)
        y = math.pow(array_lig[1] - array_prot[1], 2)
        z = math.pow(array_lig[2] - array_prot[2], 2)
        ret = math.sqrt(x + y + z)
        return ret

    def create_hash(self, residuo, distance, lig_atom, hashCarboHidrato):
        # A table is generated with the residues which are found less than x A
        if residuo not in hashCarboHidrato:
            hashCarboHidrato[residuo] = [distance, lig_atom]
        else:
            if distance < hashCarboHidrato[residuo][0]:
                hashCarboHidrato[residuo] = [distance, lig_atom]
        return hashCarboHidrato

    def get_residues(self, mol_target, mol_query, longitud_maxima):
        """
            To obtain the residues between a ligand and a protein
        """
        hash_carbo_hidrato = {}
        for i in mol_query.coords:
        
            aux_query = re.sub(' +',' ',i).strip().split(" ")
            for j in mol_target.coords:
                
                aux_target = re.sub(' +',' ',j).strip().split(" ")
                dist = self.distancia([float(aux_query[5]),float(aux_query[6]),float(aux_query[7])],[float(aux_target[5]),float(aux_target[6]),float(aux_target[7])]  )

                if dist < longitud_maxima:
                    if mol_target.name == "Protein" or mol_target.name == "DNA":
                        hash_carbo_hidrato = self.create_hash( aux_target[3] + "_" + aux_target[4], dist, aux_query[2], hash_carbo_hidrato )
                    else:
                        hash_carbo_hidrato = self.create_hash(aux_target[3] , dist, aux_query[2],
                                                              hash_carbo_hidrato)
        return hash_carbo_hidrato

    def get_about_residues(self, mol_target, mol_query):
        #
        # This function starts with 20 A. After that, the residues closest to the ligand are searched
        # until they are more than 64, which is the maximum group of gromacs
        #
        print("Entro")
        cad_residues = ""

        for i in range(self.cfg.DIST_MIN_RES):
            cad_residues = mol_query.name
            aux = self.get_residues(mol_target, mol_query, self.cfg.DIST_MIN_RES - i)
            for j in aux:
                cad_residues += " " + j
            # There is a limit of 64 groups
            if len(cad_residues.split(" ")) < self.cfg.NUM_MAX_GROUPS:
                break
        return cad_residues
    
    def cp_file(self, src, dest):
        shutil.copy(src, dest)

    def get_groups_target_queries(self):
        """
            An index is created, the ligand's itps are searched at the top, and then they are searched in the index by name
        """
        lst_molecules = []

        index = self.cfg.prefix_molec + '_index.ndx.tmp'
        cmd = f"echo q | {self.cfg.gromacs} make_ndx{self.cfg.mpi} -f {self.cfg.gro_md} -o {index}.tmp.ndx > {index} 2> /dev/null"
        
        self.execute.run(cmd)
        
        if self.cfg.profile != "TARGET":
            cmd = f'cat {self.cfg.top} | grep \\#include | grep -v force_field | grep -v _[a-z].itp | grep -v porse | grep -v posre | grep -v ^\; | grep -v forcefield | grep -v ions | grep -v tip3p | grep itp'
            lst = self.execute.run(cmd).split("\n")
            for i in lst:
                if i != "":
                    aux = i.split("\"")
                    print(aux[1])
                    if os.path.isfile(aux[1]):
                        aux[1]=aux[1]
                    else:
                        print(self.cfg.folder)
                        print(self.cfg.folder+'molecules/'+aux[1])
                        aux[1]=self.cfg.folder+'molecules/'+aux[1]
                
                    cmd = f'cat {aux[1]} | grep "\\[ atoms \\]" -A3  | tail -n 1 | awk \'{{print $4}}\''
                    print(cmd)
                    n_query = self.execute.run(cmd).strip()
                    print(n_query)
                    cmd = f'cat {index} | grep {n_query} | grep -v "non-Protein" | grep -v "non-Water" | head -1 | awk \'{{print $1}}\''
                    g_query = self.execute.run(cmd).strip()
                    cmd = f'cat {aux[1]} | grep "\\[ moleculetype \\]" -A2  | tail -n 1 | awk \'{{print $1}}\''
                    name_query = self.execute.run(cmd).strip()
                    if name_query != "":
                        lst_molecules.append(Molecule(n_query, g_query, name_query))
        cmd = f'cat {index}  | grep " DNA " | tail -1 | awk \'{{print $1}}\''
        g_dna = self.execute.run(cmd).strip()
        if g_dna != "":
            lst_molecules.insert(0, Molecule('DNA', g_dna, "DNA"))

        cmd = f'cat {index}  | grep " Protein " | tail -1 | awk \'{{print $1}}\''
        g_target = self.execute.run(cmd).strip()
        if g_target != "":
            #
            #   When the protein is the first group
            #
            lst_molecules.insert(0, Molecule('Protein', g_target, self.cfg.name_target ))

        return lst_molecules

    def split_queries(self):
        """
            Generates an array with the protein and a dictionary with the ligands with their lines extracted from the pdv file

        """
        f = open(self.cfg.pdb)
        target = self.cfg.lst_molecules[0]
        for i in f:
            lig = False
            aux = re.sub(' +', ' ', i).strip().split(" ")

            if len(aux) > 0:    # Discarding solvents and ions
                if aux[0] == "ATOM" and aux[3] != self.cfg.name_solvent and aux[3] not in self.cfg.name_ions:
                    for mol in self.cfg.lst_molecules:
                        if aux[3] == mol.name:
                            mol.coords.append(i.strip())
                            lig = True
                            break
                    if not lig:
                        target.coords.append(i.strip())

    def generate_distribution_xvg(self, xvg_in, xvg_out):
        cmd = f'{self.cfg.gromacs} {self.cfg.graph}analyze{self.cfg.mpi} -f {xvg_in} -dist {xvg_out} -bw {self.cfg.DISTRIBUTION_STEP}'
        self.cfg.tools.execute.run(cmd)

    def generate_xvg(self, n_group_1, n_group_2, type_graph, out_xtc):
        cmd = f'echo  "{n_group_1} {n_group_2}" | {self.cfg.gromacs} {self.cfg.graph}{type_graph}{self.cfg.mpi} -f {self.cfg.xtc_md} -s {self.cfg.tpr_min} -o {out_xtc}'
        self.cfg.tools.execute.run(cmd)

    def join_image(self, img_a, img_b, out):
        """
            join 2 images
        """
        comando = f"convert -size 2048x1536 {img_a} -thumbnail 800x600 {img_a}.png"
        self.execute.run(comando)
        comando = f"convert -size 2048x1536 {img_b} -thumbnail 800x600 {img_b}.png"
        self.execute.run(comando)
        comando = f"montage  -geometry +3+2  {img_a}.png {img_b}.png -tile 2x1 {out}"
        self.execute.run(comando)
        comando = f"rm {img_a}.png {img_b}.png"
        self.execute.run(comando)

    @staticmethod
    def check_directory(path):
        if os.path.exists(path):
            shutil.rmtree(path)
            os.mkdir(path)
        else:
            os.mkdir(path)
