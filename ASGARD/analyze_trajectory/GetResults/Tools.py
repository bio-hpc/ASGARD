#
#	Class with several functions about distance, checking data, split ligand-protein...
#
import re, math, os
import shutil

from .Molecule import Molecule
from .Execute import Execute

class Tools(object): 
    def __init__(self,cfg):
        self.cfg = cfg
        self.lstTemplate = []
        self.execute = Execute(self.cfg)

    def distancia(self,array_lig, array_prot):
        """
            Distance between 2 points
        """
        x = math.pow(array_lig[0] - array_prot[0], 2)
        y = math.pow(array_lig[1] - array_prot[1], 2)
        z = math.pow(array_lig[2] - array_prot[2], 2)
        ret = math.sqrt(x + y + z)
        return ret

    def create_hash(self,residuo, distance, lig_atom,hashCarboHidrato):
        # A table is generated with the residues which are found less than x A
        if residuo not in hashCarboHidrato:
            hashCarboHidrato[residuo] = [distance, lig_atom]
        else:
            if distance < hashCarboHidrato[residuo][0]:
                hashCarboHidrato[residuo] = [distance, lig_atom]
        return hashCarboHidrato

    def get_residues(self,mol_target, mol_query, longitud_maxima):
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
        print ("Entro")
        cad_residues = ""

        for i in range(self.cfg.DIST_MIN_RES):
            cad_residues = mol_query.name
            aux = self.get_residues(mol_target, mol_query, self.cfg.DIST_MIN_RES - i)
            for j in aux:
                cad_residues += " " + j
            # hasta que los grupos no sean meores que nomMaximoGrupos (Exite un limite en gromacs de 64 grupos)
            if len(cad_residues.split(" ")) < self.cfg.NUM_MAX_GROUPS:
                break
        return cad_residues
    def cp_file(self, src, dest):
        shutil.copy(src, dest)

    def get_groups_target_queries(self):
        """
            Se crea un indice, se busca en el top los itps del ligando y luego con el nombre se buscan en el indice
        """
        lst_molecules = []

        index = self.cfg.prefix_molec+'_index.ndx.tmp'
        cmd = "echo q | {0} make_ndx{1} -f {2} -o {3} > {4} 2> /dev/null".format(self.cfg.gromacs, self.cfg.mpi, self.cfg.gro_md, index+".tmp.ndx", index)

        self.execute.run(cmd)
        if self.cfg.profile != "TARGET":
            cmd = 'cat {} |grep \#include |grep -v force_field | grep -v _[a-z].itp |grep -v porse |grep -v ^\; |grep -v forcefield |grep -v ions |grep -v tip3p'.format(self.cfg.top)
            lst = self.execute.run(cmd).split("\n")
            for i in lst:
                if i != "":
                    aux = i.split("\"")
                    cmd = 'cat {0} |grep  \"\[ atoms \]\" -A3  |tail -n 1 |awk \'{{print $4}}\''.format(aux[1])
                    n_query = self.execute.run(cmd).strip()
                    cmd = 'cat {} | grep {} | grep -v "non-Protein" | grep -v "non-Water" | head -1 |awk \'{{print $1}}\''.format(index, n_query)  # Buscamos el numquery del ligando (a veces se genera un mol2 con nombre non
                    g_query = self.execute.run(cmd).strip()                                                                                        # Por eso eliminamos non-Water y non-Protein para que no interfiera
                    cmd = 'cat {0} |grep  \"\[ moleculetype \]\" -A2  |tail -n 1 |awk \'{{print $1}}\''.format(aux[1])
                    name_query = self.execute.run(cmd).strip()
                    if name_query != "":
                      lst_molecules.append(Molecule(n_query, g_query, name_query ))
                      #lst_molecules.append(Molecule(13, g_query, 'L01' ))
        cmd = 'cat {}  |grep \" DNA \" |tail -1 |awk \'{{print $1}}\''.format(index)  # buscamos el numgrupo del ADN
        g_dna = self.execute.run(cmd).strip()
        if g_dna != "":
            lst_molecules.insert(0, Molecule('DNA', g_dna, "DNA"))

        cmd = 'cat {}  |grep \" Protein \" |tail -1 |awk \'{{print $1}}\''.format(index)  # buscamos el numgrupo de la proteina
        g_target = self.execute.run(cmd).strip()
        if g_target != "":
            #
            #   En caso de que existe prtoeina siempre sera la primera
            #
            lst_molecules.insert(0, Molecule('Protein', g_target, self.cfg.name_target ))

        return lst_molecules

    def split_queries(self):
        """
            generamos un array con la prot y un dict con los ligandos con sus lineas sacadas del ficehro pdv

        """
        f = open(self.cfg.pdb)
        target = self.cfg.lst_molecules[0]
        for i in f:
            lig = False
            aux = re.sub(' +', ' ', i).strip().split(" ")

            if len(aux) > 0:	# Discarding solvents and ions
                if aux[0] == "ATOM" and aux[3] != self.cfg.name_solvent and aux[3] not in self.cfg.name_ions:
                    for mol in self.cfg.lst_molecules:
                        if aux[3] == mol.name:
                            mol.coords.append(i.strip())
                            lig = True
                            break
                    if not lig:
                        target.coords.append(i.strip())


    def generate_distribution_xvg(self, xvg_in, xvg_out):
        cmd = '{0} {1}{2}{3} -f {4}  -dist {5} -bw {6}'.format(
            self.cfg.gromacs,
            self.cfg.graph,
            "analyze",
            self.cfg.mpi,
            xvg_in,
            xvg_out,
            self.cfg.DISTRIBUTION_STEP
        )
        self.cfg.tools.execute.run(cmd)

    def generate_xvg(self, n_group_1, n_group_2, type_graph, out_xtc):
        cmd = 'echo  \"{0} {1}\" | {2} {3}{4}{5} -f {6} -s {7} -o {8} '.format(
            n_group_1,
            n_group_2,
            self.cfg.gromacs,
            self.cfg.graph,
            type_graph,
            self.cfg.mpi,
            self.cfg.xtc_md,
            self.cfg.tpr_min,
            out_xtc
        )
        self.cfg.tools.execute.run(cmd)

    def join_image(self, img_a, img_b, out):
        """
            join 2 images
        """
        comando="convert -size 2048x1536 "+img_a+" -thumbnail 800x600 "+img_a+".png"
        self.execute.run(comando)
        comando="convert -size 2048x1536 "+img_b+" -thumbnail 800x600 "+img_b+".png"
        self.execute.run(comando)
        comando="montage  -geometry +3+2  "+img_a+".png "+img_b+".png -tile 2x1 "+out
        self.execute.run(comando)
        comando="rm "+img_a+".png "+img_b+".png"
        self.execute.run(comando)

    #def get_molecule(self, mol):
    #    for i in self.cfg.lst_molecules:
    #        if i.name == mol:
    #            return i
    #    return None

    @staticmethod
    def check_directory(path):
        """
            Si existe, Borra el directorio y lo crea
            si no lo crea
        :param path: folder
        """
        if os.path.exists(path):
            shutil.rmtree(path)
            os.mkdir(path)
        else:
            os.mkdir(path)






