import os
from os.path import join
import subprocess


class ConfigHolder():
    def __init__(self, target, dir_queries, profile, gmx):
        #
        # Paths and scripts
        #
        self.path = os.getcwd() + "/"
        self.path_external_sw = join(self.path, "ASGARD/external_sw/")
        self.path_external_sw_gr = join(self.path_external_sw, 'gromacs')
        self.check_protein = join("ASGARD/extra_ASGARD/used_by_ASGARD/check_target.py")
        self.mod_diedrals_phosphate = join(self.path_external_sw_gr, 'topology', 'mod_itp_phosphate.py')
        
        #if profile == "DNA_QUERY":
        #    #   DNA
        #    self.force_field = "amber99bsc1.ff/"
        #    self.option_force_field = 1
        #else:
        #    #   Protein
        self.force_field = "amber99sb.ff/" #demomento solo tenemos 1 campo de fuerza
        self.option_force_field = 1
        
        self.path_force_field = join(self.path_external_sw_gr, 'force_field')
        self.path_forece_files_choice = join(self.path_external_sw_gr, 'force_field',self.force_field)
        self.acpype = join( self.path_external_sw, 'amber14', 'bin', 'acpype', 'acpype.py')
        
        # Ficheros topologials
        #
        self.target = target
        self.dir_queries = dir_queries
        self.dir_target = os.path.dirname(target)+"/"
        self.itp_target_chain = []  # cuando hay varas cadenas en la proteina se almacenen sus itps
        self.file_name_target, self.file_ext_target = os.path.splitext(target)
        self.name_target = os.path.basename(self.file_name_target) if self.file_name_target != "" else "queries"
        self.target_file_gro = join(self.path, self.file_name_target + '.gro') if self.file_name_target != "" else None
        self.target_file_top = join(self.path, self.file_name_target + '.top') if self.file_name_target != "" else None
        self.target_file_itp = join(self.path, self.file_name_target + '_porse.itp') if self.file_name_target != "" else None
        self.out_complex = os.path.basename(os.path.splitext(self.file_name_target)[0]) if self.file_name_target != "" else profile.lower()
        #
        # AmberHome para convertir el ligando
        #
        try:
            cmd  = 'cat {}/config.cfg |grep g_amber_home'.format(join(self.path, "ASGARD"))
            a_h  = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
            if a_h != "":
                os.environ["AMBERHOME"] = a_h.split(":")[1].strip()
                os.environ["PATH"] = os.environ["AMBERHOME"] + ":" + os.environ["PATH"]
                os.environ["PATH"] = os.environ["AMBERHOME"] + "bin/:" + os.environ["PATH"]
                """
                print os.environ["AMBERHOME"]
                print os.environ["PATH"]
                print os.environ["PATH"]
                """
        except :
            print ("Error en el fichero de configuracion")
            exit()
        #
        #   Gromcas  enviroment varible
        #
        os.environ["GMX_MAXBACKUP"] = "-1"  # se eliminan backups de gromacs
        #
        #   Parameters
        #
        self.prefix_sumulation = '_complex'
        self.num_queries = 1
        self.python_run = 'python'

        self.gromacs_run = 'gmx'
        self.solvent = 'tip3p'
        self.ignh = '-ignh'
        self.query_ext = '.mol2'

        #
        #   Formats
        #
        self.format_gro ='{:>5}{:5}{:>5}{:>5}{:>8.3f}{:>8.3f}{:>8.3f}'
        # "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
        self.format_itp = '{0:>6}{1:>6}{2:>6}{3:>6}{4:>6}'
        self.format_out_2 = '{:30}{:50}'
        self.format_out_3 = '{:30}{:50}{:100}'
        #
        #   Profiles
        #   1 target 1 query por defecto
        #   1 target varias queries
        #   1 target
        #   1 o varias queries
        self.profile = profile.upper()
        self.lst_profiles = ['TARGET_QUERY',  'TARGET_QUERIES', 'TARGET', 'QUERIES', 'TARGET_ONE_QUERY', 'DNA_QUERY', 'BIPHSIC_SYSTEMS']



