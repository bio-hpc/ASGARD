import os
from .Tools import Tools
from .Profiles import Profiles
from .TemplateJob import TemplateJob
import subprocess
import re

class ConfigHolder(object):
    format_3 = '  {1:18} {2:<20}'
    format_2 = '  {:18}{:}'

    DISTRIBUTION_STEP = 0.001  # to generate the distribution plot (rmsd, distance)
    MAX_WARNINGS = 5
    DIST_MIN_RES = 20  # It starts with 20 and lowers until 0 to find the residues (the maximum number is 50 residues due to you only can include 64 groups in GROMACS)
    NUM_MAX_GROUPS = 52  # Max number of the groups considered (residues + gromacs groups)
    MAX_TABLE_MULTIMOL_DISTANCE = 7
    ENERGIES_DISCARD = 1  # Energy average from the simulation that is discarded for the proces_interactins_gromacs y graphs interaction gromacs scripts

    def setattr(self, k, v):
        setattr(self, k, v)

    def set_profile_cfg(self, profile):
        """
        Busca el perfil
        """
        self.profiles = Profiles(self)
        self.profiles.set_profile_cfg(profile)

    def check_simulation(self, out_xtc):
        out_check = os.path.splitext(out_xtc)[0] + ".chk"
        cmd = '{} check -f {} > {} 2>&1 '.format(
            self.gromacs,
            out_xtc,
            out_check
        )
        self.tools.execute.run(cmd)
        cmd = "cat {} |grep Step".format(out_check)
        step = self.tools.execute.run(cmd)
        step = re.sub(' +', ' ', step).strip().split(" ")[1]
        return int(step)

    def center_simulation(self, xtc):
        self.xtc_md = os.path.splitext(xtc)[0] + "_complex_no_end_center.xtc"
        tmp_file = "/" + self.xtc_md
        if not os.path.isfile(tmp_file.strip()):
            cmd = 'bash {0} {1} {2} {3} {4} {5}'.format(
                self.script_center_simulation,
                self.tpr_min,
                xtc,
                self.profile,
                self.xtc_md,
                self.gromacs
            )
            self.tools.execute.run(cmd)

            self.step_md = str(self.check_simulation(self.xtc_md))
            cmd = 'bash {0} {1} {2} {3} {4} {5}'.format(
                self.script_crete_pdb,
                self.xtc_md,
                self.gro_md,
                self.folder_molec + self.sufijo + "_no_end.pdb",
                self.step_md,
                self.gromacs
            )
            self.tools.execute.run(cmd)

    def checkFile(self, file):
        tmp_file = "/" + file
        if os.path.isfile(tmp_file.strip()):
            return file
        else:
            print("\n")
            print("ERROR: you must enter the absolute path of the file, this usually happens"
                  " because you do not pass the file with the absolute path: ")
            print("No exsis: " + tmp_file)
            exit("\n")

    def __init__(self, prefix_molec, profile, gromacs):
        #
        #   Input-output folders
        #
        self.gromacs = gromacs + " "

        tmp = os.path.dirname(prefix_molec)
        self.folder = tmp[:tmp.rfind("/")] + "/"
        self.sufijo = os.path.basename(prefix_molec).strip()

        self.folder_molec = os.path.join(self.folder, 'molecules/')
        self.folder_grids = os.path.join(self.folder, 'grids/')
        self.folder_out = os.path.join(self.folder, 'out/')
        folder_template = "jobs/"

        self.log = self.folder + "resu_dm.txt"
        if os.path.isfile(self.log):
            os.remove(self.log)

        self.profile = profile
        self.tools = Tools(self)
        self.set_profile_cfg(profile)
        cmd = "ls " + self.folder_molec + self.sufijo + "_npt_*.gro"
        self.gro_md = self.tools.execute.run(cmd).strip()

        self.template_job = TemplateJob(self)
        #
        #
        #   gromacs different options
        #
        self.cmd_check = self.gromacs + " check "
        self.mpi = ""
        self.graph = ""
        self.gpu = ""
        self.threds = "-ntomp " + self.template_job.cores
        #
        #

        # Plot generation
        self.path = os.getcwd() + "/"

        self.script_center_simulation = os.path.join(self.path, "ASGARD/cluster_nodes/execute_scripts/scriptGR/center_simulation.sh")
        self.script_crete_pdb = os.path.join(self.path, "ASGARD/cluster_nodes/execute_scripts/scriptGR/create_pdb.sh")

        folder_graphs = self.path + "ASGARD/analyze_trajectory/Graphs/"
        self.standar_graph_xvg = folder_graphs + "standar_graph_xvg.py"  # standard graphs

        self.graph_gyrate_helicity = folder_graphs + "graph_gyrate_helicity.py"  #
        self.graph_step_fluctuation = folder_graphs + "graph_step_fluctuation.py"  # rmsd f for simulacion step

        self.process_interactions_gromacs = folder_graphs + "process_interactions_gromacs.py"  # protein-ligand energy interaction
        self.graph_interactions_gromacs = folder_graphs + "graph_interactions_gromacs.py"  # energies graph
        self.graph_sasa = folder_graphs + "graph_sasa.py"  # SASA graph
        self.dssp = folder_graphs + "dssp-2.0.4-linux-amd64"  # DSSP software
        self.graph_dssp = folder_graphs + "graph_dssp.py"  # DSSP graph

        self.logo_bio_hpc = self.path + "ASGARD/analyze_trajectory/extra/logo_biohpc.png"
        self.md_diagram = self.path + "ASGARD/analyze_trajectory/extra/md_diagram.png"

        try:
            cmd = 'cat {}/config.cfg |grep g_mmpbsa'.format(os.path.join(self.path, "ASGARD"))
            a_h = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True).decode('utf-8')
            self.g_mmpbsa = a_h.split(":")[1].strip()  # "ASGARD/analyze_results/Simulation_gromacs/analyze_trajectory/extra/g_mmpbsa"
        except:
            print("Error en el fichero de configuracion \"g_mmpbsa\"")
            exit()

        self.json_text_tex = self.path + "ASGARD/analyze_trajectory/extra/text_document_tex.json"
        self.graph_hbonds = folder_graphs + "graph_hbonds.py"
        self.graph_mmpbsa = folder_graphs + "graph_mmpbsa.py"
        self.grap_rmsd = folder_graphs + "graph_rmsd.py"  # rmsd ligand and distance graphs
        self.mdrum = "mdrum"

        #
        # Directorios utilizados para el analisis
        #
        self.index = self.folder_grids + self.sufijo + "_index.ndx"

        # self.tpr_min = self.checkFile(self.folder_molec+self.sufijo+"_min.tpr")
        self.tpr_min = self.checkFile(self.folder_molec + self.sufijo + "_md.tpr")
        self.tpr_pre = self.checkFile(self.folder_molec + self.sufijo + "_pre_md.tpr")

        self.prefix_molec = self.folder_molec + self.sufijo
        self.grids = self.folder_grids + self.sufijo
        self.prefix_templates = self.folder + folder_template + self.sufijo
        self.out_aux = self.folder_out + self.sufijo
        self.top = self.checkFile(self.folder_molec + self.sufijo + ".top")

        #
        #   Opciones del target y queries
        #
        self.name_target = os.path.basename(prefix_molec).split("_")[2]
        self.python_run = "python"
        self.lst_molecules = self.tools.get_groups_target_queries()

        self.g_energy = self.gromacs + " energy "
        self.grompp = '{} {} '.format(self.gromacs, 'grompp')
        self.mdrun = '{} {} '.format(self.gromacs, 'mdrun')

        if os.path.isfile(self.folder_molec + self.sufijo + "_md_center.xtc"):
            self.xtc_md = self.checkFile(self.folder_molec + self.sufijo + "_md_center.xtc")
            self.gro_md = self.checkFile(self.folder_molec + self.sufijo + "_md.gro")
            self.pdb = self.checkFile(self.folder_molec + self.sufijo + "_md.pdb")
        else:
            self.checkFile(self.folder_molec + self.sufijo + "_md.xtc")
            cmd = "ls " + self.folder_molec + self.sufijo + "_npt_*.gro"
            self.gro_md = self.tools.execute.run(cmd).strip()

            self.center_simulation(self.folder_molec + self.sufijo + "_md.xtc")
            self.pdb = self.checkFile(self.folder_molec + self.sufijo + "_no_end.pdb")

        ##self.trr_md = self.checkFile(self.folder_molec+self.sufijo+"_md.trr")
        self.edr_md = self.checkFile(self.folder_molec + self.sufijo + "_md.edr")

        cmd = "ls {}Resume_VS_GR_*".format(self.folder)
        self.resume_file = self.checkFile(self.tools.execute.run(cmd).strip())
        cmd = "ls " + self.folder_molec + self.sufijo + "_npt_*.gro"
        self.last_equilibration = self.tools.execute.run(cmd).strip()
        cmd = "ls {}*.mdp".format(self.folder_grids)
        self.config_files = self.tools.execute.run(cmd)
        for i in self.config_files.split("\n"):
            if i != "":
                aux = os.path.splitext(i)[0]
                name = aux[aux.rfind('_') + 1:]
                if name == "md":
                    self.config_file_md = i
        self.folder_results = os.path.join(self.folder, 'results')
        self.folder_results_xvg = os.path.join(self.folder_results, 'xvg/')
        self.folder_results_png = os.path.join(self.folder_results, 'png/')
        self.prefix_results = os.path.join(self.folder_results, self.sufijo)
        self.prefix_results_xvg = os.path.join(self.folder_results_xvg, self.sufijo)
        self.prefix_results_png = os.path.join(self.folder_results_png, self.sufijo)

        if not os.path.isdir(self.folder_results):
            os.mkdir(self.folder_results)
        if not os.path.isdir(self.folder_results_xvg):
            os.mkdir(self.folder_results_xvg)
        if not os.path.isdir(self.folder_results_png):
            os.mkdir(self.folder_results_png)

        #
        #  File outputs and inputs
        #

        self.folder_fluctuation = self.prefix_results_xvg + "_heatmap/"
        self.out_graph_fluctuation = '{}_step_fuctuation_protein.png'.format(self.prefix_results_png)
        self.script_graph_fluctuation = '{}_step_fuctuation.sh'.format(self.prefix_templates)
        self.script_rerun_md = '{}_rerun_md.sh'.format(self.prefix_templates)
        self.out_xvg_ddsp = '{}_dssp.xvg'.format(self.prefix_results_xvg)
        self.out_xmp_ddsp = '{}_dssp.xmp'.format(self.prefix_results_xvg)
        self.out_eps_ddsp = '{}_dssp.eps'.format(self.prefix_results_xvg)
        self.out_png_ddsp = '{}_dssp.png'.format(self.prefix_results_png)
        self.out_png_ddsp_ss = '{}_dssp_ss.png'.format(self.prefix_results_png)
        self.script_ddsp = '{}_dssp.sh'.format(self.prefix_templates)
        self.out_xvg_odg_sasa = '{}_odg_sasa.xvg'.format(self.prefix_results_png)
        self.out_xvg_o_sasa = '{}_o_sasa.xvg'.format(self.prefix_results_png)
        self.out_png_sasa = '{}_sasa.png'.format(self.prefix_results_png)
        self.script_sasa = '{}_sasa.sh'.format(self.prefix_templates)
        self.out_png_stabilization_1 = '{}_estabilization_1.png'.format(self.prefix_results_png)
        self.out_png_stabilization_2 = '{}_estabilization_2.png'.format(self.prefix_results_png)
        self.folder_helicity = self.prefix_results_xvg + "_helicity/"
        self.out_xvg_helicity = self.prefix_results_xvg + "_helicity.xvg"
        self.out_png_helicity = self.prefix_results_png + "_helicity.png"
        self.out_xvg_gyrate = self.prefix_results_xvg + "_girate.xvg"
        self.out_png_gyrate = self.prefix_results_png + "_girate.png"
        self.table_multimolecule = self.folder_results + "/" + self.sufijo + "_table_multi_molecule.tex"
        self.document_tex = self.folder_results + "/" + self.sufijo + "_documnet.tex"
        self.table_stabilization = self.prefix_results + "_table_stabilization.tex"

        self.format_file_name_hbonds_xvg = '{}_{}_{}_hbond.xvg'
        self.format_file_name_hbonds_png = '{}_{}_{}_hbond.png'
        self.format_file_name_g_mmpbsa_xvg = '{}_{}_{}_mmpbsa.xvg'
        self.format_file_name_g_mmpbsa_png = '{}_{}_{}_mmpbsa.png'
        self.format_file_name_hbonds_mmpbsa = '{}_{}_{}_hbond_mmpbsa.png'
        self.format_folder_interactions = '{}_{}_{}_interactions'  # gromacs energy interactions
        self.format_folder_dssp = '{}_dssp'  # dssp
        self.format_script_interactions = '{}_{}_{}_interactions.sh'

        self.f_molecule_rmsd_xvg = '{}_{}_rmsd.xvg'
        self.f_molecule_rmsd_distribution_xvg = '{}_{}_rmsd_distribution.xvg'
        self.f_molecule_rmsd_png = '{}_{}_rmsd.png'
        self.f_molecule_rmsd_distribution_png = '{}_{}_rmsd_distribution.png'
        self.f_molecule_rmsd_f_xvg = '{}_{}_rmsd_fluct.xvg'
        self.f_molecule_rmsd_f_png = '{}_{}_rmsd_fluct.png'
        self.f_rmsd_png = '{}_rmsd.png'
        self.f_rmsd_distribution_png = '{}_rmsd_distribution.png'
        self.f_molecule_distance_xvg = '{}_{}_{}_distance.xvg'
        self.f_molecule_distance_distribution_xvg = '{}_{}_{}_distance_distribution.xvg'
        self.f_molecule_distance_png = '{}_{}_{}_distance.png'
        self.f_molecule_distance_distribution_png = '{}_{}_{}_distance_distribution.png'
        self.f_distance_png = '{}_distance.png'
        self.f_distance_distribution_png = '{}_distance_distribution.png'

        #
        # Search parameters
        #
        f = open(self.resume_file)
        for i in f:

            if i.find("-- Command:") != -1:
                aux = i.split(" ")
                for j in range(len(aux)):
                    if aux[j] == "-step_md":
                        self.step_md = aux[j + 1].strip()
                    if aux[j] == "-step_npt":
                        self.step_npt = aux[j + 1]
                    if aux[j] == "-step_nvt":
                        self.step_nvt = aux[j + 1]
                    if aux[j] == "-step_min":
                        self.step_min = aux[j + 1]
                    if aux[j] == "-force_field":
                        self.force_field = aux[j + 1]
                    if aux[j] == "-solvent":
                        self.solvent = aux[j + 1]
                    if aux[j] == "-temp":
                        self.temp = aux[j + 1]
                    if aux[j] == "-bt":
                        self.type_grid = aux[j + 1]
                    if aux[j] == "typeGrid=":
                        self.padding = aux[j + 1]
                    if aux[j] == "-write_data":
                        self.write_data = aux[j + 1]
                    if aux[j] == "-padding_grid":
                        self.padding_grid = aux[j + 1]
        f.close()

        self.ph = "9"
        self.integration_step = 0.002
        self.time_simulation = int(int(self.step_md) * self.integration_step)
        self.ensemble = "NPT2"
        self.coulomb_type = "PME"
        self.t_coupl = "V-rescale"
        self.p_coupl = "Parrinello-Rahman"
        self.name_solvent = "SOL"
        self.name_ions = ["NA", "CL"]
        self.lst_jobs = []
        self.tools.split_queries()


