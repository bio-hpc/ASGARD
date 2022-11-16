import os
class GraphHelix():
    def __init__(self, cfg):
        self.cfg = cfg
        if self.cfg.p_graph_helicity:
            self.graph_helix()


    def graph_helix(self):
        #
        #	Genera grafica de helices (RAPIDA)
        #
        self.cfg.tools.check_directory (self.cfg.folder_helicity)
        os.chdir(self.cfg.folder_helicity)
        cmd = 'echo {0} | {1} {2}{3}{4} -f {5} -s {6} -n {7}'.format(
            self.cfg.lst_molecules[0].group,
            self.cfg.gromacs,
            self.cfg.graph,
            "helix",
            self.cfg.mpi,
            self.cfg.xtc_md,
            self.cfg.tpr_min,
            self.cfg.index_target_query
        )
        print cmd
        self.cfg.tools.execute.run(cmd)
        os.rename("helicity.xvg", self.cfg.out_xvg_helicity)
        os.chdir(self.cfg.path)
        cmd = '{} {} {}'.format(self.cfg.python_run, self.cfg.graph_gyrate_helicity, self.cfg.out_png_helicity )
        self.cfg.tooles.execute.run(cmd)

        #print (cmd)
        #################################################################################
        #	Grafica auxiliar solo para TMI Borrar Para otrsa dinamicas
        #
        #filename, file_extension = os.path.splitext(self.cfg.graoGR)
        #comando = self.cfg.python + " " + filename + "Aux.py " + self.cfg.outTxt + "_Helix.xvg 110 140"
        #self.execute.run(comando)
        #comando = self.cfg.python + " " + filename + "Aux.py " + self.cfg.outTxt + "_Helix.xvg 385 415"
        #self.execute.run(comando)
        #self.cfg.tools.joinImage(self.cfg.outTxt + "_Helix" + "_110_140_TMI.png",
        #                         self.cfg.outTxt + "_Helix" + "_385_415_TMI.png",
        #                         self.cfg.outTxt + "_Helix_TMI.png")
        #shutil.rmtree(helixfolder)
        ###################################################################################
#comando = self.cfg.python + " " + self.cfg.graph_gyrate_helicity + " " + self.cfg.outTxt + "_Helix.xvg"
#self.execute.run(comando)
#helixfolder = self.cfg.outTxt + "Helix/"
#self.checkDirectorio(helixfolder)
#comando = "echo " + self.cfg.protein[
#    self.cfg.nomProtein] + " | " + self.cfg.gromacs + " " + self.cfg.graph + "helix" + self.cfg.mpi + " -f " + self.cfg.xtcSimu + " -s " + self.cfg.tprMin + " -n " + self.cfg.grids + "_index_prot_lig_sol.ndx"
#self.execute.run(comando)

# shutil.rmtree(helixfolder)  #decomentar cuando no sea TMI