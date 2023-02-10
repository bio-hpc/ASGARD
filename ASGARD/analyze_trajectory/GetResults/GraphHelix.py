import os
class GraphHelix():
    def __init__(self, cfg):
        self.cfg = cfg
        if self.cfg.p_graph_helicity:
            self.graph_helix()


    def graph_helix(self):
        #
        #	Generate helix graph
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
