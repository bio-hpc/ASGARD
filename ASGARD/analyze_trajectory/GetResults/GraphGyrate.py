class GraphGyrate():
    def __init__(self, cfg):
        self.cfg = cfg
        if self.cfg.p_graph_gyrate:
            cmd = 'echo {0} | {1} {2}{3}{4} -f {5} -s {6} -o {7}'.format(
                self.cfg.lst_molecules[0].group,
                self.cfg.gromacs,
                self.cfg.graph,
                "gyrate",
                self.cfg.mpi,
                self.cfg.xtc_md,
                self.cfg.tpr_min,
                self.cfg.out_xvg_gyrate
            )
            self.cfg.tools.execute.run(cmd)
            cmd = '{} {} {} {}'.format(self.cfg.python_run, self.cfg.graph_gyrate_helicity, self.cfg.out_xvg_gyrate, self.cfg.out_png_gyrate)
            self.cfg.tools.execute.run(cmd)
