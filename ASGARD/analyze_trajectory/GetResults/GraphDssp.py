class GraphDssp():
    def __init__(self, cfg):
        self.cfg = cfg
        if self.cfg.p_graph_dssp:

            folder_xvg_out = self.cfg.format_folder_dssp.format(self.cfg.prefix_results_xvg )
            self.cfg.tools.check_directory(folder_xvg_out)
            lst_cmd = []
            lst_cmd.append('cd {} '.format(folder_xvg_out))
            lst_cmd.append('export DSSP={}'.format(self.cfg.dssp))
            lst_cmd.append('echo {0} | {1} do_dssp{2} -s {3} -f {4} -sc {5} -o {6}'.format(
                self.cfg.lst_molecules[0].group,
                self.cfg.gromacs,
                self.cfg.mpi,
                self.cfg.tpr_min,
                self.cfg.xtc_md,
                self.cfg.out_xvg_ddsp,
                self.cfg.out_xmp_ddsp
                )
            )
            lst_cmd.append('{0} xpm2ps{1} -f {2} -o {3} -by 50 -skip 5 -size 1024'.format(
                self.cfg.gromacs,
                self.cfg.mpi,
                self.cfg.out_xmp_ddsp,
                self.cfg.out_eps_ddsp
                )
            )
            lst_cmd.append('convert {} {}'.format(self.cfg.out_eps_ddsp, self.cfg.out_png_ddsp_ss))
            lst_cmd.append('cd {}'.format(self.cfg.path))
            lst_cmd.append('{} {} {} {}'.format(self.cfg.python_run, self.cfg.graph_dssp, self.cfg.out_xvg_ddsp, self.cfg.out_png_ddsp))
            self.cfg.template_job.execute_job(self.cfg.script_ddsp, lst_cmd)
