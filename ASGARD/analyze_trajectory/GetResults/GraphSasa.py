#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Generates the rmsd graphs for protein and ligand (if it exists)
#
#


class GraphSasa(object):
    def __init__(self, cfg):
        self.cfg = cfg
        if self.cfg.p_graph_sasa:
            lst_cmd = []
            lst_cmd.append('echo {0} | {1} {2}{3}{4} -f {5} -s {6} -o {7}'.format(
                self.cfg.lst_molecules[0].group,
                self.cfg.gromacs,
                self.cfg.graph,
                "sasa",
                self.cfg.mpi,
                self.cfg.xtc_md,
                self.cfg.tpr_min,
                self.cfg.out_xvg_o_sasa
                )
            )
            lst_cmd.append('echo {0} | {1} {2}{3}{4} -f {5} -s {6} -odg {7}'.format(

                self.cfg.lst_molecules[0].group,
                self.cfg.gromacs,
                self.cfg.graph,
                "sasa",
                self.cfg.mpi,
                self.cfg.xtc_md,
                self.cfg.tpr_min,
                self.cfg.out_xvg_odg_sasa
            )
            )

            lst_cmd.append('{} {} {} {} {}'.format(self.cfg.python_run, self.cfg.graph_sasa,self.cfg.out_xvg_o_sasa, self.cfg.out_xvg_odg_sasa, self.cfg.out_png_sasa))
            self.cfg.template_job.execute_job(self.cfg.script_sasa, lst_cmd)








