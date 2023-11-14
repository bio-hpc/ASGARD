#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Generates a fluctuation step by step graph

N_RESOLUTION_GRAPH = 100

class GraphStepFluctuation:

    def __init__(self, cfg):
        self.cfg = cfg
        if self.cfg.p_step_fluctuation:
            self.generate_graph()

    def generate_graph(self):
        """
        Generates xvg files with fluctuation data each x steps and create a graph after that
        """

        step = int(self.cfg.write_data) / N_RESOLUTION_GRAPH
        write_each = int(self.cfg.time_simulation) / int(self.cfg.write_data)
        step = write_each * step if write_each * step != 0 else 1

        lst_cmd = []
        lst_cmd.append("mkdir {} ".format(self.cfg.folder_fluctuation))
        cmd_gromacs = '"{} {}{}"'.format(self.cfg.gromacs, self.cfg.graph + "rmsf", self.cfg.mpi)
        lst_cmd.append('{0} {1}\\\n\t {2}\\\n\t {3}\\\n\t {4}\\\n\t {5}\\\n\t {6}\\\n\t {7}\\\n\t {8}'
                       .format(self.cfg.python_run, self.cfg.graph_step_fluctuation,
                               self.cfg.folder_fluctuation,
                               self.cfg.tpr_min,
                               self.cfg.xtc_md,
                               cmd_gromacs,
                               step,
                               N_RESOLUTION_GRAPH,
                               self.cfg.out_graph_fluctuation)
                       )
        lst_cmd.append("cd {}".format(self.cfg.path))
        # lst_cmd.append("rm -r {}".format(self.cfg.folder_fluctuation))
        self.cfg.template_job.execute_job(self.cfg.script_graph_fluctuation, lst_cmd)
