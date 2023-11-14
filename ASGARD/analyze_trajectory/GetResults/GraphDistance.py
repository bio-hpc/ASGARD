class GraphDistance:

    def __init__(self, cfg):
        self.cfg = cfg
        if self.cfg.p_graph_distance or self.cfg.p_table_multimolecule:
            self.mol_target = self.cfg.lst_molecules[0]
            self.distance_molecules_xvg()
            self.distance_molecules_png()

    def distance_molecules_xvg(self):
        #
        #   The first molecule is always taken as reference
        #   If there is only a molecule, the graph is not generated
        #

        for i in range(1, len(self.cfg.lst_molecules)):
            mol_query = self.cfg.lst_molecules[i]
            out_xvg = self.cfg.f_molecule_distance_xvg.format(self.cfg.prefix_results_xvg, self.mol_target.original_name,
                                                             mol_query.original_name)
            out_xvg_distribution = self.cfg.f_molecule_distance_distribution_xvg.format(self.cfg.prefix_results_xvg,
                                                                                        self.mol_target.original_name,
                                                                                        mol_query.original_name)
            cmd = '{0} {1}{2}{3} -f {4} -s {5} -oall {6} -select \'com of group {7} plus com of group {8}\''.format(
                self.cfg.gromacs,
                self.cfg.graph,
                "distance",
                self.cfg.mpi,
                self.cfg.xtc_md,
                self.cfg.tpr_min,
                out_xvg,
                self.mol_target.group,
                mol_query.group
            )
            self.cfg.tools.execute.run(cmd)
            self.cfg.tools.generate_distribution_xvg(out_xvg, out_xvg_distribution)

    def generate_graph(self, lst, out, lst_names):
        cmd = '{} {} {} {} {}'.format(
            self.cfg.python_run,
            self.cfg.grap_rmsd,
            " ".join(lst),
            '\"Distance ' + ' '.join(lst_names) + ' \"',
            out
        )
        self.cfg.tools.execute.run(cmd)

    def distance_molecules_png(self):
        #
        #   A graph with all the distances of the molecules is created and the distribution is done
        #
        lst_mols = [self.cfg.f_molecule_distance_xvg.format(self.cfg.prefix_results_xvg, self.mol_target.original_name,
                                                           i.original_name) for i
                    in self.cfg.lst_molecules]
        lst_names = [i.original_name for i in self.cfg.lst_molecules]
        lst_mols.pop(0)
        self.generate_graph(lst_mols, self.cfg.f_distance_png.format(self.cfg.prefix_results_png), lst_names)

        lst_mols = [self.cfg.f_molecule_distance_distribution_xvg.format(self.cfg.prefix_results_xvg,
                                                                         self.mol_target.original_name, i.original_name)
                    for i in self.cfg.lst_molecules]
        lst_names = [i.original_name for i in self.cfg.lst_molecules]
        lst_mols.pop(0)
        self.generate_graph(lst_mols, self.cfg.f_distance_distribution_png.format(self.cfg.prefix_results_png),
                            lst_names)
        if len(self.cfg.lst_molecules) > 2:
            for i in range(1, len(self.cfg.lst_molecules)):
                mol_query = self.cfg.lst_molecules[i]
                in_xvg = [self.cfg.f_molecule_distance_xvg.format(self.cfg.prefix_results_xvg,
                                                                  self.mol_target.original_name,
                                                                  mol_query.original_name)]
                lst_names = [self.mol_target.original_name, mol_query.original_name]
                self.generate_graph(in_xvg,
                                    self.cfg.f_molecule_distance_png.format(self.cfg.prefix_results_png,
                                                                              self.mol_target.original_name,
                                                                              mol_query.original_name), lst_names)

                in_xvg = [self.cfg.f_molecule_distance_distribution_xvg.format(self.cfg.prefix_results_xvg,
                                                                               self.mol_target.original_name,
                                                                               mol_query.original_name)]
                self.generate_graph(in_xvg,
                                    self.cfg.f_molecule_distance_distribution_png.format(self.cfg.prefix_results_png,
                                                                                          self.mol_target.original_name,
                                                                                          mol_query.original_name),
                                    lst_names)
