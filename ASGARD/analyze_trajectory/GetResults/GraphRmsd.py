class GraphRmsd():

    def __init__(self, cfg):
        self.cfg = cfg
        if self.cfg.p_graph_rmsd:
            self.rmsd_molecules_xvg()
            self.rmsd_molecules_png()

    def rmsd_molecules_xvg(self):
        #
        #   Generate all the xvg files with rmsd and rmsf data
        #
        for mol in self.cfg.lst_molecules:
            #   rmsd
            out_xvg = self.cfg.f_molecule_rmsd_xvg.format(self.cfg.prefix_results_xvg, mol.original_name )
            out_xvg_distribution = self.cfg.f_molecule_rmsd_distribution_xvg.format(self.cfg.prefix_results_xvg, mol.original_name )
            self.cfg.tools.generate_xvg(mol.group, mol.group, "rms", out_xvg)
            self.cfg.tools.generate_distribution_xvg(out_xvg, out_xvg_distribution)

            #  rmsdf
            xvg_f_out = self.cfg.f_molecule_rmsd_f_xvg.format(self.cfg.prefix_results_xvg, mol.original_name)
            png_f_out = self.cfg.f_molecule_rmsd_f_png.format(self.cfg.prefix_results_png, mol.original_name)
            self.cfg.tools.generate_xvg(mol.group, mol.group, "rmsf", xvg_f_out)
            cmd = '{} {} {} {} {}'.format(
                self.cfg.python_run,
                self.cfg.grap_rmsd,
                xvg_f_out,
                "\"RMSD Fluctuacion: " + mol.original_name + "\"",
                png_f_out
            )
            self.cfg.tools.execute.run(cmd)

    def generate_graph(self, lst, out, title):
        cmd = '{} {} {} {} {}'.format(
            self.cfg.python_run,
            self.cfg.grap_rmsd,
            " ".join(lst),
            title,
            out

        )
        self.cfg.tools.execute.run(cmd)


    def rmsd_molecules_png(self):
        #
        #   Generate a graph with all the rmsd of the molecules and calculate their distributions
        #
        lst_mols = [self.cfg.f_molecule_rmsd_xvg.format(self.cfg.prefix_results_xvg, i.original_name) for i
                    in self.cfg.lst_molecules]
        self.generate_graph(lst_mols,   self.cfg.f_rmsd_png.format(self.cfg.prefix_results_png), '\"Rmsd All Molecules\"')


        lst_mols = [self.cfg.f_molecule_rmsd_distribution_xvg.format(self.cfg.prefix_results_xvg, i.original_name) for i
                    in self.cfg.lst_molecules]

        self.generate_graph(lst_mols, self.cfg.f_rmsd_distribution_png.format(self.cfg.prefix_results_png), '\"Rmsd All Molecules\"')
        if len( self.cfg.lst_molecules) > 1:
            for mol in self.cfg.lst_molecules:
                lst_mols = [ self.cfg.f_molecule_rmsd_xvg.format(self.cfg.prefix_results_xvg, mol.original_name) ]
                self.generate_graph(lst_mols, self.cfg.f_molecule_rmsd_png.format(self.cfg.prefix_results_png,mol.original_name), '\"Rmsd ' + mol.original_name + ' \"')
                lst_mols = [self.cfg.f_molecule_rmsd_distribution_xvg.format(self.cfg.prefix_results_xvg, mol.original_name)]
                self.generate_graph(lst_mols,
                                    self.cfg.f_molecule_rmsd_distribution_png.format(self.cfg.prefix_results_png, mol.original_name),
                                    '\"Rmsd Distribution ' + mol.original_name + ' \"')

