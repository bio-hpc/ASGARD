class GraphRmsd:

    def __init__(self, cfg):
        self.cfg = cfg
        if self.cfg.p_graph_rmsd:
            self.rmsd_molecules_xvg()
            self.rmsd_molecules_png()

    def rmsd_molecules_xvg(self):
        for mol in self.cfg.lst_molecules:
            # rmsd
            out_xvg = self.cfg.f_molecule_rmsd_xvg.format(self.cfg.prefix_results_xvg, mol.original_name)
            out_xvg_distribution = self.cfg.f_molecule_rmsd_distribution_xvg.format(self.cfg.prefix_results_xvg,
                                                                                    mol.original_name)
            self.cfg.tools.generate_xvg(mol.group, mol.group, "rms", out_xvg)
            self.cfg.tools.generate_distribution_xvg(out_xvg, out_xvg_distribution)

            # rmsdf
            xvg_f_out = self.cfg.f_molecule_rmsd_f_xvg.format(self.cfg.prefix_results_xvg, mol.original_name)
            png_f_out = self.cfg.f_molecule_rmsd_f_png.format(self.cfg.prefix_results_png, mol.original_name)
            self.cfg.tools.generate_xvg(mol.group, mol.group, "rmsf", xvg_f_out)
            cmd = f'{self.cfg.python_run} {self.cfg.grap_rmsd} {xvg_f_out} "RMSD Fluctuacion: {mol.original_name}" {png_f_out}'
            self.cfg.tools.execute.run(cmd)

    def generate_graph(self, lst, out, title):
        cmd = f'{self.cfg.python_run} {self.cfg.grap_rmsd} {" ".join(lst)} {title} {out}'
        self.cfg.tools.execute.run(cmd)

    def rmsd_molecules_png(self):
        if len(self.cfg.lst_molecules) == 1:
          for i in self.cfg.lst_molecules:
            original_name_target=i.original_name
          lst_mols = [self.cfg.f_molecule_rmsd_xvg.format(self.cfg.prefix_results_xvg, i.original_name) for i in self.cfg.lst_molecules]
          self.generate_graph(lst_mols, self.cfg.f_molecule_rmsd_png.format(self.cfg.prefix_results_png,original_name_target), '"RMSD All Molecules"')
          lst_mols = [self.cfg.f_molecule_rmsd_distribution_xvg.format(self.cfg.prefix_results_xvg, i.original_name) for i in self.cfg.lst_molecules]
          self.generate_graph(lst_mols, self.cfg.f_molecule_rmsd_distribution_png.format(self.cfg.prefix_results_png,original_name_target), '"RMSD All Molecules"')
          
        lst_mols = [self.cfg.f_molecule_rmsd_xvg.format(self.cfg.prefix_results_xvg, i.original_name) for i in self.cfg.lst_molecules]
        self.generate_graph(lst_mols, self.cfg.f_rmsd_png.format(self.cfg.prefix_results_png), '"RMSD All Molecules"')
        lst_mols = [self.cfg.f_molecule_rmsd_distribution_xvg.format(self.cfg.prefix_results_xvg, i.original_name) for i in self.cfg.lst_molecules]
        self.generate_graph(lst_mols, self.cfg.f_rmsd_distribution_png.format(self.cfg.prefix_results_png), '"RMSD All Molecules"')

        if len(self.cfg.lst_molecules) > 1:
            for mol in self.cfg.lst_molecules:
                lst_mols = [self.cfg.f_molecule_rmsd_xvg.format(self.cfg.prefix_results_xvg, mol.original_name)]
                self.generate_graph(lst_mols, self.cfg.f_molecule_rmsd_png.format(self.cfg.prefix_results_png, mol.original_name), f'"RMSD {mol.original_name}"')
                lst_mols = [self.cfg.f_molecule_rmsd_distribution_xvg.format(self.cfg.prefix_results_xvg, mol.original_name)]
                self.generate_graph(lst_mols, self.cfg.f_molecule_rmsd_distribution_png.format(self.cfg.prefix_results_png, mol.original_name), f'"RMSD Distribution {mol.original_name}"')