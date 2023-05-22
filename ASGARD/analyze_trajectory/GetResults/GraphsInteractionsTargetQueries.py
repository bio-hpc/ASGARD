#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import fileinput


class GraphsInteractionsTargetQueries(object):
    """
        Generates the coordinates xvg files for the outData graphs
    """
    def __init__(self, cfg):
        self.cfg = cfg
        self.input_file_tmp = self.cfg.grids + "_tmp.tmp.ndx"
        self.index_all_residues = self.cfg.grids + "_all_residues.ndx"
        self.rerun_md = self.cfg.prefix_molec + "_rerun_md"
        self.rerun_md_log = self.cfg.grids + "_rerun_md.log"
        self.generate_graph()


    def generate_graph(self):
        """
            Creates an index with all the residues of the protein (without "Protein_" and rerun with the rest of the molecules
        """
        
        if self.cfg.p_interactions_gromacs:
            self.generate_index()
            self.generate_rerun()
        
        dependencies = 0  
        if not self.cfg.p_sequential:
            dependencies = self.cfg.lstJobs[len(self.cfg.lstJobs)-1]

        mol_target = self.cfg.lst_molecules[0] #the protein is always 0
        for i_mol in range(1, len(self.cfg.lst_molecules)):
            mol_query = self.cfg.lst_molecules[i_mol]
            print(mol_query)
            print(mol_target)
            self.generate_graph_hbonds(mol_target, mol_query)
            self.generate_graph_van_elect(mol_target, mol_query)

            if self.cfg.p_interactions_gromacs:
               
                folder_xvg_out = self.cfg.format_folder_interactions.format(self.cfg.prefix_results_xvg, mol_target.original_name, mol_query.original_name)
                self.cfg.tools.check_directory(folder_xvg_out)

                cad_residues = self.cfg.tools.get_about_residues(mol_target, mol_query)

                mpd_file = '{}_{}_{}.mdp'.format(os.path.splitext(self.cfg.config_file_md)[0], mol_target.original_name, mol_query.original_name)
                self.cfg.tools.cp_file(self.cfg.config_file_md, mpd_file)

                # add the group line which check the energies in rerun step
                f_mdp = open(mpd_file, "a") 
                f_mdp.write("energygrps\t     =  "+cad_residues +"\n")
                f_mdp.close()        
                prefix_query = '{}_{}_{}_rerun'.format(self.cfg.prefix_results, mol_target.original_name, mol_query.original_name)
                prefix_query_out = '{}_{}_{}_rerun'.format(self.cfg.prefix_results_png, mol_target.original_name, mol_query.original_name)
                lst_cmd = self.rerun_target_query(mol_target, mol_query, mpd_file, prefix_query)

                lst_cmd.append(
                    '{0} {1}\\\n\t {2}\\\n\t {3}\\\n\t {4}\\\n\t {5}\\\n\t {6}\\\n\t {7}\\\n\t {8}\\\n\t'
                    ' {9}\\\n\t {10}\\\n\t {11}\\\n\t {12}\\\n\t {13}\\\n\t {14}\\\n\t {15}\\\n\t'
                    ' {16}\\\n\t {17}\\\n\t {18}'.format(
                        self.cfg.python_run,
                        self.cfg.process_interactions_gromacs,

                        self.cfg.gromacs,
                        "\""+self.cfg.graph+"\"",
                        "\""+self.cfg.mpi+"\"",
                        self.rerun_md,
                        prefix_query,
                        self.cfg.tpr_min,
                        folder_xvg_out,
                        self.cfg.name_solvent,
                        self.cfg.graph_interactions_gromacs,
                        prefix_query_out,
                        mol_query.name,
                        mol_query.original_name,
                        mol_target.name,
                        mol_target.original_name,
                        "\""+self.cfg.g_energy+"\"",
                        self.cfg.pdb,
                        self.cfg.ENERGIES_DISCARD
                    )
                )
                #lst_cmd.append('rm {}_{}_{}_rerun*'.format(self.cfg.prefix_results, mol_target.original_name, mol_query.original_name))
                self.cfg.template_job.execute_job_with_dependency(
                    self.cfg.format_script_interactions.format(self.cfg.prefix_templates, mol_target.original_name, mol_query.original_name),
                    lst_cmd, dependencies)

    def generate_rerun(self):
        self.cfg.template_job.execute_job(self.cfg.script_rerun_md, self.def_rerun_md())

    def generate_index(self):
        aux=""
        for i in self.cfg.lst_molecules[1:]:
            aux += aux + " " + i.name

        f = open(self.input_file_tmp,'w')
        f.write('splitres 1\n')
        f.write('"Water_and_ions'+aux.replace(" ","_")+'"''\n')
        f.write("q\n")
        f.close()

        cmd = '{} {} -f {} -o {} -n {} < {}'.format(self.cfg.gromacs, self.cfg.graph + "make_ndx" + self.cfg.mpi,
                                              self.cfg.gro_md, self.index_all_residues,  self.cfg.index, self.input_file_tmp)
        self.cfg.tools.execute.run(cmd)
        self.cfg.tools.execute.run('rm {}'.format(self.input_file_tmp))
        cmd = "sed -i 's/" + self.cfg.lst_molecules[0].name + "_//g' " + self.index_all_residues  # remove the protein 
        self.cfg.tools.execute.run(cmd)
        
        f = open(self.input_file_tmp,'w')
        f.write('1 | 13\n')
        f.write("q\n")
        f.close()
        
        cmd = '{} {} -f {} -o {} -n {} < {}'.format(self.cfg.gromacs, self.cfg.graph + "make_ndx" + self.cfg.mpi,
                                              self.cfg.gro_md, self.index_all_residues, self.index_all_residues, self.input_file_tmp)
                                              
        self.cfg.tools.execute.run(cmd)
        self.cfg.tools.execute.run('rm {}'.format(self.input_file_tmp))
        
    def join_van_elect_hbonds(self, mol_target, mol_query):
        if self.cfg.p_graph_mmpbsa and self.cfg.p_graph_hbonds:
            self.cfg.tools.join_image(
                self.cfg.format_file_name_g_mmpbsa_png.format(self.cfg.prefix_results_png, mol_target.original_name, mol_query.original_name),
                self.cfg.format_file_name_hbonds_png.format(self.cfg.prefix_results_png, mol_target.original_name, mol_query.original_name),
                self.cfg.format_file_name_hbonds_mmpbsa.format(self.cfg.prefix_results_png, mol_target.original_name, mol_query.original_name))

    def generate_graph_van_elect(self, mol_target, mol_query):
        """
            Generate van elec graph, first xvg, after plot and after join it with the hbonds graph
            :param n_lig: ligand group number
            :param nom_lig: ligand name
        """
        if self.cfg.p_graph_mmpbsa:
            xvg_file = self.cfg.format_file_name_g_mmpbsa_xvg.format(self.cfg.prefix_results_xvg, mol_target.original_name, mol_query.original_name)
            png_file = self.cfg.format_file_name_g_mmpbsa_png.format(self.cfg.prefix_results_png, mol_target.original_name, mol_query.original_name)
            cmd = 'echo {0} {1} | {2} -f {3} -s {4} -mm {5} '.format(
                mol_target.group,
                mol_query.group,
                self.cfg.g_mmpbsa,
                self.cfg.xtc_md,
                self.cfg.tpr_min,
                xvg_file
            )
            self.cfg.tools.execute.run(cmd)
            cmd = '{0} {1} {2} {3} {4}'.format(
                self.cfg.python_run,
                self.cfg.graph_mmpbsa,
                xvg_file,
                "\"Mmpbsa " + mol_target.original_name + " " + mol_query.original_name + "\"",
                png_file
            )
            self.cfg.tools.execute.run(cmd)

    def generate_graph_hbonds(self, mol_target, mol_query ):
        """
           Generate hydrogen bonds graphs
        :param n_lig:   ligand group number
        :param nom_lig: ligand name
        """
        if self.cfg.p_graph_hbonds:
            xvg_file = self.cfg.format_file_name_hbonds_xvg.format(self.cfg.prefix_results_xvg, mol_target.original_name, mol_query.original_name)
            png_file = self.cfg.format_file_name_hbonds_png.format(self.cfg.prefix_results_png, mol_target.original_name, mol_query.original_name)
            cmd = 'echo {0} {1} | {2} {3} -f {4} -s {5} -num {6}'.format(
                mol_target.group,
                mol_query.group,
                self.cfg.gromacs,
                self.cfg.graph + "hbond" + self.cfg.mpi,
                self.cfg.xtc_md,
                self.cfg.tpr_min,
                xvg_file
                )
            self.cfg.tools.execute.run(cmd)
            cmd = '{0} {1} {2} {3} {4}'.format(
                self.cfg.python_run,
                self.cfg.graph_hbonds,
                xvg_file,
                "\"HBnum " + mol_target.original_name+ " " + mol_query.original_name + "\"",
                png_file
            )
            self.cfg.tools.execute.run(cmd)




    def rerun_target_query(self, mol_target, mol_query, mpd_file, out):
        """
            rerun for target_query
        """
        lst_cmd = []

        lst_cmd.append(
            '{0}\\\n\t -f {1}\\\n\t -c {2}\\\n\t -p {3}\\\n\t -o {4}\\\n\t -n {5}\\\n\t -maxwarn {6}'.format(
            self.cfg.grompp,
            mpd_file,
            self.cfg.last_equilibration,
            self.cfg.top,
            out,
            #self.cfg.index,
            self.index_all_residues,
            self.cfg.MAX_WARNINGS
            )
        )
        lst_cmd.append(
            '{0}\\\n\t  -deffnm {1}\\\n\t  {2}\\\n\t  -g {3}\\\n\t  -rerun {4} '.format(
            self.cfg.mdrun,
            out,
            self.cfg.threds,
            out+".log",
            self.cfg.xtc_md
            )
        )
        return lst_cmd

    def def_rerun_md(self):
        """
            Rerun sequentally 
        """
        mpd_file = '{}_{}.mdp'.format(os.path.splitext(self.cfg.config_file_md)[0], self.cfg.lst_molecules[0].name)
        self.cfg.tools.cp_file(self.cfg.config_file_md, mpd_file)
        energies = " ".join([ i.name for i in self.cfg.lst_molecules] )
        f_mdp = open(mpd_file, "a") 
        f_mdp.write("energygrps\t     =  "+energies +"\n")
        f_mdp.close()        
        
        lst_cmd = []
        '{0} {1}\\\n\t {2}\\\n\t {3}\\\n\t {4}\\\n\t {5}\\\n\t {6}\\\n\t {7}\\\n\t {8}'
        lst_cmd.append('{0}\\\n\t  -f {1}\\\n\t  -c {2}\\\n\t  -p {3}\\\n\t  -o {4}\\\n\t -n {5}\\\n\t  -maxwarn {6}'.format(
            self.cfg.grompp,
            mpd_file,
            self.cfg.last_equilibration,
            self.cfg.top,
            self.rerun_md+".tpr",
            self.cfg.index,
            self.cfg.MAX_WARNINGS
        ))
        lst_cmd.append('{0}\\\n\t  -deffnm {1}\\\n\t  {2}\\\n\t  -g {3}\\\n\t  -rerun {4} '.format(
            self.cfg.mdrun,
            self.rerun_md,
            self.cfg.threds,
            self.rerun_md_log,
            self.cfg.xtc_md

        ))


        return lst_cmd


