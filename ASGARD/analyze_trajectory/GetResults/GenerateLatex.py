#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import json

f_1 = "    {}\n"
f_2 = "        {}\n"
f_3 = "            {}\n"
f_4 = "                {}\n"
f_5 = "                    {}\n"
NANOSECOND = 1000


class GenerateLatex(object):
    def __init__(self, cfg):
        self.cfg = cfg
        if self.cfg.p_document_latex:
            with open(self.cfg.json_text_tex) as f:
                self.text = json.load(f)

            self.document = open(self.cfg.document_tex, 'w')
            self.generate_head(self.document)
            self.generate_setup(self.document)
            self.generate_stability(self.document)
            self.generate_system_dinamics(self.document)
            self.generate_interactions(self.document)
            self.generate_footer(self.document)
            print(self.cfg.document_tex)

    def put_text(self, document, profile_check, text, format_x):
        if profile_check and self.cfg.p_text_in_document:
            document.write(format_x.format(r'\normalsize{'))
            document.write(format_x.format(text))
            document.write(format_x.format('}'))

    def put_figure_0_5(self, document, figure, format_x, center=None):
        document.write(format_x.format('\\begin{figure}[H]'))
        if center:
            document.write(format_x.format('\centering'))
        document.write(format_x.format('\includegraphics[width=0.5\textwidth] {' + figure + '} '))
        document.write(format_x.format(r'\end{figure}'))

    def put_figure(self, document, figure, format_x, capt=None):
        document.write(format_x.format(r'\begin{figure}[H]'))
        document.write(format_x.format('\centering'))
        document.write(format_x.format('\includegraphics[width=0.8 \\textwidth] {' + figure + '} '))
        if capt:
            document.write(format_x.format('\caption{' + capt + '}'))
        document.write(format_x.format('\\end{figure}'))

    def put_figure_2(self, document, figure_a, figure_b, capt_a, capt_b, format_x):
        document.write(format_x.format(r'\begin{figure}[H]'))
        document.write(format_x.format(r'\begin{minipage}{0.48\textwidth''}'))
        document.write(format_x.format('\centering'))
        document.write(format_x.format('\includegraphics[width=1\linewidth]{' + figure_a + '}'))
        if capt_a != "":
            document.write(format_x.format('\caption{' + capt_a + '}\label{Fig:Data1}'))
        document.write(format_x.format('\end{minipage}\hfill'))
        document.write(format_x.format(r'\begin{minipage}{0.48\textwidth}'))
        document.write(format_x.format(r'\centering'))
        document.write(format_x.format('\includegraphics[width=1\linewidth]{' + figure_b + '}'))
        if capt_b != "":
            document.write(format_x.format('\caption{' + capt_b + '}\label{Fig:Data2}'))
        document.write(format_x.format('\end{minipage}'))
        document.write(format_x.format('\end{figure}'))

    def chemical_propertie(self, document):
        if self.cfg.p_graph_stabilization:
            document.write(f_2.format('\subsection{Physical chemical properties}'))
            self.put_text(document, self.cfg.p_graph_stabilization, self.text['stabilization'], f_3)
            document.write(f_4.format('\input ' + self.cfg.table_stabilization))
            document.write(f_4.format(
                '\captionof{table}{\small{Average physical chemical properties for the simulation time, presented with standard deviation values.}}'))
            self.put_figure(document, self.cfg.out_png_stabilization_1, f_4)
            self.put_figure(document, self.cfg.out_png_stabilization_2, f_4)

    def rmsd(self, document):
        if self.cfg.p_graph_rmsd:
            document.write(f_2.format('\subsection{Root Mean Square Deviation}'))
            self.put_text(document, self.cfg.p_graph_stabilization, self.text['rmsd'], f_3)
            self.put_figure_2(document, self.cfg.f_rmsd_png.format(self.cfg.prefix_results_png),
                              self.cfg.f_rmsd_distribution_png.format(self.cfg.prefix_results_png), "", "", f_4)
            for mol in self.cfg.lst_molecules:
                img_a = self.cfg.f_molecule_rmsd_png.format(self.cfg.prefix_results_png, mol.original_name)
                img_b = self.cfg.f_molecule_rmsd_distribution_png.format(self.cfg.prefix_results_png,
                                                                          mol.original_name)
                self.put_figure_2(document, img_a, img_b, "", "", f_4)

    def rmsdf(self, document):
        if self.cfg.p_graph_rmsdf:
            document.write(f_2.format('\subsection{System flexibility}'))
            self.put_text(document, self.cfg.p_graph_stabilization, self.text['rmsdf'], f_3)
            if self.cfg.p_step_fluctuation:
                self.put_figure(document, self.cfg.out_graph_fluctuation.format(self.cfg.prefix_results_png), f_4)
            if self.cfg.p_graph_rmsdf:
                lst_figures_put = []  # create a list with the figures obtained
                for i in range(len(self.cfg.lst_molecules)):
                    if i + 1 < len(self.cfg.lst_molecules) and not i in lst_figures_put:
                        lst_figures_put.append(i)
                        lst_figures_put.append(i + 1)
                        img_a = self.cfg.f_molecule_rmsd_f_png.format(self.cfg.prefix_results_png,
                                                                     self.cfg.lst_molecules[i].original_name)
                        img_b = self.cfg.f_molecule_rmsd_f_png.format(self.cfg.prefix_results_png,
                                                                     self.cfg.lst_molecules[i + 1].original_name)
                        self.put_figure_2(document, img_a, img_b, "", "", f_4)

                    elif i not in lst_figures_put:
                        lst_figures_put.append(i)
                        img_a = self.cfg.f_molecule_rmsd_f_png.format(self.cfg.prefix_results_png,
                                                                     self.cfg.lst_molecules[i].original_name)
                        self.put_figure_0_5(document, img_a, f_4)
            if self.cfg.p_graph_gyrate:
                self.put_figure(document, self.cfg.out_png_gyrate, f_4)

    def distance(self, document):
        if self.cfg.p_graph_distance:
            document.write(f_2.format('\subsection{Distance center of mass}'))
            self.put_figure_2(document, self.cfg.f_distance_png.format(self.cfg.prefix_results_png),
                              self.cfg.f_distance_distribution_png.format(self.cfg.prefix_results_png), "", "", f_4)
            if self.cfg.p_desglose:
                mol_target = self.cfg.lst_molecules[0]
                if len(self.cfg.lst_molecules) > 2:  # if there are only 2 graphs
                    for i in range(1, len(self.cfg.lst_molecules)):
                        mol_query = self.cfg.lst_molecules[i]
                        img_a = self.cfg.f_molecule_distance_png.format(self.cfg.prefix_results_png,
                                                                        mol_target.original_name,
                                                                        mol_query.original_name)
                        img_b = self.cfg.f_molecule_distance_distribution_png.format(self.cfg.prefix_results_png,
                                                                                        mol_target.original_name,
                                                                                        mol_query.original_name)
                        self.put_figure_2(document, img_a, img_b, "", "", f_4)

    def sasa(self, document):
        if self.cfg.p_graph_sasa:
            document.write(f_2.format('\subsection{Solvent Accessible Surface}'))
            self.put_text(document, self.cfg.p_graph_stabilization, self.text['sasa'], f_3)
            self.put_figure(document, self.cfg.out_png_sasa.format(self.cfg.prefix_results_png), f_4)

    def generate_stability(self, document):
        if self.cfg.p_step_fluctuation or self.cfg.p_graph_rmsd or self.cfg.p_graph_rmsdf \
                or self.cfg.p_graph_stabilization or self.cfg.p_graph_gyrate \
                or self.cfg.p_graph_sasa:
            document.write(f_1.format('\section{System stability}'))
            self.chemical_propertie(document)
            self.rmsd(document)
            self.rmsdf(document)
            self.distance(document)
            self.sasa(document)

    def dssp(self, document):
        if self.cfg.p_graph_dssp:
            self.put_figure(document, self.cfg.out_png_ddsp.format(self.cfg.prefix_results_png), f_4,
                            "Number of amino acid residues in each secondary structure type, as defined by DSSP, along the simulation time.")
            self.put_figure(document, self.cfg.out_png_ddsp_ss, f_4,
                            'Evolution of secondary structure, calculated by DSSP, as a function of both simulation time and amino acid residue position on the polypeptidic chain')

    def generate_system_dinamics(self, document):
        document.write(f_1.format('\section{System dynamics}'))
        self.put_text(document, self.cfg.p_graph_stabilization, self.text['protein_dinamics'], f_3)
        self.dssp(document)

    def all_interactions(self, document):
        mol_target = self.cfg.lst_molecules[0]  # if there is a protein, it always is 0
        for i in range(1, len(self.cfg.lst_molecules)):
            mol_query = self.cfg.lst_molecules[i]
            document.write(f_2.format('\subsection{' + mol_target.original_name + ' ' + mol_query.original_name + '}'))
            if self.cfg.p_graph_hbonds or self.cfg.p_graph_mmpbsa:
                img_a = self.cfg.format_file_name_g_mmpbsa_png.format(self.cfg.prefix_results_png,
                                                                    mol_target.original_name,
                                                                    mol_query.original_name)
                img_b = self.cfg.format_file_name_hbonds_png.format(self.cfg.prefix_results_png,
                                                                    mol_target.original_name, mol_query.original_name)

                self.put_figure_2(document, img_a, img_b, "", "", f_4)

            if self.cfg.p_interactions_gromacs:
                results_prefix_root = '{}_{}_{}_rerun'.format(self.cfg.prefix_results, mol_target.original_name,
                                                              mol_query.original_name)
                results_prefix_png = '{}_{}_{}_rerun'.format(self.cfg.prefix_results_png, mol_target.original_name,
                                                              mol_query.original_name)
                # from graph_interactions_gromacs scripts
                n_g_global_hist_res = results_prefix_png + "_hist_global_energy_res.png"
                n_g_join_hist = results_prefix_png + "_hist_split_energy_res.png"
                n_g_global_line_res = results_prefix_png + "_line_global_energy_res.png"
                n_g_poseview = results_prefix_png + "_poseview.png"
                n_t_latex = results_prefix_root + "_table_interations.tex"

                self.put_figure_2(document, n_g_global_hist_res, n_g_join_hist, "", "", f_4)
                self.put_figure_2(document, n_g_global_line_res, n_g_poseview, "", "", f_4)
                document.write(f_4.format('\centering'))
                document.write(f_5.format('\input ' + n_t_latex))
                document.write(f_5.format(
                    '\captionof{table}{\small{Average ligand-receptor interaction energy for the simulation time, presented with standard deviation values, as well as the specific contribution of the main amino acid residues in the interaction.}}'))

    def table_multiligand(self, document):
        if self.cfg.p_table_multimolecule:
            document.write(f_4.format('\centering'))
            document.write(f_5.format('\input ' + self.cfg.table_multimolecule))

    def generate_interactions(self, document):
        if self.cfg.p_graph_hbonds or self.cfg.p_graph_mmpbsa or self.cfg.p_interactions_gromacs or self.cfg.p_table_multimolecule:
            document.write(f_1.format('\section{System interactions}'))
            self.all_interactions(document)
            self.table_multiligand(document)

    def generate_setup(self, document):

        with open(self.cfg.gro_md) as f:
            lines_gro = f.readlines()
        num_atoms = lines_gro[1].strip()
        num_sol = str(len([i for i in lines_gro if 'SOL' in i]))
        num_ions = str(len([i for i in lines_gro if self.cfg.name_ions[0] in i or self.cfg.name_ions[1] in i]))
        num_molecules = len(self.cfg.lst_molecules)
        if self.cfg.time_simulation > NANOSECOND:
            time_sim = str(self.cfg.time_simulation / NANOSECOND) + " ns"
        else:
            time_sim = str(self.cfg.time_simulation) + " ps"
        document.write(f_1.format('\section{Simulation setup}'))

        document.write(f_2.format(r'\footnotesize {'))

        # document.write(f_3.format('\\hspace{-0.5cm}'))
        document.write(f_4.format(r'\begin{tabular}{ l r | l r | l r }'))
        document.write(f_5.format('{} & {} & {} & {} & {} & {}  \\\\'.format(
            'Force field:', self.cfg.force_field,
            'Cut-off:', self.cfg.padding_grid,
            'Simulation time:', time_sim,
        )))
        document.write(f_5.format('{} & {} & {} & {} & {} & {}  \\\\'.format(
            'Number of sol:', num_sol,
            'Water model:', self.cfg.solvent,
            'Integration step:', self.cfg.integration_step,
        )))
        document.write(f_5.format('{} & {} & {} & {} & {} & {}  \\\\'.format(
            'Number of atoms:', num_atoms,
            'Temperature', self.cfg.temp,
            'Box type:', self.cfg.type_grid,
        )))
        document.write(f_5.format('{} & {} & {} & {} & {} & {}  \\\\'.format(
            'Number of ions:', num_ions,
            'Non-bonded int:', self.cfg.coulomb_type,
            'P coupling:', self.cfg.p_coupl,
        )))
        document.write(f_5.format('{} & {} & {} & {} & {} & {}  \\\\'.format(
            'Number of molecules:', num_molecules,
            'Ensemble:', self.cfg.ensemble,
            'T coupling', self.cfg.t_coupl,
        )))
        document.write(f_4.format('\\end{tabular}'))
        # document.write(f_3.format('\\hspace{-0.5cm}'))
        document.write(f_2.format('}'))

    def generate_head(self, document):
        footer = " ".join([i.original_name for i in self.cfg.lst_molecules])
        document.write('% !TeX TS-program = xelatex \n')
        document.write('\documentclass[12pt,a4paper]{article} \n')
        # document.write('\usepackage[a4paper,letterpaper,margin=1in]{geometry} \n')

        document.write(r'\usepackage{geometry} \geometry { a4paper, total = {170mm, 235mm}, left = 20mm, top = 10mm''} \n')

        document.write(r'\usepackage{fontspec''} \n')
        document.write(r'\usepackage{polyglossia''} \n')
        document.write(r'%\setdefaultlanguage{english''} \n')
        document.write(r'\usepackage{amsfonts''} % for the \checkmark command  \n')
        document.write(r'\usepackage{textpos''} \n')
        document.write(r'\usepackage{sectsty''} \n')
        document.write(r'\usepackage{multicol''} \n')
        document.write(r'\usepackage{enumitem''} \n')
        document.write(r'\usepackage{graphicx''} \n')
        document.write(r'\usepackage{color''} \n')
        document.write(r'\usepackage{titlepic''} \n')
        document.write(r'\usepackage {capt - of''} \n')
        document.write(r'\usepackage{float''}  \n')  # %for the position of the figures in the document
        document.write(r'\setlength{\columnseprule''}{1pt''} \n')
        document.write(r'\usepackage{fancyhdr''} \n')

        document.write(r'\setlength{\headheight''}{80pt}  \n')
        document.write(r'%\usepackage{lastpage''}  \n')  # %ultima pagina
        document.write(r'\pagestyle{fancy''}                    %Add footer \n''')
        document.write(r'\renewcommand{\headrulewidth''}{0pt} \\n')  # %remove header
        document.write(r'\fancyhf''{''}			\n')  # %remove header text

        document.write(r'\fancyfoot[FR,FR]{\color{blue''} \small ' + footer + '} \n')
        document.write(r'\fancyfoot[C]{\thepage''} %num pagina \n')

        document.write(r'\def\columnseprulecolor{\color{black}''} \n')
        document.write(r'\lhead{\includegraphics[width = 4cm]{' + self.cfg.logo_bio_hpc + '}} \n')
        document.write(r'\title{Simulation Report ' + footer + '} \n')
        # document.write(r'\author{bio-hpc} \n')
        document.write(r'\author{bio-hpc''} \n')
        document.write(
            r'\titlepic{ \vspace{8cm} \includegraphics[width =\textwidth]{' + self.cfg.md_diagram + '}}\n')
        document.write(r'\begin{document''} \n')
        document.write(f_1.format('\maketitle'))
        document.write(f_1.format('\\newpage'))

    def generate_footer(self, document):
        document.write('\end  {document} \n')
