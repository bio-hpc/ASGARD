import re
import numpy as np

class TableMultiMolecule:
    def __init__(self, cfg):
        self.cfg = cfg
        self.format_table = '{0:>10} & {1:>10} & {2:>10} & {3:>10} & {4:>10} & {5:>10} & {6:>10} & {7:>10} & {8:>10} & {9:>10} & {10:>10}\\\\ \n'
        if self.cfg.p_table_multimolecule:
            self.generate_table_multimol()

    def generate_table_multimol(self):
        data_print = []
        mol_target = self.cfg.lst_molecules[0]
        for i in range(1, len(self.cfg.lst_molecules)):
            mol_query = self.cfg.lst_molecules[i]
            distance_checks = []
            for cnt in range(5):
                check = "\times"
                for m_qc in mol_query.coords:
                    for m_tc in mol_target.coords:
                        aux_q = re.sub(' +', ' ', m_qc).strip().split(" ")
                        aux_t = re.sub(' +', ' ', m_tc).strip().split(" ")
                        distance = self.cfg.tools.distancia(
                            [float(aux_q[5]), float(aux_q[6]), float(aux_q[7])],
                            [float(aux_t[5]), float(aux_t[6]), float(aux_t[7])]
                        )
                        if distance < self.cfg.MAX_TABLE_MULTIMOL_DISTANCE - cnt:
                            check = "\checkmark"
                            break
                distance_checks.append(check)

            aux_tmp = []
            with open(self.cfg.f_molecule_distance_xvg.format(self.cfg.prefix_results_xvg, mol_target.original_name, mol_query.original_name)) as f:
                for line in f.read().splitlines():
                    if not line.startswith("@") and not line.startswith("#"):
                        aux_tmp.append(float(re.sub(' +', ' ', line).strip().split(" ")[1]))

            data_print.append(self.format_table.format(
                mol_target.original_name,
                mol_query.original_name,
                round(np.mean(aux_tmp), 3),
                round(np.var(aux_tmp), 3),
                round(aux_tmp[0], 3),
                round(aux_tmp[len(aux_tmp) - 1], 3),
                distance_checks[0],
                distance_checks[1],
                distance_checks[2],
                distance_checks[3],
                distance_checks[4]
            ))

        self.generate_table(data_print)

    def generate_table(self, data_print):
        with open(self.cfg.table_multimolecule, "w") as f:
            f.write('\\begin{center} ' + "\n")
            f.write('  \\begin{tabular}{ l r r r r r r r r r r }' + "\n")
            f.write('\t\multicolumn{11}{c} {Distance Protein} \\\\' + "\n")
            f.write('\t\hline \n')
            f.write("\t\t"+self.format_table.format(
                "Receptor", "Ligand", "Mean", "Variance", "Inicio", "Final", "7A", "6A", "5A", "4A", "3A"))
            f.write('\t\hline \n')
            for i in data_print:
                f.write("\t\t"+i)
            f.write('\t\hline \n')
            f.write('  \\end{tabular}' + "\n")
            f.write('\\end{center}' + "\n")
