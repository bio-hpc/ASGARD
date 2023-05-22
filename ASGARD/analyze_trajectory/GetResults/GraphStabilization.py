#
#   Generates a table with the stabilization data
#
#
import numpy as np
import re


class GraphStabilization(object):

    def __init__(self, cfg):

        self.cfg = cfg
        if self.cfg.p_graph_stabilization:
            data = ["Potential", "Kinetic-En.", "Total-Energy", "Temperature", "Pressure", "Density", "Volume"]
            out_data = ["_potencial", "_kinetic-en", "_total-energy", "_temperature", "_presure", "_density", "_volume"]
            datos = []
            label = []

            for i in range(len(data)):
                cmd = 'echo {0} | {1} {2}{3}{4} -f {5} -s {6} -o {7}'.format(
                    data[i],
                    self.cfg.gromacs,
                    self.cfg.graph,
                    "energy",
                    self.cfg.mpi,
                    self.cfg.edr_md,
                    self.cfg.tpr_min,
                    self.cfg.prefix_results_xvg+out_data[i] + ".xvg"

                )
                self.cfg.tools.execute.run(cmd)
                x, j = self.read_file(self.cfg.prefix_results_xvg+out_data[i] + ".xvg")
                datos.append(x)
                label.append(j)
                cmd = '{} {} {} {} {}'.format(
                    self.cfg.python_run,
                    self.cfg.standar_graph_xvg,
                    self.cfg.prefix_results_xvg + out_data[i] + ".xvg",
                    out_data[i][1:],
                    self.cfg.prefix_results_png + out_data[i] + ".png",
                )
                self.cfg.tools.execute.run(cmd)
                if i == 2:
                    self.join_image(self.cfg.prefix_results_png + out_data[i - 2] + ".png", self.cfg.prefix_results_png + out_data[i - 1] + ".png",
                                   self.cfg.prefix_results_png + out_data[i] + ".png", self.cfg.out_png_stabilization_1)
                if i == 5:
                    self.join_image(self.cfg.prefix_results_png + out_data[i - 2] + ".png", self.cfg.prefix_results_png + out_data[i - 1] + ".png",
                                   self.cfg.prefix_results_png + out_data[i] + ".png", self.cfg.out_png_stabilization_2)

            f = open(self.cfg.table_stabilization, "w")
            f.write('\\begin{center} ' + "\n")
            f.write('\t\\begin{tabular}{ l r l l }' + "\n")
            f.write('\t\t\multicolumn{3}{c} {Stabilization} \\\\' + "\n")
            f.write('\t\t\hline \n')
            f.write('\t\t\t{:<15} & {:<15} &  {:<15}      \\\\\n'.format("Properties" ,"Mean","Variance"))
            f.write('\t\t\hline \n')
            print(self.cfg.prefix_results + "_table_stabilization.tex")


            for i in range(len(datos)):
                f.write('\t\t\t{:<15} & {:<15} & $\\pm$ {:<15} \\\\ \n'.format(out_data[i][1:], round(np.mean(datos[i]), 1),round(np.std(datos[i]), 1)) )
            f.write('\t\t\hline \n')
            f.write('\t\\end{tabular}' + "\n")
            f.write('\\end{center}' + "\n")
            f.close()

    def read_file(self, file):
        lst = []
        label = ""
        f = open(file)
        for i in f:
            if i.startswith("@    yaxis  label "):
                i = re.sub(' +', ' ', i).strip()  # c
                aux = i.split(" ")
                label = aux[3].strip()[2:-2]  # removes parenthesis and double quotes
                label = label.replace("^", "\\^")
            if not i.startswith("#") and not i.startswith("@"):
                i = re.sub(' +', ' ', i).strip()  # removes double spaces at the beginning and end
                aux = i.split(" ")  
                lst.append(float(aux[1]))  # Coulomb energy
        f.close()
        return lst, label

    def join_image(self, img_a, img_b, img_c, out):
        comando="convert -size 2048x1536 "+img_a+" -thumbnail 800x600 "+img_a+".png"
        self.cfg.tools.execute.run(comando)
        comando="convert -size 2048x1536 "+img_b+" -thumbnail 800x600 "+img_b+".png"
        self.cfg.tools.execute.run(comando)
        comando="convert -size 2048x1536 "+img_c+" -thumbnail 800x600 "+img_c+".png"
        self.cfg.tools.execute.run(comando)
        comando="montage  -geometry +3+2  "+img_a+".png "+img_b+".png "+img_c+".png -tile 3x1 "+out
        self.cfg.tools.execute.run(comando)
        comando="rm "+img_a+".png "+img_b+".png "+img_c+".png"
        self.cfg.tools.execute.run(comando)
