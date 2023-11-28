#!/usr/bin/env python3
# coding: utf-8
# version: 2
#
# Generates the graphs of the dynamics, the arrays, titles, colors, etc., will be passed to it
# _______________________________________________________________________________________________________
import re
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Set the backend to 'Agg' for non-GUI environments
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.font_manager import FontProperties
import os

# Set matplotlib to use Agg backend
matplotlib.use('Agg')

FONT_SIZE = 10
EACH_XTICS = 1000  # 1 ps
PERCENT_MARGIN = 0.1  # 10%

class GenerateGraph:
    def __init__(self):
        self.dpi = 200  # dpis de las graficas

    def read_xvg(self, fichero):
        title = ""
        y_title = ""
        x_title = ""
        subtitle = ""
        x = []
        y = []
        with open(fichero) as f:
            for i in f:
                if "xaxis" in i:
                    x_title = re.sub(' +', ' ', i).strip().split("\"")[1]  # remove initial and final double spaces
                elif "yaxis" in i:
                    y_title = re.sub(' +', ' ', i).strip().split("\"")[1]  # remove initial and final double spaces
                elif " title" in i:
                    title = re.sub(' +', ' ', i).strip().split("\"")[1]  # remove initial and final double spaces
                elif "subtitle" in i:
                    subtitle = re.sub(' +', ' ', i).strip().split("\"")[1]  # remove initial double spacese and finalaux[1]
                elif not i.startswith("@") and not i.startswith("#") and i.strip() != "":
                    aux = re.sub(' +', ' ', i).strip().split(" ")  # remove initial and final double spaces
                    x.append(float(aux[0]))
                    y.append(float(aux[1]))
        return x, y, title, x_title, y_title, subtitle

    def get_fig_ax(self, xTitulo, yTitulo, titulo, ylabel2):
        fig, ax = plt.subplots()
        return fig, ax

    def get_color(self, cnt):
        clrs = ["b", "r", "g", "y", "k", "m", "c"]
        return clrs[cnt % len(clrs)] if cnt != 0 else clrs[0]

    def plot_xy_Line(self, legend, x, y, out_png, x_label, y_label, title, ylabel2):
        line_width = 1
        self.set_plt(x, y, x_label, y_label, title)
        if len(y) > 0:
            if isinstance(y[0], list):  # If it is a list type, it means it will have more data
                ax = plt.gca()
                ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=7))
                for i in range(len(y)):
                    if not isinstance(x[0], list):  # If the x's are not a list, it is assumed that they are time steps of the simulation
                        plt.plot(x, y[i], linestyle="-", linewidth=line_width, color=self.get_color(i))
                    else:
                        plt.plot(x[i], y[i], linestyle="-", linewidth=line_width, color=self.get_color(i))
            else:
                plt.plot(x, y, linestyle="-", linewidth=line_width, color=self.get_color(0))
            if legend != "":
                self.put_legend(legend,6)
            self.save_graph(out_png)

    def ticsx(self, min, max):
        diff = abs(max - min)
        m = diff * PERCENT_MARGIN
        return [min - m, max + m]

    def set_plt(self, x, y, x_label, y_label, title):
        plt.rc('font', size=10)
        # plt.axes((0.15, 0.1, 0.8, 0.8))
        plt.axes((0.15, 0.1, 0.8, 0.8))
        if isinstance(x[0], list):
            if x[0] != 0:
                min_x = min(x[0])
                max_x = max(x[0])
            else:
                min_x = 0
                max_x = max(x[0])

            if max_x < 50:
                plt.xticks(x[0], size='small', rotation=45)
        elif isinstance(x, list):
            if x[0] != 0:
                min_x = min(x)
                max_x = max(x)
            else:
                min_x = 0
                max_x = max(x)
            if max_x < 50:
                plt.xticks(x,size='small', rotation=45)
        if isinstance(y[0], list):
            min_y = min(min(val) for val in y)  # Calculate the minimum across all sublists
            max_y = max(max(val) for val in y)  # Calculate the maximum across all sublists
        else:
            min_y = min(y)
            max_y = max(y)

        plt.xlim(min_x, max_x)
        plt.ylim(self.ticsx(min_y, max_y))
        plt.xlabel(x_label)
        plt.title(title)
        plt.ylabel(y_label)

    def put_legend(self, legend, size):
        ax = plt.gca()
        if len(legend) > 0:
            for i in np.arange(len(legend)):
                ax.plot(1, 1.5, color=self.get_color(i), linewidth=2.5, linestyle="-", label=legend[i])
            ax.legend(loc='upper left', prop={'size': size})

    def graph_doble_line(self, leyenda, datosX, datosY, outPut, x_title, y_title, title, ylabel2):
        self.plot_xy_Line(leyenda, datosX, datosY, outPut, x_title, y_title, title, ylabel2)

    def line_graph(self, leyenda, datosX, datosY, outPut, xTitulo, yTitulo, titulo, ylabel2):
        self.plot_xy_Line(leyenda, datosX, datosY, outPut, xTitulo, yTitulo, titulo, ylabel2)

    def generate_histogram(self, legend, x, y, outPut, xTitulo, yTitulo, titulo):
            fig, ax = plt.subplots()
            mas_menos_graph = 0.5
            ind = np.arange(1)  # the x locations for the groups
            width = 0.1  # the width of the bars
            for i in np.arange(0, len(y)):
                ax.bar(ind + i * width, float(y[i]), width, color=self.get_color(0), edgecolor="none")
            ax.set_ylim(0, max(np.array(y, dtype=float)) + mas_menos_graph)
            ax.set_ylabel("Number of hydrogen bonds")
            ax.set_xlabel("MD Steps")
            ax.set_title(titulo)
            xmax = 1 + float(len(x)) / 10
            ax.set_xlim(-1, xmax)
            ax.axhline(y=0, xmin=-1, xmax=xmax, linewidth=0.2, color='k')

            if legend:
                self.put_legend(legend,12)

            # Set y-axis ticks to integers
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            
            self.save_graph(outPut)

    def generate_multiple_bar(self, legend, datos, resName, outPut, yTitulo, titulo):
            width = 0.2  # the width of the bars
            fig, ax = plt.subplots()
            # fig.subplots_adjust(bottom=0.35)
            
            fig.set_size_inches(13, 10)  # You can adjust these numbers as needed
            ax.set_ylabel(yTitulo,fontsize=15)
            ax.set_title(titulo,fontsize=20)
            # print("resName is",resName)
            # horizontal line at the y coordinate y 0
            ax.axhline(y=0, xmin=-1, xmax=len(datos[0]), linewidth=0.2, color='k')
            mini = []  # sirve para buscar el minimo y el maximo de los datos
            maxi = []
            tam_linea_borde_bar = [0.5]
            #
            # Cargo bars
            #
            for i in range(len(datos)):
                mini.append(min(datos[i]))
                maxi.append(max(datos[i]))
                ind = np.arange(len(datos[i]))  # list of lists
                ax.bar(ind + (width * i) + (width / 2), datos[i], width, color=self.get_color(i), edgecolor="none")
            min_y = min([min(i) for i in datos])
            max_y = max([max(i) for i in datos])
            ax.set_ylim((self.ticsx(min_y, max_y)))
            font_prop = FontProperties(size=20)
            ax.set_xticklabels(resName,fontweight='bold',fontproperties=font_prop)
            ax.set_xlabel('Residues', fontsize=14)  # Font size for axis label
            ax.tick_params(axis='both', which='minor')
            
            legend_fontsize = 10  # Choose a font size that works for your graph
            legend_columns = 1  # Adjust the number of columns as needed
            legend_properties = {'weight': 'bold', 'size': 12, 'family': 'sans-serif'}
            plt.legend(loc='upper right', prop=legend_properties, ncol=legend_columns, fancybox=True, framealpha=1, shadow=True, borderpad=1)
            plt.yticks(fontsize=12,fontweight='bold')
            self.put_legend(legend,12)
            self.putLegendBars(ax, resName, ind, width)
            self.save_graph(outPut)

    def putLegendBars(self, ax, xSubTitle, ind, width):
        ax.set_xticks(ind + width + (width / 2))
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(12)
            tick.label.set_rotation(90)

    def generate_heatmap(self, data, steps, out, title, x_label, y_label, x_ticks, yTicks, max_steps):
        plt.xlim(0, max_steps)
        plt.ylim(0, data.shape[0] - 1)
        plt.xlabel(x_label)
        plt.title(title)
        plt.ylabel(y_label)
        aspect_ratio = 0.3
        if data.shape[0] < 50:
            aspect_ratio = 5
        heatmap = plt.imshow(data, interpolation='nearest', aspect=aspect_ratio)
        plt.colorbar(heatmap, aspect=50)
        aux = max(steps) / 10
        put_step = [i * (aux) for i in range(0, 10)]
        pos_labels = range(0, max_steps, 10)
        rounded_labels = [round(val, 1) for val in put_step]
        plt.xticks(pos_labels, rounded_labels, size='small', rotation=45)
        self.save_graph(out)

    def save_graph(self, out):
        plt.savefig(out, dpi=self.dpi)

def main():
    # Example usage of GenerateGraph class
    graph = GenerateGraph()
    x, y, title, x_title, y_title, subtitle = graph.read_xvg("sample.xvg")
    legend = ["Legend1", "Legend2"]
    graph.line_graph(legend, x, y, "output.png", x_title, y_title, title, "")

if __name__ == "__main__":
    main()
