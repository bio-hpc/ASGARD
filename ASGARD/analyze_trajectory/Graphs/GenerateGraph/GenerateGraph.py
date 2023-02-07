#!/usr/bin/env python2
# coding: utf-8
#    version: 2
#
#   Genera las graficas de las dinamicas, se le pasaran los arrys titulos colores etc
# _______________________________________________________________________________________________________
import re
import sys
import matplotlib
#from pylab import *
matplotlib.use('Agg')# Force matplotlib to not use any Xwindows backend.
import numpy as np
import matplotlib.pyplot as plt
import os
python_packages_path = os.path.join(os.path.dirname(__file__), '../../../../../external_sw/python2-packages/')
python_packages_path = os.path.realpath(python_packages_path)
os.environ['PYTHONPATH'] = python_packages_path
sys.path.append(python_packages_path)
FONT_SIZE = 10
EACH_XTICS = 1000 #1  ps
PERCENT_MARGIN = 0.1 # 10%
#try:
#
#	#import seaborn as sns
#except ImportError as e:
#	print (e)
#	pass


# Tamaño de letra
#plt.rc('font', size = 10)
# plt.xlim(-5,370) # Los valores del eje y variarán entre -5 y 370
#plt.ylim(0,1.2) # Los valores del eje y variarán entre 0 y 1.2 (limites)
# Colocamos las etiq uetas, meses, en las posiciones, dí as,
# con color azul y rotadas 45º
#plt.xticks(dias, meses,size = 'small',color = 'b',rotation = 45)
#plt.axes((0.1,0.1,0.8,0.8))

class GenerateGraph(object):
    #'/home/jdelapena/lanzador/lanzador/lanzador/externalSw/gromacs/analizarResults/Graph'
    def __init__(self):
        self.dpi = 200   # dpis de las graficas
        #self.clrs = ["b","r", "g",  "y", "k","m", "c"]

    def read_xvg(self,fichero):
        title = ""
        y_title = ""
        x_title = ""
        subtitle = ""
        x = []
        y = []
        f = open(fichero)
        for i in f:
            if i.find("xaxis")!=-1:
                x_title = re.sub(' +',' ',i).strip().split("\"")[1] 		#eliminamos espacios dobles inicial y final
            elif i.find("yaxis")!=-1:
                y_title = re.sub(' +',' ',i).strip().split("\"")[1] 		#eliminamos espacios dobles inicial y final
            elif i.find(" title")!=-1:
                title = re.sub(' +',' ',i).strip().split("\"")[1] 			#eliminamos espacios dobles inicial y final
            elif i.find("subtitle")!=-1:
                subtitle = re.sub(' +',' ',i).strip().split("\"")[1] 		#eliminamos espacios dobles inicial y finalaux[1]
            elif not i.startswith("@") and not i.startswith("#") and i != "":
                aux = re.sub(' +',' ',i).strip().split(" ") 		#eliminamos espacios dobles inicial y final
                x.append(float(aux[0]))
                y.append(float(aux[1]))
        f.close()
        return x, y, title, x_title, y_title, subtitle

    def get_fig_ax(self, xTitulo,yTitulo,titulo,ylabel2):
        #fig = plt.figure()
        fig, ax = plt.subplots()
        """
        ax = fig.add_subplot(111)
        if ylabel2!="":
            ax2 =  fig.add_subplot(111, sharex=ax, frameon=False)#//fig.add_subplot(111)
            ax2.set_ylabel(ylabel2, color="green")
            ax2.yaxis.tick_right()
            ax2.yaxis.set_label_position("right")
            ax2.tick_params(
                axis='y',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                left='off',        # ticks along the bottom edge are off
                right='off',       # ticks along the top edge are off
                labelright='off')  # labels along the bottom edge are off
        ax.set_xlabel(xTitulo, color="blue")
        ax.set_ylabel(yTitulo, color="blue")
        title(titulo)
        """
        return fig,ax

    def get_color(self, cnt):
        clrs = ["b", "r", "g", "y", "k", "m", "c"]
        return clrs[ cnt % len(clrs) ] if cnt != 0 else clrs[0]

    def plot_xy_Line(self,legend, x, y, out_png, x_label, y_label, title, ylabel2):
        #y = [ i  for i in y if len(i) >0 ]

        line_width = 1
        self.set_plt(x, y, x_label, y_label, title)
        if len(y) > 0:
            if type(y[0]) is list:				#si es un tipo lista quiere decir que tendra mas datos
                for i in range(len(y)):
                    if not type (x[0]) is list: #si las x no son lista se supone que son laps  pasos de la simulacion
                        plt.plot(x,y[i],  linestyle="-", linewidth=line_width, color=self.get_color(i) )#, label=legend[i])
                    else:
                        plt.plot(x[i], y[i],  linestyle="-", linewidth=line_width, color=self.get_color(i))
            else:
                plt.plot(x, y, linestyle="-", linewidth=line_width, color=self.get_color(0))
            if legend != "":
                self.put_legend(legend)
            self.save_graph(out_png)

    def ticsx(self, min, max):
        diff = abs(max - min)
        m = diff * PERCENT_MARGIN
        return [min - m , max + m]

    def set_plt(self, x, y, x_label, y_label, title):

        plt.rc('font', size=10)
        plt.axes((0.1, 0.1, 0.8, 0.8))
        if isinstance(x[0], list):
            if x[0] != 0:
                min_x = min (x[0])
                max_x = max (x[0])
            else:
                min_x = 0
                max_x = max(x[0])

            if max_x < 50:
                plt.xticks(x[0],  size='small',  rotation=45)
        elif isinstance(x, list):
            if x[0] != 0:
                min_x = min(x)
                max_x = max(x)

            else:
                min_x = 0
                max_x = max(x)
            if max_x < 50:
                plt.xticks(x, size='small',  rotation=45)              
        if isinstance(y[0], list):
            min_y = np.min(y)
            max_y = np.max(y)
            if isinstance(min_y, list):

                min_y = min (min_y)
                max_y = max(max_y)
        else:
            min_y = min(y)
            max_y = max(y)


        plt.xlim( min_x, max_x)
        plt.ylim(self.ticsx(min_y, max_y))
        plt.xlabel(x_label)
        plt.title(title)
        plt.ylabel(y_label)

    def put_legend(self, legend):
        """
            Leyenda
            datos. array con los datos
            colors:array colors=["g","r","b","y","k"] #los colores seran del primero al ultimo
        """
        ax =  plt.gca()
        if len(legend) > 0:
            for i in np.arange(len(legend)):
                ax.plot(1, 1.5, color=self.get_color(i), linewidth=2.5, linestyle="-", label=legend[i])
            ax.legend(loc='upper left', prop={'size': 6})


    def graph_doble_line(self, leyenda, datosX, datosY, outPut, x_title, y_title, title,ylabel2):
        self.plot_xy_Line(leyenda, datosX, datosY, outPut, x_title, y_title, title, ylabel2)

    def line_graph(self, leyenda, datosX,datosY, outPut, xTitulo,yTitulo,titulo,ylabel2):
        self.plot_xy_Line(leyenda, datosX, datosY, outPut, xTitulo, yTitulo, titulo, ylabel2)



    def generate_histogram(self, legend, x, y, outPut, xTitulo,yTitulo,titulo):
        fig, ax = plt.subplots()
        mas_menos_graph = 0.5
        ind = np.arange(1)  # the x locations for the groups
        width = 0.1       # the width of the bars
        for i in np.arange(0, len(y)):
            ax.bar(ind + i * width, float(y[i]), width, color=self.get_color(0), edgecolor="none")
        ax.set_ylim(min(np.array(y, dtype=float))- mas_menos_graph, max(np.array(y, dtype=float)) + mas_menos_graph)
        ax.set_ylabel(yTitulo)
        ax.set_title(titulo)
        xmax = 1+float(len(y))/10
        ax.set_xlim(-1, xmax)
        ax.axhline(y=0, xmin=-1, xmax=xmax, linewidth=0.2, color='k')
        ax.set_xticklabels("")

        if legend:
            self.put_legend( legend)
        self.save_graph(outPut)


    def generate_multiple_bar(self, legend, datos, resName, outPut, yTitulo,titulo):
        """
            Grafica desglosada por atomos
            Leyenda: leyenda a mostrar tienen que corresponmder con los dtos
            datos: hastable con datos a mostra
            tipoAtomo con los tipo atomos correlaccionados con datos
            colors: array colors=["g","r","b","y","k"] #los colores seran del primero al ultimo
            outPut: ruta donde se hara la grafica
        """
        width = 0.2                               # the width of the bars
        fig, ax = plt.subplots()
        ax.set_ylabel(yTitulo)
        ax.set_title(titulo)
        # linea HORIZONTAL en la coordenada  y 0
        ax.axhline(y=0, xmin=-1, xmax=len(datos[0]), linewidth=0.2, color='k')
        mini = []  # sirve para buscar el minimo y el maximo de los datos
        maxi = []
        tam_linea_borde_bar = [0.5]
        #
        # Cargo bars
        #
        for i in range(len(datos)):
            mini.append( min (datos[i] ))
            maxi.append( max (datos[i] ))
            ind = np.arange(len(datos[i]))          # es una lista de listas
            ax.bar(ind+(width*i)+(width/2), datos[i], width, color=self.get_color(i), edgecolor="none")
        ax.set_xlim(min(ind)-(len(datos)/2),max(ind)+(len(datos)/2))

        min_y = min([ min(i) for i in datos ])
        max_y = max([ max(i) for i in datos ])
        ax.set_ylim((self.ticsx(min_y, max_y)))
        #ax.set_ylim(min(mini)-0.3, max(maxi)+(5*len(datos))) #el valor 20.3 es un poco a ojo
        ax.set_xticks(ind)
        ax.set_xticklabels(resName)
        plt.legend(loc='upper left', prop={'size': 6})
        self.put_legend( legend)
        self.putLegendBars(ax, resName,ind,width)
        self.save_graph(outPut)

    def putLegendBars(self,ax, xSubTitle,ind,width):
        ax.set_xticks(ind+width+(width/2))
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(8)
            tick.label.set_rotation(90)

    def generate_heatmap(self, data, steps, out, title, x_label, y_label, x_ticks, yTicks, max_steps):
        plt.xlim(0,max_steps )
        plt.ylim(0,data.shape[0]-1 )
        plt.xlabel(x_label)
        plt.title(title)
        plt.ylabel(y_label)
        aspect_ratio = 0.3
        if data.shape[0] < 50:
            aspect_ratio = 5
        heatmap = plt.imshow(data, interpolation='nearest', aspect=aspect_ratio)
        plt.colorbar(heatmap, aspect=50)
        aux = max(steps) / 10
        put_step = [ i* (aux) for i in range(0,10) ]
        pos_labels = range(0, max_steps, 10)
        plt.xticks(pos_labels, put_step, size='small',
                   rotation=45)
        self.save_graph(out)

#    def save_graph(self, out):
#        plt.subplots_adjust(top=0.2,bottom=0.15,left=0.2)
#        plt.savefig(out, dpi=self.dpi)
    def save_graph(self, out):
        plt.savefig(out, dpi=self.dpi)
    """
    @staticmethod
    def rotate_ticks_x(ang, ax):
    
        
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(8)
            tick.label.set_rotation(ang)

    @staticmethod
    def rotate_ticks_y(ang, ax):
        for tick in ax.yaxis.get_major_ticks():
            tick.tick2On = False
            tick.label.set_fontsize(8)
            tick.label.set_rotation(ang)

    @staticmethod
    def delete_column(datos, pos):
        
      
        for i in datos:
            i.pop(pos)
    
    @staticmethod
    def put_valors(datos, ind, width, f_size):
        #
        #	Pone valores encima de las barras, no se usa
        #
        for w, _k in zip(datos, range(len(datos))):
            w = float(w)
            if w > 0:
                annotate('%.1f' % w, (width/2+(ind+_k*width), w+.1), va='bottom', ha='center', fontsize=f_size)
            else:
                annotate('%.1f' % w, (width/2+(ind+(_k*width)), w+(-.1)), va='bottom', ha='center', fontsize=f_size)
    """
