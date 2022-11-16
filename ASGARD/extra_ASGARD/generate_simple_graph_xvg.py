#
#	Script que genera una grafica de un xvg con solo 2 columnas
#
import sys, re
import os
import matplotlib
matplotlib.use('Agg')# Force matplotlib to not use any Xwindows backend.
from pylab import *
import matplotlib.pyplot as plt

DPI = 500
LINE_WIDTH = 0.5
COLORS = ["b", "g", "r", "y", "k", "m", "c"]


def put_legend(ax, leyenda):
    if len(leyenda) > 0:
        for i in np.arange(len(leyenda)):
            ax.plot(1, 1.5, color=COLORS[i % len(COLORS)], linewidth=2.5, linestyle="-", label=leyenda[i])
        ax.legend(loc='upper left', prop={'size': 6})


def get_fig_ax(x_title, y_title, title, y_label_2):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if y_label_2 != "":
        ax2 = fig.add_subplot(111, sharex=ax, frameon=False)#//fig.add_subplot(111)
        ax2.set_ylabel(y_label_2, color="green")
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")
        ax2.tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            left='off',        # ticks along the bottom edge are off
            right='off',       # ticks along the top edge are off
            labelright='off')  # labels along the bottom edge are off

    ax.set_xlabel(x_title, color="blue")
    ax.set_ylabel(y_title, color="blue")
    plt.title(title)
    return fig, ax


def read_xvg(fichero):
    titulo = ""
    y_title = ""
    x_title = ""
    subtitle = ""
    x = []
    y = []
    f=open(fichero)
    for i in f:
        if i.find("xaxis") != -1:
            aux = re.sub(' +', ' ', i).strip().split("\"") 		#eliminamos espacios dobles inicial y final
            x_title = aux[1]
        elif i.find("yaxis") != -1:
            aux=re.sub(' +', ' ', i).strip().split("\"") 		#eliminamos espacios dobles inicial y final
            y_title = aux[1]
        elif i.find(" title") != -1:
            aux=re.sub(' +', ' ', i).strip().split("\"") 		#eliminamos espacios dobles inicial y final
            titulo = aux[1]
        elif i.find("subtitle") != -1:
            aux=re.sub(' +', ' ', i).strip().split("\"") 		#eliminamos espacios dobles inicial y final
            subtitle=aux[1]
        elif not i.startswith("@") and not i.startswith("#") and i!="":
            aux = re.sub(' +', ' ', i).strip().split(" ") 		#eliminamos espacios dobles inicial y final
            x.append(float(aux[0]))
            y.append(float(aux[1]))
    f.close()
    return x, y, titulo, x_title, y_title, subtitle


if len(sys.argv) < 3:
    print("")
    print("Parameters:")
    print("1 file xvg")
    print("2 ...")
    print("3 ...")
    print("n-2 Title")
    print("n-1 Out png")

    print ("")
    exit()


title = sys.argv[len(sys.argv)-2]
out_png = sys.argv[len(sys.argv)-1]
lst_x = []
lst_y = []
legend = []
for i in range(1, len(sys.argv)-2):
    file_xvg = sys.argv[i]
    x, y, title, x_title, y_title, _ = read_xvg(file_xvg)
    lst_x.append(x)
    lst_y.append(y)
    #legend.append(os.path.basename(os.path.splitext(file_xvg)[0]))
legend.append("Distance DNA DEPHBC")
fig, ax = get_fig_ax(x_title, y_title, title, "")
lst_aux_x=[]
lst_aux_y=[]
print len (lst_x[0])
for i in range(0,len(lst_x[0])):
    if i %1 ==0:
        lst_aux_x.append(lst_x[0][i])
        lst_aux_y.append(lst_y[0][i])

lst_x[0] = lst_aux_x
print len (lst_x[0])
lst_y[0] = lst_aux_y
if len(lst_y) > 0:
    if type(lst_y[0]) is list:				#si es un tipo lista quiere decir que tendra mas datos
        for i in range(len(lst_y)):
            if not type (lst_x[0]) is list: #si las x no son lista se supone que son laps  pasos de la simulacion
                ax.plot(lst_x[0], lst_y[i],  linestyle="-", linewidth=LINE_WIDTH, color=COLORS[i % len(COLORS)])
                print ("lst_x[0], lst_y[i]")
            else:
                ax.plot(lst_x[i], lst_y[i],  linestyle="-", linewidth=LINE_WIDTH, color=COLORS[i % len(COLORS)])
                print("lst_x[i], lst_y[i]")
    else:
        ax.plot(lst_x[0], lst_y[0],  linestyle="-", linewidth=LINE_WIDTH, color=COLORS[0])
        print("lst_x[0], lst_y[0]")
put_legend(ax, legend)
fig.savefig(out_png, dpi=DPI)
