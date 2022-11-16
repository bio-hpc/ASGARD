#
#	Lee un fichero y remplaza una palabra por su palabra y numero de ocurrencia
#
import glob
import os
import subprocess
import sys
 
file=sys.argv[1]
mol=""
#for i in glob.glob("*.mol2"):
cmd='cat {} |grep DB |uniq'.format(file)
molecules = subprocess.check_output(cmd, shell=True)
for i in molecules.split("\n"):
	if i.strip() != "":
		if mol == "":
			in_file=file
		else:
			in_file=mol+'_conf.tmp'
		mol = i
		cmd = 'awk \'{for(x=1;x<=NF;x++)if($x~/'+mol+'/){sub(/'+mol+'/,"'+mol+'_"++i)}}1\' '+in_file +' > ' +mol+'_conf.tmp'
		subprocess.check_output(cmd, shell=True)
		if "_conf.tmp" in in_file:
			cmd='rm '+in_file
			subprocess.check_output(cmd, shell=True)

f_name, f_ext = os.path.splitext(file)
cmd="mv {}  {}"	.format(mol+'_conf.tmp', f_name+"_conf"+f_ext)
subprocess.check_output(cmd, shell=True)

