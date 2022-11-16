import glob
from os.path import join
import os
from options.components.tools_topology import *
from six import unichr


class Target_queries():
    def __init__(self, cfg):
        self.cfg = cfg

    def execute(self):
        if not os.path.isdir(self.cfg.dir_queries+"/originals"):
            os.mkdir(self.cfg.dir_queries+"/originals")
        target_charge = create_topology_target(self.cfg)
        print(self.cfg.format_out_2.format('Target Charge: ', target_charge))
        print(self.cfg.format_out_2.format('Chains: ', self.cfg.target_chains))
        lst_class_queries = []
        if self.cfg.target_chains == "1":                                       # target con un a cadena
            lst_name_queries = []
            for file_query in glob.glob(self.cfg.dir_queries+'*'+self.cfg.query_ext):                   # para la list de queries
                name_query, ext_query = os.path.splitext(os.path.basename(file_query))
                lst_name_queries.append(name_query)
                lst_class_queries.append(self.queries_data(file_query, name_query))
                shutil.move(file_query, self.cfg.dir_queries + "/originals")

            name_out = ""
            for i in lst_name_queries:
                name_out += i + "-"
            name_out = name_out[:-1]
            with open(os.path.join(self.cfg.dir_queries, name_out + ".mol2"), 'a'):
                os.utime(os.path.join(self.cfg.dir_queries, name_out + ".mol2"), None)

            prefix_out = self.cfg.dir_queries + self.cfg.out_complex + '_' + name_out
            #prefix_out = self.cfg.dir_queries + self.cfg.out_complex + '_' + self.cfg.profile.lower()
            create_target_queries_gro(lst_class_queries, prefix_out, self.cfg)
            mod_target_top_queries(prefix_out, self.cfg, lst_class_queries)
            create_charge_file(target_charge, lst_class_queries, prefix_out, self.cfg)
            create_conf_file([], prefix_out, self.cfg)
            print(self.cfg.format_out_3.format('Topology generated:',
                                              os.path.basename(self.cfg.file_name_target) + "\t" + str(lst_name_queries).replace('\'','')[1:-1],
                                               bcolors.OKGREEN + 'OK' + bcolors.ENDC))
            os.remove(self.cfg.target_file_gro)
            os.remove(self.cfg.target_file_itp)
            os.remove(self.cfg.target_file_top)

        else:  # target con varias cadenas, cambia la forma de genetarar la top
            cmd = 'ls {}*porse*.itp'.format(self.cfg.file_name_target)
            itps_porse = execute_cmd(cmd)
            cmd = 'ls {}*.itp |grep -v porse'.format(self.cfg.file_name_target)
            itps = execute_cmd(cmd)


            #name_query, ext_query = os.path.splitext(os.path.basename(file_query))
            cnt_letter = 97
            self.cfg.itp_target_chain = []  # reiniciamos las cadenas de la prot
            for i in range(len(itps)):
                letter = str(unichr(cnt_letter))
                out_itp_porse = join(self.cfg.path, self.cfg.dir_queries, self.cfg.out_complex + '_' + self.cfg.profile.lower() +self.cfg.prefix_sumulation+ '_porse_' + letter + '.itp')
                out_itp = join(self.cfg.path, self.cfg.dir_queries, self.cfg.out_complex + "_" + self.cfg.profile.lower() + self.cfg.prefix_sumulation+ "_"+ letter +".itp")
                cmd = 'cp {} {}'.format(itps[i], out_itp)
                execute_cmd(cmd)
                cmd = 'cp {} {}'.format(itps_porse[i], out_itp_porse)
                execute_cmd(cmd)
                #   Se modifica el include porse dentro del fichero itp
                cmd = "sed -i \"s,^#include.*,#include \\\"" + out_itp_porse + "\\\",\" " + out_itp
                execute_cmd(cmd)
                self.cfg.itp_target_chain.append(out_itp)
                cnt_letter += 1

            lst_name_queries = []
            for file_query in glob.glob(self.cfg.dir_queries + '*' + self.cfg.query_ext):  # para la list de queries
                name_query, ext_query = os.path.splitext(os.path.basename(file_query))
                lst_name_queries.append(name_query)
                lst_class_queries.append(self.queries_data(file_query, name_query))
                shutil.move(file_query, self.cfg.dir_queries + "/originals")

            name_out = ""
            for i in lst_name_queries:
                name_out += i + "-"
            name_out = name_out[:-1]
            with open(os.path.join(self.cfg.dir_queries, name_out + ".mol2"), 'a'):
                os.utime(os.path.join(self.cfg.dir_queries, name_out + ".mol2"), None)

            prefix_out = self.cfg.dir_queries + self.cfg.out_complex + '_' + name_out
            #prefix_out = self.cfg.dir_queries + self.cfg.out_complex + '_' + self.cfg.profile.lower()
            create_target_queries_gro(lst_class_queries, prefix_out, self.cfg)
            mod_target_top_queries(prefix_out, self.cfg, lst_class_queries)
            create_conf_file([], prefix_out, self.cfg)
            create_charge_file(target_charge, lst_class_queries, prefix_out, self.cfg)
            print(self.cfg.format_out_3.format('Topology generated:',
                                               os.path.basename(self.cfg.file_name_target) + "\t" + str(
                                                   lst_name_queries).replace('\'', '')[1:-1],
                                               bcolors.OKGREEN + 'OK' + bcolors.ENDC))



            for i in range(len(itps_porse)):
                os.remove(itps[i])
                os.remove(itps_porse[i])
            os.remove(self.cfg.target_file_gro)
            os.remove(self.cfg.target_file_top)

    def queries_data(self, file_query, name_query):
        conformations = check_conformations_ligand_mol2(file_query)
        if conformations == 1:
            return ( create_ligand_topology(self.cfg, file_query) )
        else:
            print(self.cfg.format_out_3.format(bcolors.FAIL + 'Topology error',
                                              os.path.basename(self.cfg.file_name_target) + "\t" + name_query,
                                              'varias conformaaciones' + bcolors.ENDC))
            exit()


