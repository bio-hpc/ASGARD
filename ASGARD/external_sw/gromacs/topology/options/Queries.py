import glob
import os
import shutil

from options.components.tools_topology import *


class Queries():

    def __init__(self, cfg):
        self.cfg = cfg

    def execute(self):
        lst_class_queries =[]
        lst_name_queries = []
        if not os.path.isdir(self.cfg.dir_queries+"/originals"):
            os.mkdir(self.cfg.dir_queries+"/originals")


        for file_query in sorted( glob.glob(self.cfg.dir_queries + '*' + self.cfg.query_ext)):  # para la list de queries
            print  (file_query)
            name_query, ext_query = os.path.splitext(os.path.basename(file_query))
            lst_name_queries.append(name_query)
            lst_class_queries.append(self.queries_data(file_query, name_query))
            shutil.move(file_query, self.cfg.dir_queries+"/originals")

        name_out= ""
        for i in lst_name_queries:
            name_out+=i+"-"
        name_out=name_out[:-1]
        with open(os.path.join(self.cfg.dir_queries,name_out+".mol2"), 'a'):
            os.utime(os.path.join(self.cfg.dir_queries,name_out+".mol2"), None)


        prefix_out = self.cfg.dir_queries + self.cfg.out_complex + '_' +name_out

        create_target_queries_gro(lst_class_queries, prefix_out, self.cfg)

        new_topology_queries(prefix_out, self.cfg, lst_class_queries)
        create_charge_file(0, lst_class_queries, prefix_out, self.cfg)
        create_conf_file([], prefix_out, self.cfg)
        print(self.cfg.format_out_3.format('Topology generated:',
                                            str(
                                               lst_name_queries).replace('\'', '')[1:-1],
                                           bcolors.OKGREEN + 'OK' + bcolors.ENDC))


    def queries_data(self, file_query, name_query):
        conformations = check_conformations_ligand_mol2(file_query)
        if conformations == 1:
            return ( create_ligand_topology(self.cfg, file_query) )
        else:
            print(self.cfg.format_out_3.format(bcolors.FAIL + 'Topology error',
                                              os.path.basename(self.cfg.file_name_target) + "\t" + name_query,
                                              'varias conformaaciones' + bcolors.ENDC))