import glob
import shutil
from os.path import join
import os
from options.components.tools_topology import *
from six import unichr


class Target_query():
    def __init__(self, cfg):
        self.cfg = cfg

    def execute(self):
        target_charge = create_topology_target(self.cfg)
        print(self.cfg.format_out_2.format('Target Charge: ', target_charge))
        print(self.cfg.format_out_2.format('Chains: ', self.cfg.target_chains))
        if self.cfg.target_chains == "1":                                       # target con un a cadena
            for file_query in glob.glob(self.cfg.dir_queries+'*'+self.cfg.query_ext) :                   # para la list de queries
                name_query, ext_query = os.path.splitext(os.path.basename(file_query))
                self.create_topol_target_query(file_query,name_query,target_charge)


                #    self.cfg.num_queries = create_complex_topology(self.cfg, file_query )
            os.remove(self.cfg.target_file_gro)
            os.remove(self.cfg.target_file_itp)
            # os.remove(target_file_top)

        else:  # target con varias cadenas, cambia la forma de genetarar la top
            cmd = 'ls {}*porse*.itp'.format(self.cfg.file_name_target)
            itps_porse = execute_cmd(cmd)
            cmd = 'ls {}*.itp |grep -v porse'.format(self.cfg.file_name_target)
            itps = execute_cmd(cmd)

            for file_query in glob.glob(self.cfg.dir_queries + '*' + self.cfg.query_ext):  # para la list de queries
                name_query, ext_query = os.path.splitext(os.path.basename(file_query))
                cnt_letter = 97
                self.cfg.itp_target_chain = []  # reiniciamos las cadenas de la prot
                for i in range(len(itps)):
                    letter = str(unichr(cnt_letter))
                    out_itp_porse = join(self.cfg.path, self.cfg.dir_queries, self.cfg.out_complex + '_' + name_query +self.cfg.prefix_sumulation+ '_porse_' + letter + '.itp')
                    out_itp = join(self.cfg.path, self.cfg.dir_queries, self.cfg.out_complex + "_" + name_query + self.cfg.prefix_sumulation+ "_"+ letter +".itp")

                    shutil.copyfile(itps[i], out_itp)
                    shutil.copyfile(itps_porse[i], out_itp_porse)
                    #   Se modifica el include porse dentro del fichero itp

                    cmd = "sed -i \"s,^#include.*,#include \\\"" + out_itp_porse + "\\\",\" " + out_itp
                    execute_cmd(cmd)
                    self.cfg.itp_target_chain.append(out_itp)
                    cnt_letter += 1
                self.create_topol_target_query(file_query, name_query, target_charge)
                #num_queries = create_complex_topology( file_query, num_queries, self.cfg.file_name_target)

            for i in range(len(itps_porse)):
                os.remove(itps[i])
                os.remove(itps_porse[i])
            os.remove(self.cfg.target_file_gro)
            os.remove(self.cfg.target_file_top)

    def create_topol_target_query(self, file_query, name_query,target_charge):
        conformations = check_conformations_ligand_mol2(file_query)
        if conformations == 1:
            prefix_out = self.cfg.dir_queries + self.cfg.out_complex + '_' + name_query
            lst_class_queries = []
            lst_class_queries.append(create_ligand_topology(self.cfg, file_query))
            create_target_queries_gro(lst_class_queries, prefix_out, self.cfg)
            mod_target_top_queries(prefix_out, self.cfg, lst_class_queries)

            create_charge_file(target_charge, lst_class_queries, prefix_out, self.cfg)
            create_conf_file( [],prefix_out, self.cfg)
            print(self.cfg.format_out_3.format('Topology generated:',
                                               os.path.basename(self.cfg.file_name_target) + "\t" + name_query,
                                               bcolors.OKGREEN + 'OK' + bcolors.ENDC))
        else:
            print(self.cfg.ormat_out_3.format(bcolors.FAIL + 'Topology error',
                                              os.path.basename(self.cfg.file_name_target) + "\t" + name_query,
                                              'varias conformaaciones' + bcolors.ENDC))


