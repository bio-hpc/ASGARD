import shutil

from options.components.tools_topology import *
from six import unichr

class Target():
    def __init__(self, cfg):
        self.cfg = cfg

    def execute(self):

        target_charge = create_topology_target(self.cfg)
        lst_class_queries = []
        if self.cfg.target_chains == "1":  # target con un a cadena

            prefix_out = self.cfg.dir_target + self.cfg.out_complex + '_' + self.cfg.profile.lower()

            create_target_queries_gro(lst_class_queries, prefix_out, self.cfg)
            mod_target_top_queries(prefix_out, self.cfg, lst_class_queries)
            create_charge_file(target_charge, lst_class_queries, prefix_out, self.cfg)
            create_conf_file([], prefix_out, self.cfg)
            print(self.cfg.format_out_3.format('Topology generated:',
                                               os.path.basename(self.cfg.file_name_target) + "\t" ,
                                               bcolors.OKGREEN + 'OK' + bcolors.ENDC))

            os.remove(self.cfg.target_file_itp)
            os.remove(self.cfg.target_file_top)

        else:
            prefix_out = self.cfg.dir_target + self.cfg.out_complex + '_' + self.cfg.profile.lower()


            cmd = 'ls {}*porse*.itp'.format(self.cfg.file_name_target)
            itps_porse = execute_cmd(cmd)
            cmd = 'ls {}*.itp |grep -v porse'.format(self.cfg.file_name_target)
            itps = execute_cmd(cmd)
            self.cfg.itp_target_chain = []  # reiniciamos las cadenas de la prot
            cnt_letter = 97
            for i in range(len(itps)):
                letter = str(unichr(cnt_letter))
                out_itp_porse = join(self.cfg.path, self.cfg.dir_target,
                                     self.cfg.out_complex + '_' + self.cfg.profile.lower() + self.cfg.prefix_sumulation + '_porse_' + letter + '.itp')
                out_itp = join(self.cfg.path, self.cfg.dir_target,
                               self.cfg.out_complex + "_" + self.cfg.profile.lower() + self.cfg.prefix_sumulation + "_" + letter + ".itp")
                shutil.move(itps[i], out_itp)
                shutil.move(itps_porse[i],out_itp_porse)

                #   Se modifica el include porse dentro del fichero itp
                cmd = "sed -i \"s,^#include.*,#include \\\"" + out_itp_porse + "\\\",\" " + out_itp
                execute_cmd(cmd)
                self.cfg.itp_target_chain.append(out_itp)
                cnt_letter += 1


            create_target_queries_gro(lst_class_queries, prefix_out, self.cfg)
            mod_target_top_queries(prefix_out, self.cfg, lst_class_queries)
            create_conf_file([], prefix_out, self.cfg)
            create_charge_file(target_charge, lst_class_queries, prefix_out, self.cfg)
            os.remove(self.cfg.target_file_top)

            print(self.cfg.format_out_3.format('Topology generated:',
                                           os.path.basename(self.cfg.file_name_target) + "\t" +"",
                                           bcolors.OKGREEN + 'OK' + bcolors.ENDC))