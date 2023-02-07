from .Profiles import profiles
import os
class Resume():
    def __init__(self, cfg):
        print("\nProfile settings:   "+cfg.profile)

        profile = profiles[cfg.profile]
        for k, v in profile.items():
            if v == 1:
                print(cfg.format_2.format(k.title(), "True"))
            else:
                print(cfg.format_2.format(k.title(), "False"))



        print ("\nResume input params:")
        print(cfg.format_2.format('step_dm', cfg.step_md+"\t"+ str(int(int(cfg.step_md)*cfg.integration_step ))+ " ps" ) )
        print(cfg.format_2.format('step_nvt', cfg.step_nvt+"\t"+ str(int(int(cfg.step_nvt)*cfg.integration_step ))+ " ps" ) )
        print(cfg.format_2.format('step_npt', cfg.step_npt+"\t"+ str(int(int(cfg.step_npt)*cfg.integration_step ))+ " ps" ) )
        print(cfg.format_2.format('step_min', cfg.step_min+"\t"+ str(int(int(cfg.step_min)*cfg.integration_step ))+ " ps" ) )
        print(cfg.format_2.format('force_field', cfg.force_field))
        print(cfg.format_2.format('solvent', cfg.solvent))
        print(cfg.format_2.format('temp', cfg.temp))
        print(cfg.format_2.format('typeGrid', cfg.type_grid))
        print(cfg.format_2.format('write_data', cfg.write_data))
        print(cfg.format_2.format('padding_grid', cfg.padding_grid))
        print(cfg.format_2.format('temp', cfg.temp))
        print(cfg.format_2.format('ph', cfg.ph))
        print(cfg.format_2.format('integration_step', cfg.integration_step))
        print(cfg.format_2.format('ensemble', cfg.ensemble))
        print(cfg.format_2.format('coulomb_type', cfg.coulomb_type))
        print(cfg.format_2.format('t_coupl', cfg.t_coupl))
        print(cfg.format_2.format('p_coupl', cfg.p_coupl))
        print(cfg.format_2.format('name_solvent', cfg.name_solvent))
        print(cfg.format_2.format('name_ions', cfg.name_ions))

        print("\nPaths:")
        print(cfg.format_2.format('tpr_min', cfg.tpr_min))
        print(cfg.format_2.format('xtc_md', cfg.xtc_md))
        ##print(cfg.format_2.format('edr_md', cfg.edr_md))
        ##print(cfg.format_2.format('trr_md', cfg.trr_md))
        print(cfg.format_2.format('md_gro', cfg.gro_md))
        print(cfg.format_2.format('md_top', cfg.top))
        print(cfg.format_2.format('prefix_molec', cfg.prefix_molec))
        print(cfg.format_2.format('grids', cfg.folder_grids))
        print(cfg.format_2.format('prefix_templates', cfg.prefix_templates))
        print(cfg.format_2.format('out_aux', cfg.out_aux))


        print(cfg.format_2.format('pdb', cfg.pdb))
        print(cfg.format_2.format('resume_file', cfg.resume_file))
        print ("\n config files")
        for i in cfg.config_files.split('\n'):
            if i != "":
                aux = os.path.splitext(i)[0]
                name = aux[aux.rfind('_')+1: ]
                if name[0].isdigit():
                    name="npt "+name
                print(cfg.format_2.format(name, i))

        print("\nMolecules/s: \t    Groups \t"+str(len(cfg.lst_molecules)))
        for mol in cfg.lst_molecules:
            print(cfg.format_2.format(mol.name, mol.group +"\t\t" + mol.original_name ))

        print ("Queue Manager \t "+cfg.template_job.queue_manager)
        print(cfg.format_2.format('run_job', cfg.template_job.run_job))
        print(cfg.format_2.format('dependency_cmd', cfg.template_job.dependency_job_cmd))
        print(cfg.format_2.format('out_job', cfg.template_job.out_job))
        print(cfg.format_2.format('err_job', cfg.template_job.err_job))

        for i in cfg.template_job.lst_header:
            print(cfg.format_2.format(i, ''))

        print ("\nOptions:")
        print(cfg.format_2.format('gromacs', cfg.gromacs))
        print(cfg.format_2.format('cmd_check', cfg.cmd_check))
        print("\n")











