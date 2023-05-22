
profiles = {

      'DEBUGGING':      # Only debbuging part
        {
            'p_graph_hbonds': False,
            'p_graph_mmpbsa': False,
            'p_interactions_gromacs': True,
            'p_step_fluctuation': False,  # True
            'p_graph_dssp': False,
            'p_graph_sasa': False,
            'p_graph_stabilization': False,
            'p_graph_rmsd': False,
            'p_graph_rmsdf': False,
            'p_graph_distance': False,
            'p_graph_helicity': False,  
            'p_graph_gyrate': False,
            'p_table_multimolecule': False,
            'p_text_in_document': False,
            'p_document_latex': False,  # Tru
            'p_desglose': False,
            'p_sequential': False,  # false
            'p_debug': False,  # False
        },

      'TARGET_QUERY_NO_HB':
        {
            'p_graph_hbonds': False,
            'p_graph_mmpbsa': True,
            'p_interactions_gromacs': True,
            'p_step_fluctuation': True,  # True
            'p_graph_dssp': True,
            'p_graph_sasa': True,
            'p_graph_stabilization': True,
            'p_graph_rmsd': True,
            'p_graph_rmsdf': True,
            'p_graph_distance': True,
            'p_graph_helicity': False,  
            'p_graph_gyrate': True,
            'p_table_multimolecule': True,
            'p_text_in_document': True,
            'p_document_latex': True,  # Tru
            'p_desglose': True,
            'p_sequential': True,  # false
            'p_debug': False,  # False
        },

        'TARGET_QUERY':
        {
            'p_graph_hbonds': True,
            'p_graph_mmpbsa': True,
            'p_interactions_gromacs': True,
            'p_step_fluctuation': True,  # True
            'p_graph_dssp': True,
            'p_graph_sasa': True,
            'p_graph_stabilization': True,
            'p_graph_rmsd': True,
            'p_graph_rmsdf': True,
            'p_graph_distance': True,
            'p_graph_helicity': False,  
            'p_graph_gyrate': True,
            'p_table_multimolecule': True,
            'p_text_in_document': True,
            'p_document_latex': True,  # Tru
            'p_desglose': True,
            'p_sequential': True,  # false
            'p_debug': False,  # False
        },
        'TARGET':
        {
            'p_graph_hbonds': False,
            'p_graph_mmpbsa': False,
            'p_interactions_gromacs': False,
            'p_step_fluctuation': True,  # True
            'p_graph_dssp': True,
            'p_graph_sasa': True,
            'p_graph_stabilization': True,
            'p_graph_rmsd': True,
            'p_graph_rmsdf': True,
            'p_graph_distance': False,
            'p_graph_helicity': False, 
            'p_graph_gyrate': True,
            'p_table_multimolecule': False,
            'p_text_in_document': True,
            'p_document_latex': True,  # Tru
            'p_desglose': True,
            'p_sequential': True,  # false
            'p_debug': False,  # False
        },
        'TARGET_QUERIES':
        {
            'p_graph_hbonds': True,
            'p_graph_mmpbsa': True,
            'p_interactions_gromacs': True,
            'p_step_fluctuation': True,  # True
            'p_graph_dssp': True,
            'p_graph_sasa': True,
            'p_graph_stabilization': True,
            'p_graph_rmsd': True,
            'p_graph_rmsdf': True,
            'p_graph_distance': True,
            'p_graph_helicity': False,  
            'p_graph_gyrate': True,
            'p_table_multimolecule': True,
            'p_text_in_document': True,
            'p_document_latex': True,  # Tru
            'p_desglose': True,
            'p_sequential': True,  # false
            'p_debug': False,  # False
        },

        'QUERIES':
        {
            'p_graph_hbonds': True,
            'p_graph_mmpbsa': True,
            'p_interactions_gromacs': True,
            'p_step_fluctuation': False,  # True
            'p_graph_dssp': False,
            'p_graph_sasa': False,
            'p_graph_stabilization': True,
            'p_graph_rmsd': True,
            'p_graph_rmsdf': True,
            'p_graph_distance': True,
            'p_graph_helicity': False,  
            'p_graph_gyrate': False,
            'p_table_multimolecule': True,
            'p_text_in_document': True,
            'p_document_latex': True,  # Tru
            'p_desglose': True,
            'p_sequential': True,  # false
            'p_debug': False,  # False

        },
        'DNA_QUERY':
        {
            'p_graph_hbonds': True,
            'p_graph_mmpbsa': True,
            'p_interactions_gromacs': True,
            'p_step_fluctuation': True,  # True
            'p_graph_dssp': False,
            'p_graph_sasa': True,
            'p_graph_stabilization': True,
            'p_graph_rmsd': True,
            'p_graph_rmsdf': True,
            'p_graph_distance': True,
            'p_graph_helicity': False,  
            'p_graph_gyrate': False,
            'p_table_multimolecule': True,
            'p_text_in_document': True,
            'p_document_latex': True,  # Tru
            'p_desglose': True,
            'p_sequential': True,  # false
            'p_debug': False,  # False
        },
         'TEST':
        {

            'p_graph_hbonds': False,
            'p_graph_mmpbsa': False,
            'p_interactions_gromacs': True,
            'p_step_fluctuation': False,  # True
            'p_graph_dssp': False,
            'p_graph_sasa': False,
            'p_graph_stabilization': False,
            'p_graph_rmsd': False,
            'p_graph_rmsdf': False,
            'p_graph_distance': False,
            'p_graph_helicity': False, 
            'p_graph_gyrate': False,
            'p_table_multimolecule': False,
            'p_text_in_document': False,
            'p_document_latex': False,  # Tru
            'p_desglose': False,
            'p_sequential': True,  # false
            'p_debug': False,  # False
        },
        'BIPHSIC_SYSTEMS':{
            'p_graph_hbonds': True,
            'p_graph_mmpbsa': True,
            'p_interactions_gromacs': True,
            'p_step_fluctuation': False,  # True
            'p_graph_dssp': False,
            'p_graph_sasa': False,
            'p_graph_stabilization': True,
            'p_graph_rmsd': True,
            'p_graph_rmsdf': True,
            'p_graph_distance': True,
            'p_graph_helicity': False,  
            'p_graph_gyrate': False,
            'p_table_multimolecule': True,
            'p_text_in_document': True,
            'p_document_latex': True,  # Tru
            'p_desglose': True,
            'p_sequential': True,  # false
            'p_debug': False,  # False
        }

        }
class Profiles(object):
    def __init__(self, cfg):
        self.cfg=cfg

    def set_profile_cfg(self, profile_input):
        if profile_input in profiles:
            self.profile = profiles[profile_input]
            for k, v in self.profile.items():
                self.cfg.setattr(k, v)
        else:
            print ("PROFILE ERROR, the available profiles are: ")
            for k, v in profiles.items():
                print(k)
            exit()
