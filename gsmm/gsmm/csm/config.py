# Define model names and their paths
# model_paths = {
#     "ModelE": "/Users/karthik/Desktop/PHCCO IISc Internship/Models/Epithelial_csm.xml",
#     "ModelM": "/Users/karthik/Desktop/PHCCO IISc Internship/Models/Mesenchymal_csm.xml",
#     "ModelMF": "/Users/karthik/Desktop/PHCCO IISc Internship/Models/mesenchymal_fasting_integrated_csm.xml",
# }

# base_model_path = "base_model.xml"  # Path to save the filtered SBML model, i,e base model to carry out reconstruction

metabolites_of_interest = [
    "ala__L_c", "arg__L_c", "asn__L_c", "asp__L_c", "cys__L_c", "gln__L_c",
    "glu__L_c", "gly_c", "his__L_c", "ile__L_c", "leu__L_c", "lys__L_c",
    "met__L_c", "phe__L_c", "pro__L_c", "ser__L_c", "thr__L_c", "trp__L_c",
    "tyr__L_c", "val__L_c", "damp_c", "dcmp_c", "dgmp_c", "dtmp_c",
    "cmp_c", "gmp_c", "ump_c", "amp_c", "glygn2_c", "sphmyln_hs_c",
    "chsterol_c", "xolest_hs_c", "mag__hs_c", "dag_hs_c", "pail_hs_c",
    "pe_hs_c", "ps_hs_c", "pchol_hs_c", "lpchol_hs_c", "clpn_hs_c",
    "pa_hs_c", "hdcea_c", "hdca_c", "ocdcea_c", "ocdca_c", "ptrc_c",
    "spmd_c", "sprm_c", "gthrd_c", "nad_c", "nadp_c", "Q10_c",
    "paps_c", "thbpt_c", "crn_c", "atp_c", "adp_c", "pi_c",
    "h2o_c", "h_c",
]

gene_id_column = "Gene_ID"  # Column name in expression data containing gene IDs

flux_filepath = 'flux_data.pkl'  # Filepath for storing flux data

sink_flux_filepath = 'sink_flux_data.pkl'  # Filepath for storing sink flux data

output_dir__reaction_deletion_results = 'reaction_deletion_results'  # Directory for storing reaction deletion results

output_dir = '.'  # Output directory for various analysis outputs