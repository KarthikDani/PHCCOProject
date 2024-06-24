
Following datafiles are supplied:

---------------------------------------------------------------------------------------------------------------------

Models				Inclusions
---------------------------------------------------------------------------------------------------------------------

1) EGFR_SN.mat			EGFR network in COBRA format with GPRs included from which EGFR_E and EGFR_M originates in each cell lines

2) D492.mat 			EGFR_E and EGFR_M models, 
				VE and VM in D492 cells

3) HMLE_slug.mat		EGFR_E_hmle_slug and EGFR_M models, 
				VE_slug and VM VM_slug in HMLE_slug cells 

4) HMLE_snail.mat		EGFR_E_hmle_snail and EGFR_M_hmle_snail models 
				VE_snail and VM_snail in HMLE_snail cells
 
5) HMLE_twist.mat		EGFR_E_hmle_twist and EGFR_M_hmle_twist models 
				VE_twist and VM_twist in HMLE_twist cells

6) MCF10A_snail.mat		EGFR_E_MCF10A_snail and EGFR_M_MCF10A_snail models 
				VE_MCF10A_snail and VM_MCF10A_snail in MCF10A_snail cells

7) MCF10A_tgfb.mat		EGFR_E_MCF10A_tgfb and EGFR_M_MCF10A_tgfb models
				VE_MCF10A_tgfb and VM_MCF10A_tgfb in MCF10A_tgfb cells

8) MCF7_mirna.mat		EGFR_E_MCF7_mirna and EGFR_M_MCF7_mirna models
				VE_mirna and VM_mirna in MCF7_mirna cells

9) MCF7_snail.mat		EGFR_E_MCF7_snail and EGFR_M_MCF7_snail models
				VE_MCF7_snail and VM_MCF7_snail in MCF7_snail cells

10) ecoli_core_model.mat	E.coli core model to check optimization algorithm

11) example.xlsx		Example expression data to be used in createE&M.m

VE and VM are the mean flux distribution of epithelial and mesenchymal cells

.mat files of each cell lines also has flux distribution at each sample points: sampleE_points and sampleM_points for epithelial and mesenchymal models respectively.
These were used for determining the mean flux distribution of each reaction (VE and VM)


------------------------------------------------------------------------------------------------------------------




