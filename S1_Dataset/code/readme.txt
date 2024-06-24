Following files are present in this folder:

--------------------------------------------------

readme.txt			This file
createE&M.m			Creates EGFR_E and EGFR_M from EGFR_SN model provided in the 'data' folder for the expression data provided
affectedRxns_inhibition.m	Helper function for createE&M.m used for solving the boolean of GPRs and suggests affected reactions based on the gene 					list.	
relax_rxns.m			Optimization code for adding exchanges and minimization of distance between EGFR_E and EGFR_M
example.m			For testing optimization algorithm on e.coli model for adding exchanges

---------------------------------------------------

Instructions:
---------------------------------------------------

EGFR_E and EGFR_M for each cell lines are provided in the 'data' folder. 

createE&M.m is a script that can be used to create a new EGFR_E and EGFR_M model from EGFR_SN model provided in the 'data' folder, by constraining it with expression data. 
Example of expression data 'example.xlsx' is provided in the 'data' folder.

relax_rxns.m is the optimization algorithm used for adding exchanges and minimization of euclidean distance. 

example.m can be used for testing the optimization algorithm in e.coli model for adding exchanges provided in the 'data' folder. Exchanges present in the e.coli model can be removed to check whether the algorithm correctly adds the exchanges. In addition to adding correct exchanges in the model, this algorithm also removes the dead ends initially present in the e.coli model.

For minimzation of Euclidean distance, following command can be run if EGFR_M has to be reverted to EGFR_E.



[relaxed_model, ~, d, ~, rxns_relaxed] = relax_rxns(EGFR_M,[],[],[],VE,0.9);