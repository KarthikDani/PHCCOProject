
load('../data/ecoli_core_model.mat')
changeCobraSolver('gurobi5');
% modelR = model;
modelR = removeRxns(model, model.rxns(strmatch('EX_', model.rxns))); % remove all the exchanges in e.coli model
modelN = modelR;
modelN.lb(findRxnIDs(modelN,modelN.rxns)) = 1; %lower bound set to 1
metList = modelN.mets;


for i=1:length(metList)
modelN = addExchangeRxn(modelN, metList(i), 0, 0);  %Exchanges were added for every reacting species  
end  

idx_relax = length(modelR.rxns):length(modelN.rxns);

[relaxed_model,r,d,f,rxns_relaxed,v] = relax_rxns(modelN,idx_relax,[], [], [], 0);
