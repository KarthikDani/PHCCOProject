load('..\data\EGFR_SN.mat')
[~,~,RAWdata] = xlsread('..\data\example.xlsx'); % this is an example file

% this script needs a expression data, which is called here RAWdata in the
%format: {'O43609' 19.6484780000000  0}
% uniprot_id, expression_data and p-value

% this script also requires function affectedRxns_inhibitions, provided
% separately.

for i =1:length(EGFR_SN.genes) % remove duplicates
    try
ida=find(ismember(RAWdata(:,1),EGFR_SN.genes(i)));
if isempty(RAWdata(ida(find(cell2mat(RAWdata(ida,2)) == max(abs(unique(cell2mat(RAWdata(ida, 2))))))),2))
    RAWdata(ida(find(cell2mat(RAWdata(ida,2)) ~= -max(abs(unique(cell2mat(RAWdata(ida, 2))))))),:) = [];
else
RAWdata(ida(find(cell2mat(RAWdata(ida,2)) ~= max(abs(unique(cell2mat(RAWdata(ida, 2))))))),:) = [];
end
    end
end
 


mesen = [];
epi = [];
no = [];

% parse gene according to cut-off 2 for mesenchmal cells and -2 (or whatever is decided) for
% epithelial cells and p-value <=0.05



for i =1:length(EGFR_SN.genes)
    ida=find(ismember(RAWdata(:,1), EGFR_SN.genes(i)));
    indmesen =cell2mat(RAWdata(ida,2)) >= 2 & cell2mat(RAWdata(ida,3)) <= 0.05;
    mesen =[mesen;ida(indmesen == 1)];
    indepi =cell2mat(RAWdata(ida,2)) <= -2 & cell2mat(RAWdata(ida,3)) <= 0.05; 
    epi = [epi;ida(indepi == 1)];
    no = [no;ida(indmesen == 0);ida(indepi == 0)];
end
    
mesengenes = unique(RAWdata(mesen,1)); % genes with high expression in mesenchymal cells    
epigenes = unique(RAWdata(epi,1)); % genes with high expression in epithelial cells   
nochangegenes = unique(RAWdata(no,1)); % genes with nochange 

commonRows = ismember(nochangegenes,mesengenes,'rows');
nochangegenes(commonRows,:) = [];

commonRows = ismember(nochangegenes,epigenes,'rows');

[upregulatedE, inhibitedE]= affectedRxns_inhibition(EGFR_SN, epigenes); %reactions upregulated and inhibited in epithelial model
[upregulatedM, inhibitedM]= affectedRxns_inhibition(EGFR_SN, mesengenes); %reactions upregulated and inhibited in mesenchymal model

[minF maxF] = fastFVA(EGFR_SN);
EGFR_E = EGFR_SN;
EGFR_M = EGFR_SN;

% regulate flux bounds of both the models

EGFR_M.lb(upregulatedE) = minF(upregulatedE)/100;
EGFR_E.lb(upregulatedM) = minF(upregulatedM)/100;
EGFR_M.ub(upregulatedE) = maxF(upregulatedE)/100;
EGFR_E.ub(upregulatedM) = maxF(upregulatedM)/100;
    
   
EGFR_E.ub(inhibitedE) = EGFR_E.ub(inhibitedE)/100
EGFR_E.lb(inhibitedE) = EGFR_E.lb(inhibitedE)/100
EGFR_M.ub(inhibitedM) = EGFR_M.ub(inhibitedM)/100
EGFR_M.lb(inhibitedM) = EGFR_M.lb(inhibitedM)/100



