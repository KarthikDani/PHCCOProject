function [idxAffectedRxns, idinh]=affectedRxns_inhibition(model, geneList)
%[idxAffectedRxns, idinh]=affectedRxns_mo_inh(model, geneList)

% Returns a list of the reactions affected by the genes in geneList
%
% Based on the deleteModelGenes function.

% Find gene indices in model
[isInModel,geneInd] = ismember(geneList,model.genes);

geneInd = geneInd( find( geneInd ) );

% Find rxns associated with this gene
rxnInd = find(any(model.rxnGeneMat(:,geneInd),2));
if (~isempty(rxnInd))
    x = true(size(model.genes));
    x(geneInd) = false;
    
    constrainRxn = false(length(rxnInd),1);
    
    % Figure out if any of the reaction states is changed
    for j = 1:length(rxnInd)
        rulesplit=strread(model.rules{rxnInd(j)},'%s','delimiter','~');
        %            
        try
            if (~eval(rulesplit{1}))
                constrainRxn(j) = true;
            end
        end
    end
    idxAffectedRxns=rxnInd(constrainRxn);
end
if (~isempty(rxnInd))
    x = false(size(model.genes));
    x(geneInd) = true;
    constrainRxn1 = false(length(rxnInd),1);
    
    for j = 1:length(rxnInd)
        rulesplit=strread(model.rules{rxnInd(j)},'%s','delimiter','~');
        try
            if (eval(rulesplit{2}))
                %            
                constrainRxn1(j) = true;
            end
        end
    end
    idinh=rxnInd(constrainRxn1);

end
end