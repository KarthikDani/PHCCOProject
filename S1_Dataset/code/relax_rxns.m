function [relaxed_model,r,d,f,rxns_relaxed,v] = relax_rxns(model,idx_relax,min_obj,max_changes,target_flux,alpha)
% relax_rxns: Relax reactions in the model to get a feasible model (alpha = 0) 
% or bring it closer to a target flux state (0 < alpha < 1)
% 
% [relaxed_model] = relax_rxns(model) minimally relaxes bounds of reactions to make a model feasible.
%
% [relaxed_model,r] = relax_rxns(model) also returns total relaxation r.
%
% [...] = relax_rxns(model,idx_relax,min_obj,max_changes,target_flux,alpha)
% has 5 optional parameters.
%
% idx_relax: index of reactions that are to be relaxed. (default: all reactions in model)
% 
% min_obj: desired minimum value of model objective function; for eg.
% minimum biomass. (default = 1e-5), only used in models with objective
% function and with alpha > 0.
%
% max_changes: maximum number of reactions that can be relaxed. (default: length of reactions in model)
% 
% target_flux: desired flux state.
% 
% alpha: weighting parameter in the objective function (minimize (alpha)*d + (1-alpha)*r)).
% (default = 0)
%       alpha = 0: minimizes number of exchange reactions for relaxation. 
%       0 < alpha < 1: minimizes number of exchange reactions for relaxation plus distance between source and target flux state.
%       alpha = 1: minimizes distance between source and target flux state.
% 
% [relaxed_model,r,d,f,rxns_relaxed] = relax_rxns(...) gives additional
% outputs d, distance (model_flux - target_flux); f, value of objective
% function and rxns_relaxed, reactions which are relaxed in the new model
% with the new bounds.
% 
% EXAMPLE: 
% load ecoli_core_model.mat
% model.lb(find(ismember(model.rxns,'EX_glc(e)'))) = 0;
% changeCobraSolver('gurobi5');

tol = 1e-5;

if ~exist('alpha','var') || isempty(alpha)
    alpha = 0;
end

if ~exist('max_changes','var') || isempty(max_changes)
    max_changes = length(model.rxns);
end

if ~exist('min_obj','var') || isempty(min_obj)
    min_obj = tol;
end

if ~exist('idx_relax','var') || isempty(idx_relax)
    idx_relax = 1:length(model.rxns);
end

if ~exist('target_flux','var') || isempty(target_flux)
    target_flux = zeros(length(model.rxns),1);
end

% Setup optimization problem

nrxns = length(model.rxns);
lb=model.lb;
ub=model.ub;
S=model.S;

if alpha == 0 && ~isempty(find(model.c))
    sol=optimizeCbModel(model);
    assert(sol.f < tol,'Model is already feasible');
    c=zeros(nrxns,1);
    c(find(model.c))=1;
end

if max(ub)>1000
    rMax = 1000;
else
    rMax = max(ub);
end

nrelax = length(idx_relax);
idx_norelax = setdiff(1:nrxns, idx_relax);

qflag=cvx_quiet(true);
cvx_begin
    cvx_solver gurobi
    variables v(nrxns) rpos(nrelax) rneg(nrelax)
    variable rxn_switch(nrelax) binary

    if alpha == 0
        sum(rxn_switch)<= max_changes
        rpos <= rxn_switch.*rMax
        rneg <= rxn_switch.*rMax
        minimize(sum(rxn_switch))
    else

        sum(rxn_switch)<= max_changes
        rpos <= rxn_switch.*rMax
        rneg <= rxn_switch.*rMax
        
        minimize (alpha*norm(target_flux - v,2) + (1-alpha)*(sum(rxn_switch)))
    end

    subject to
    %
    S*v == 0

    lb(idx_norelax) <= v(idx_norelax) <= ub(idx_norelax)
    % Relaxing the bounds of reactions which are to be relaxed
    lb(idx_relax) - rneg <= v(idx_relax) <= ub(idx_relax) + rpos

    if alpha == 0 && ~isempty(find(model.c))
        c'*v >= min_obj;
    end
  
    rpos >= 0
    rneg >= 0

cvx_end
cvx_quiet(qflag)

if ~strcmpi(cvx_status,'Solved')
    error('Relaxtion procedure failed.')
end

rpos(abs(rpos) < tol)=0;
rneg(abs(rneg) < tol)=0;

% Relax model
relaxed_model = model;
relaxed_model.lb(idx_relax) = relaxed_model.lb(idx_relax) - rneg*1.01; % Add just a tiny bit extra
relaxed_model.ub(idx_relax) = relaxed_model.ub(idx_relax) + rpos*1.01; % to avoid infeasibilities later on

% Check if relaxed model is feasible
if alpha == 0
    solnew = optimizeCbModel(relaxed_model);
    if solnew.stat ~= 1
        fprintf('\nModel is infeasible after relaxation :-(\n')
    else
        fprintf('\nFeasible model obtained after relaxation\n')
        fprintf('Objective value (sum of relaxations) %1.4f\n', cvx_optval)
    end
end

fprintf('Number of constraints relaxed: %d (out of %d)\n', sum(rpos ~=0 | rneg ~=0), nrelax)
fprintf('\nReaction\tOriginal Model\tRelaxed Model\tv\n')
fprintf('---------------------------------------------------------\n')

rxns_relaxed = [];
for i=1:nrelax
    if rpos(i) > 0 || rneg(i) > 0
        rxns_relaxed = cat(1,rxns_relaxed,cat(2,model.rxns(idx_relax(i)),model.lb(idx_relax(i)) - rneg(i), model.ub(idx_relax(i)) + rpos(i)));
        fprintf('%s\t[%1.2f,%1.2f]\t[%1.2f,%1.2f]\t[%1.2f]\n', model.rxns{idx_relax(i)}, ...
            model.lb(idx_relax(i)), model.ub(idx_relax(i)), ...
            model.lb(idx_relax(i)) - rneg(i), model.ub(idx_relax(i)) + rpos(i), ...
            v(idx_relax(i))),
    end
end

f = cvx_optval;
d = norm(target_flux - v,2);
r = (sum(rpos) + sum(rneg));

