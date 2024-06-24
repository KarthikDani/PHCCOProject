import cobra
from cobra import Model

def affected_rxns_inhibition(model, gene_list):
    """
    Returns a list of the reactions affected by the genes in gene_list.

    Parameters:
    - model: The COBRApy model.
    - gene_list: List of gene IDs.

    Returns:
    - idx_affected_rxns: Indices of affected reactions.
    - id_inh: Indices of inhibited reactions.
    """
    # Find gene indices in the model
    gene_indices = [i for i, gene in enumerate(model.genes) if gene.id in gene_list]
    
    # Find reactions associated with these genes
    rxn_indices = [i for i, rxn in enumerate(model.reactions) if any(gene_indices)]
    
    # Initialize lists to store affected and inhibited reactions
    idx_affected_rxns = []
    id_inh = []

    for rxn_idx in rxn_indices:
        rxn = model.reactions[rxn_idx]
        gene_in_rxn = [gene.id in gene_list for gene in rxn.genes]
        
        # Evaluate the reaction's gene-reaction rule
        if rxn.gene_reaction_rule:
            affected = eval_gene_reaction_rule(rxn.gene_reaction_rule, gene_list)
            if affected:
                idx_affected_rxns.append(rxn_idx)
            if not affected:
                id_inh.append(rxn_idx)

    return idx_affected_rxns, id_inh

def eval_gene_reaction_rule(rule, gene_list):
    """
    Evaluates a gene reaction rule against the provided gene list.
    
    Parameters:
    - rule: The gene reaction rule.
    - gene_list: List of gene IDs.
    
    Returns:
    - Boolean: True if the rule evaluates to true, False otherwise.
    """
    # Replace gene IDs in rule with their presence in the gene_list
    for gene in gene_list:
        rule = rule.replace(gene, 'True')
    rule = rule.replace(' and ', ' & ').replace(' or ', ' | ').replace(' not ', ' ~ ')
    
    try:
        return eval(rule)
    except Exception as e:
        print(f"Error evaluating rule {rule}: {e}")
        return False

# Example usage
model = cobra.io.read_sbml_model('path_to_model.xml')
gene_list = ['gene1', 'gene2', 'gene3']
idx_affected_rxns, id_inh = affected_rxns_inhibition(model, gene_list)

print("Affected reactions:", idx_affected_rxns)
print("Inhibited reactions:", id_inh)
