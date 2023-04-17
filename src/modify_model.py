def add_all_EX_rxns(basis_model, model):
    EX_rxns = [r for r in basis_model.reactions if len(r.products)==0]
    EX_rxns_model = [r.id for r in model.reactions if len(r.products)==0]
    missing_rxns = [r.copy() for r in EX_rxns if r.id not in EX_rxns_model]
    model.add_reactions(missing_rxns)
    return model
    
def remove_compartment(model, compartment):
    # remove compartment from model
    rxns_n = [r for r in model.reactions if compartment in r.compartments]
    mets_n = [m for m in model.metabolites if compartment in m.compartment]
    model.remove_reactions(rxns_n)
    model.remove_metabolites(mets_n)
    return model

def fix_compartment_dict(model):
    # add values to dictionary keys (necessary for valid sbml model)
    new_compartment_dict = model.compartments.copy()
    for k in list(new_compartment_dict.keys()):
        new_compartment_dict[k] = k 
    model.compartments = new_compartment_dict
    return model

def add_rxns2groups(model):
    # add all reactions to their subsystem group
    g_dict = {g.name:g.id for g in model.groups}
    for r in model.reactions:
        model.groups.get_by_any(g_dict[r.subsystem])[0].add_members([r])
    return model

def add_id_suffix(model, suffix):
    # to rxns and mets
    for r in model.reactions:
        r.id = r.id + suffix
    for m in model.metabolites:
        m.id = m.id + suffix
        m.compartment = m.compartment + suffix
    return model

def close_PR_EX(model):
    # close PR exchange reactions  
    rxns = [r for r in model.reactions if '_PR' in r.id if len(r.products) == 0]
    for r in rxns:
        r.bounds=(0,0)
    return model

def close_EX(model):
    # close all exchange reactions  
    rxns = [r for r in model.reactions if len(r.products) == 0]
    for r in rxns:
        r.bounds=(0,0)
    return model

def open_RPE_EX_ub(model):
    # open RPE exchange reactions upper boundaries (allowing things to move out of the system) 
    rxns = [r for r in model.reactions if '_RPE' in r.id if len(r.products) == 0]
    for r in rxns:
        r.bounds=(0,1000)
    return model

def open_PR_EX(model,lb,ub):
    # close PR exchange reactions  
    rxns = [r for r in model.reactions if '_PR' in r.id if len(r.products) == 0]
    for r in rxns:
        r.bounds=(lb,ub)
    return model

def create_permutation_dicts(rxn_id_bounds_dict):
    import itertools
    keys, values = zip(*rxn_id_bounds_dict.items())
    permutations_dicts = [dict(zip(keys, v)) for v in itertools.product(*values)]
    return permutations_dicts

def change_bounds(model,rxn_bounds_dict):  
    for d in dictionaries:
        for k in d:
            model.reactions.get_by_id(k).bounds = d[k]
    return model

