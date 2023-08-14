def add_all_EX_rxns(basis_model, model):
    """
    
    Parameters
    ----------
    basis_model : cobra.core.Model
        

    model : cobra.core.Model


    """
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

def set_exchange_bounds(model, ex_dict, RPE_PR = 'RPE'): 
    """
    set bounds for exchange reactions in model based on dict
    first look for reactions that end with _RPE
    then look for reactions that end with _PR
    then look for reactions that end with _eRPE_PR
    then look for reactions that do not have any of these endings

    Parameters
    ----------
    model : cobra.Model
        model to set bounds for
    ex_dict : dict
        dictionary with exchange reaction ids as keys and bounds as values
    RPE_PR : str, optional
        string to look for in exchange reaction ids, by default 'RPE'
        but can be set to 'PR' to open up exchange with PR

    Returns
    -------
    model : cobra.Model
        model with updated bounds

    """
    if RPE_PR == 'RPE':
        for ex in ex_dict.keys():
            if ex + '_RPE' in [r.id for r in model.reactions]:
                model.reactions.get_by_id(ex + '_RPE').bounds = ex_dict[ex] # set bounds for RPE exchange
            elif ex + '_PR' in [r.id for r in model.reactions]:
                print('no RPE exchange for ' + ex + ' in model: ' + model.id + '. PR exchange is opened instead')
                model.reactions.get_by_id(ex + '_PR').bounds = ex_dict[ex] # if no RPE exchange, set bounds for PR exchange
            elif ex + '_eRPE_PR' in [r.id for r in model.reactions]:
                print('no RPE exchange for ' + ex + ' in model: ' + model.id + '. PR exchange is opened instead')
                model.reactions.get_by_id(ex + '_eRPE_PR').bounds = ex_dict[ex] # if no RPE or PR exchange, set bounds for eRPE_PR exchange
            elif ex in [r.id for r in model.reactions]:
                print('no PR or RPE exchange reaction for ' + ex + ' in model: ' + model.id + '. Generic exchange is opened instead.')
                model.reactions.get_by_id(ex).bounds = ex_dict[ex] # if no RPE, PR, or eRPE_PR exchange, set bounds for exchange (generic)
            else:  
                print('no exchange reaction for ' + ex + ' in model: ' + model.id) # if no exchange reaction, print message

    elif RPE_PR == 'PR':
        for ex in ex_dict.keys():
            if ex + '_PR' in [r.id for r in model.reactions]:
                model.reactions.get_by_id(ex + '_PR').bounds = ex_dict[ex] # set bounds for RPE exchange
            elif ex + '_eRPE_PR' in [r.id for r in model.reactions]:
                model.reactions.get_by_id(ex + '_eRPE_PR').bounds = ex_dict[ex] # if no RPE or PR exchange, set bounds for eRPE_PR exchange
            elif ex + '_RPE' in [r.id for r in model.reactions]:
                print('no PR exchange reaction for ' + ex + ' in model: ' + model.id + '. RPE exchange is opened instead.') # if no PR exchange reaction, print message
                model.reactions.get_by_id(ex + '_PR').bounds = ex_dict[ex] # if no RPE exchange, set bounds for PR exchange
            elif ex in [r.id for r in model.reactions]:
                print('no PR or RPE exchange reaction for ' + ex + ' in model: ' + model.id + '. Generic exchange is opened instead.')
                model.reactions.get_by_id(ex).bounds = ex_dict[ex] # if no RPE, PR, or eRPE_PR exchange, set bounds for exchange (generic)
            else:  
                print('no exchange reaction for ' + ex + ' in model: ' + model.id) # if no exchange reaction, print message
    else:
        print('third input (RPE_PR) must be either "RPE" or "PR"')
    return model 
