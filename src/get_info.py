# functions for MAKE DF WITH MET / RXN INFO TO ADD FBA RESULTS TO
def get_annotation_ids(x, annotation_keys):
    """"
    A function to get the annotation IDs of a reaction or metabolite.

    Parameters
    ----------
    x : cobra.core.Reaction or cobra.core.Metabolite
        A reaction or metabolite
    annotation_keys : list
        A list of annotation keys

    Returns
    -------
    annotation : list
        A list of annotation IDs

    Examples
    --------
    >>> import cobra
    >>> from src.get_info import get_annotation_ids
    >>> mod = cobra.io.read_sbml_model('models/Human1.xml')
    >>> get_annotation_ids(mod.reactions.get_by_id('EX_glc__D_e'), ['bigg.reaction', 'rhea', 'kegg.reaction', 'metanetx.reaction', 'sabiork.reaction', 'seed.reaction', 'biocyc'])

    """
    annotation = []
    for k in annotation_keys:
        if k in x.annotation.keys():
            annotation.append(x.annotation[k])
        else:
            annotation.append('')
    return annotation

def add_compartment2rxn(r):
    """
    A function to add the compartment to a reaction.
    
    Parameters
    ----------
    r : cobra.core.Reaction
        A reaction

    Returns
    -------
    rxn : str
        A reaction string with the compartment added to the metabolite IDs

    """

    met_list = [m.id for m in r.model.metabolites]
    rxn = " ".join([r.model.metabolites.get_by_id(x).name+'['+r.model.metabolites.get_by_id(x).compartment+']'\
    if x in met_list else x for x in r.reaction.split()])
    return rxn

# make RXN df 

def make_rxn_df(model):
    """
    A function to make a dataframe with the reactions in a model.
    
    Parameters
    ----------
    model : cobra.core.Model
        A model
    
    Returns
    -------
    rxns : pandas.core.frame.DataFrame
        A dataframe with the reactions in the model
        the dataframe has the following columns:
        'cell' : cell compartment
        'lb' : lower bound
        'ub' : upper bound
        'name' : reaction name
        'subsystem' : reaction subsystem
        'compartments' : compartments the reaction is in
        'reaction' : reaction string with metabolite IDs
        'reaction (IDs)' : reaction string with metabolite IDs
        'reaction (names)' : reaction string with metabolite names
        'GPR' : gene-protein-reaction rule
    """

    import pandas as pd
    # make df with rxns in model
    l = [list(r.annotation.keys()) for r in  model.reactions]
    annotation_keys = list(set([item for sublist in l for item in sublist]))
    rxns = pd.DataFrame([[get_cell_rxn(r),\
                          r.lower_bound,r.upper_bound]+get_annotation_ids(r,annotation_keys)+\
                         [r.name,r.subsystem,add_compartment2rxn(r),r.reaction,r.gpr]\
                      for r in model.reactions],index=[r.id for r in model.reactions],\
                        columns=['cell','lb','ub']+annotation_keys+\
                        ['name','subsystem','reaction','met_IDs','GPR']) 
    return rxns

# make RXN df 
def make_compact_rxn_df(model):

    """
    A function to make a dataframe with the reactions in a model.

    Parameters
    ----------
    model : cobra.core.Model
        A model
    
    Returns
    -------
    rxns : pandas.core.frame.DataFrame  
        A dataframe with the reactions in the model
        the dataframe has the following columns:
        'rxn_ID' : reaction ID
        'cell' : cell compartment
        'lb' : lower bound
        'ub' : upper bound
        'name' : reaction name
        'subsystem' : reaction subsystem
        'compartments' : compartments the reaction is in
        'reaction' : reaction string with metabolite IDs
    """


    import pandas as pd
    # make df with rxns in model
    l = [list(r.annotation.keys()) for r in  model.reactions]
    annotation_keys = list(set([item for sublist in l for item in sublist]))
    rxns = pd.DataFrame([[r.id,get_cell_rxn(r),r.lower_bound,r.upper_bound,\
                         r.name,r.subsystem,add_compartment2rxn(r)]\
                      for r in model.reactions],index=[r.id for r in model.reactions],\
                        columns=['rxn_ID','cell','lb','ub','name','subsystem','reaction']) 
    return rxns


def get_vmh_id(reaction):
    """
    
    A function to get the VMH ID of a reaction.

    Parameters
    ----------
    reaction : cobra.core.Reaction

    Returns
    -------
    vmh_id : str
        The VMH ID of the reaction


    """
    vmh_id = ''
    if 'vmhreaction' in reaction.annotation.keys(): 
        vmh_id = reaction.annotation['vmhreaction']
    return  vmh_id

def get_cell_rxn(r):

    """

    A function to get the cell compartment of a reaction.

    Parameters
    ----------
    r : cobra.core.Reaction

    Returns
    -------
    cell : str
        The cell compartment of the reaction
    """

    import numpy as np
    m1 = r.id.endswith('_RPE')
    m2 = r.id.endswith('_PR')
    m3 = r.id.endswith('_RPE_PR')
    m4 = r.id.endswith('_PR_RPE')
    m5 = r.id.endswith('_eRPE_PR')
    cell = np.select([m1, m2, m3, m4, m5], ['RPE','PR','RPE_PR','PR_RPE','RPE_PR_IF'], default=None)
    return str(cell)

def get_cell(model):
    """
    A function to get the cell compartment of a reaction.

    Parameters
    ----------
    model : cobra.core.Model

    Returns
    -------
    cell_l : list
        A list of cell compartments of the reactions in the model
    """ 

    import numpy as np
    m1 = [r.id.endswith('_RPE')  for r in  model.reactions]
    m2 = [r.id.endswith('_PR')  for r in  model.reactions]
    m3 = [r.id.endswith('_RPE_PR')  for r in  model.reactions]
    m4 = [r.id.endswith('_PR_RPE')  for r in  model.reactions]
    m5 = [r.id.endswith('_eRPE_PR')  for r in  model.reactions]
    cell_l = list(np.select([m1, m2, m3, m4, m5], ['RPE','PR','RPE_PR','PR_RPE','RPE_PR_IF'], default=None))
    return cell_l

# make RXN df 
def RPE_PR_rxn_df(model):
    """    
    A function to make a dataframe with the reactions in the RPE_PR model.

    Parameters
    ----------
    model : cobra.core.Model
        RPE_PR model
        
    Returns
    ------- 
    rxns : pandas.core.frame.DataFrame
        A dataframe with the reactions in the RPE_PR model
        the dataframe has the following columns:
        'Human1.reaction' : reaction ID
        'cell' : cell compartment
        'lb' : lower bound
        'ub' : upper bound
        'name' : reaction name
        'subsystem' : reaction subsystem
        'compartments' : compartments the reaction is in
        'reaction' : reaction string with metabolite IDs
        'reaction (IDs)' : reaction string with metabolite IDs
        'reaction (names)' : reaction string with metabolite names
        'GPR' : gene-protein-reaction rule
    
    Examples
    --------
    >>> import cobra
    >>> from src.combine_RPE_PR import RPE_PR_rxn_df
    >>> mod_RPE_PR = cobra.io.read_sbml_model('models/Human1_RPE_PR.xml')
    >>> rxns = RPE_PR_rxn_df(mod_RPE_PR)
    >>> rxns
    
    """
    import pandas as pd
    # make df with rxns in model
    l = [list(r.annotation.keys()) for r in  model.reactions]
    annotation_keys = list(set([item for sublist in l for item in sublist]))
    rxns = pd.DataFrame([[r.id,r.id.split('_')[1].replace('eRPE','RPE_PR_interface'),\
                          r.lower_bound,r.upper_bound]+get_annotation_ids(r,annotation_keys)+\
                     [r.name,r.subsystem,r.compartments,add_compartment2rxn(r),\
                      r.reaction,r.build_reaction_string(use_metabolite_names = True),r.gpr]\
                      for r in model.reactions],index=[r.id for r in model.reactions],\
                    columns=['Human1.reaction','cell','lb','ub',]+annotation_keys+\
                    ['name','subsystem','compartments','reaction','reaction (IDs)','reaction (names)','GPR']) 
    return rxns


def get_met_names(model,reaction):
    """
    get list of metabolite names from reaction id

    Parameters
    ----------
    model : cobra.Model
        model to get metabolite names from
    reaction : str
        reaction id

    Returns
    -------
    met_names : list
        list of metabolite names

    """
    met_names = []
    for m in model.reactions.get_by_id(reaction).reactants:
        met_names = met_names + [m.name]
    # if only single metabolite, return string instead of list
    if len(met_names) == 1:
        return met_names[0]
    return met_names