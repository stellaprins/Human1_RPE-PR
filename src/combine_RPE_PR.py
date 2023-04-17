# # make RPE_PR model using Human1 models
def make_RPE_PR_model(mod_RPE, mod_PR):
    from src.modify_model import fix_compartment_dict, add_rxns2groups, add_id_suffix
    from src.combine_RPE_PR import add_interface_RPE, add_interface_PR, fuse_RPE_PR, del_RPE_PR_dupl_rxns
    mod_PR = add_id_suffix(mod_PR, '_PR')
    mod_RPE = add_id_suffix(mod_RPE, '_RPE')
    mod_RPE = add_interface_RPE(mod_RPE, 'e_RPE')
    mod_PR = add_interface_PR(mod_PR, 'e_PR')
    mod_RPE_PR = fuse_RPE_PR(mod_RPE,mod_PR)
    mod_RPE_PR = fix_compartment_dict(mod_RPE_PR)
    mod_RPE_PR = add_rxns2groups(mod_RPE_PR)
    mod_RPE_PR = del_RPE_PR_dupl_rxns(mod_RPE_PR)
    return mod_RPE_PR

# make RPE_PR model using Recon3D models
def make_RPE_PR_model_Recon3D(mod_RPE, mod_PR):
    from src.modify_model import fix_compartment_dict, add_rxns2groups, add_id_suffix
    from src.combine_RPE_PR import add_interface_RPE, add_interface_PR, fuse_RPE_PR
    mod_PR = add_id_suffix(mod_PR, '_PR')
    mod_RPE = add_id_suffix(mod_RPE, '_RPE')
    mod_RPE = add_interface_RPE(mod_RPE, '[e]_RPE')
    mod_PR = add_interface_PR(mod_PR,'[e]_PR')
    mod_RPE_PR = fuse_RPE_PR(mod_RPE,mod_PR)
    mod_RPE_PR = fix_compartment_dict(mod_RPE_PR)
    mod_RPE_PR = add_rxns2groups(mod_RPE_PR)
    mod_RPE_PR = del_RPE_PR_dupl_rxns(mod_RPE_PR)
    return mod_RPE_PR

def add_interface_RPE(mod_RPE, e_RPE):   
    # add RPE_PR interface reactions to RPE model (rxn_RPE_PR = RPE <--> RPE_PR interface)
    # e_RPE is the string 'e_RPE' used in Human1 for the met IDs, for Recon3D it is '[e]_RPE'
    for rxns in  [r for r in mod_RPE.reactions if e_RPE in r.reaction]:     
        rxn = rxns.copy()
        rxn.id = rxn.id.replace('_RPE','_RPE_PR')
        for p in rxn.products:
            if e_RPE in p.id:
                p.id = p.id.replace(e_RPE,'e_RPE_PR')
                p.compartment =  'e_RPE_PR'
        for r in rxn.reactants:
            if e_RPE in r.id:
                r.id = r.id.replace(e_RPE,'e_RPE_PR')
                r.compartment =  'e_RPE_PR'
        mod_RPE.add_reactions([rxn])
    return mod_RPE

def add_interface_PR(mod_PR, e_PR):  
    # add RPE_PR interface reactions to PR model (rxn_PR-RPE = PR <--> RPE_PR interface)
        # e_PR is the string 'e_PR' used in Human1 for the met IDs, for Recon3D it is '[e]_PR'
    for rxn in [r for r in mod_PR.reactions if e_PR in r.reaction]:
        rxn.id = rxn.id.replace('_PR', '_PR_RPE')
        for p in rxn.products:
            if e_PR in p.id:
                p.id = p.id.replace(e_PR,'e_RPE_PR')
                p.compartment = 'e_RPE_PR'
        for r in rxn.reactants:
            if e_PR in r.id:
                r.id = r.id.replace(e_PR,'e_RPE_PR')
                r.compartment = 'e_RPE_PR'
    return mod_PR

# fuse RPE and PR model 
def fuse_RPE_PR(mod_RPE, mod_PR):  
    import pandas as pd
    # fuse models
    mod_RPE_PR = mod_RPE.copy()
    mod_PR_copy = mod_PR.copy()
    mod_RPE_PR.id = mod_RPE.id + '_' + mod_PR.id
    mod_RPE_PR.add_reactions(mod_PR_copy.reactions)
    return mod_RPE_PR

def del_RPE_PR_dupl_rxns(mod_RPE_PR): 
    import pandas as pd
    # list reaction IDs to be deleted (because their reactions are duplicated)
    df = pd.DataFrame([[r.id, r.reaction] for r in mod_RPE_PR.reactions \
                       if r.compartments == {'e_RPE_PR'}], columns=['id','reaction'])
    l  = list(df['id'][df['reaction'].duplicated()])
    # delete duplicated reactions and change reaction IDs (suffix: _eRPE_PR) 
    # of the reactions that only involve the e_RPE_PR compartment (interface)
    for r in [r for r in mod_RPE_PR.reactions if r.compartments=={'e_RPE_PR'}]: 
        if r.id in l:
            mod_RPE_PR.remove_reactions([r.id])
        elif '_RPE_PR' in r.id:
            r.id = r.id.replace('_RPE_PR','_eRPE_PR')
        elif '_PR_RPE' in r.id:
            r.id = r.id.replace('_PR_RPE','_eRPE_PR')
    return mod_RPE_PR
    
# SET EXCHANGE BOUNDS
# 'id' col with IDs and an 'upper_bound' col with bounds
def set_bounds_RPE(model,bounds):
    
    # open RPE exchange reactions
    for i in range(len(bounds)): 
        try:
            model.reactions.get_by_id(bounds['id'][i]+'_RPE').bounds = (-bounds['upper_bound'][i],1000)
        except KeyError:
            model.reactions.get_by_id(bounds['vmh_id'][i]+'_RPE').bounds = (-bounds['upper_bound'][i],1000)
    
    # close PR exchange reactions  
    rxns = [r for r in model.reactions if '_PR' in r.id if len(r.products) == 0]
    for r in rxns:
        r.bounds=(0,0)
    return model

# make RXN df 
def RPE_PR_rxn_df(model):
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