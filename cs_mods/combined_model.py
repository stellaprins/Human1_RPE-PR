
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