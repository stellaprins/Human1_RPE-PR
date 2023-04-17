# functions for MAKE DF WITH MET / RXN INFO TO ADD FBA RESULTS TO
def get_annotation_ids(x, annotation_keys):
    annotation = []
    for k in annotation_keys:
        if k in x.annotation.keys():
            annotation.append(x.annotation[k])
        else:
            annotation.append('')
    return annotation

def add_compartment2rxn(r):
    met_list = [m.id for m in r.model.metabolites]
    rxn = " ".join([r.model.metabolites.get_by_id(x).name+'['+r.model.metabolites.get_by_id(x).compartment+']'\
    if x in met_list else x for x in r.reaction.split()])
    return rxn

# make RXN df 
def make_rxn_df(model):
    import pandas as pd
    # make df with rxns in model
    l = [list(r.annotation.keys()) for r in  model.reactions]
    annotation_keys = list(set([item for sublist in l for item in sublist]))
    rxns = pd.DataFrame([[r.id,get_cell_rxn(r),\
                          r.lower_bound,r.upper_bound]+get_annotation_ids(r,annotation_keys)+\
                         [r.name,r.subsystem,add_compartment2rxn(r),r.reaction,r.gpr]\
                      for r in model.reactions],index=[r.id for r in model.reactions],\
                        columns=['rxn_ID','cell','lb','ub']+annotation_keys+\
                        ['name','subsystem','reaction','met_IDs','GPR']) 
    return rxns


def get_vmh_id(reaction):
    vmh_id = ''
    if 'vmhreaction' in reaction.annotation.keys(): 
        vmh_id = reaction.annotation['vmhreaction']
    return  vmh_id

def get_cell_rxn(r):
    import numpy as np
    m1 = r.id.endswith('_RPE')
    m2 = r.id.endswith('_PR')
    m3 = r.id.endswith('_RPE_PR')
    m4 = r.id.endswith('_PR_RPE')
    m5 = r.id.endswith('_eRPE_PR')
    cell = np.select([m1, m2, m3, m4, m5], ['RPE','PR','RPE_PR','PR_RPE','RPE_PR_IF'], default=None)
    return str(cell)

def get_cell(model):
    import numpy as np
    m1 = [r.id.endswith('_RPE')  for r in  model.reactions]
    m2 = [r.id.endswith('_PR')  for r in  model.reactions]
    m3 = [r.id.endswith('_RPE_PR')  for r in  model.reactions]
    m4 = [r.id.endswith('_PR_RPE')  for r in  model.reactions]
    m5 = [r.id.endswith('_eRPE_PR')  for r in  model.reactions]
    cell_l = list(np.select([m1, m2, m3, m4, m5], ['RPE','PR','RPE_PR','PR_RPE','RPE_PR_IF'], default=None))
    return cell_l

    
# df= pd.DataFrame([[get_vmh_id(r),r.name,r.subsystem,r.reaction,r.build_reaction_string(use_metabolite_names = True)] \
#              for r in mod_RPE_PR.reactions],\
#                 index = [r.id  for r in mod_RPE_PR.reactions],
#             columns=['vmh id', 'name', 'subsystem', 'reaction (id)', 'reaction (name)'])

## match VMH_IDs from bounds file to Human1 IDs to create 'R3D301_EX_rxns_RPE_PR.xlsx'

## import excel file with info on bounds
#bounds = pd.read_excel(folder / 'rxn_bounds' / 'files PL' / 'R3D301_EX_rxns.xlsx')

## vmh_id lookup table + look up RPE rxn
#df_vmh_id = pd.DataFrame([[get_vmh_id(r),r.id,r.id+'_RPE'] for r in mod.reactions if 'EX_' in get_vmh_id(r)],columns=['vmh_id','human1_id','id'])
#merged_df = df_vmh_id.merge(bounds, how = 'right', on = ['vmh_id'])
#merged_df.to_clipboard(excel=True, sep=None)
#merged_df