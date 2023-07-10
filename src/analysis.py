# perform FVA and FBA on a model setting RPE_ATP demand to fixed values

def fba_analysis(model, boundary_dicts, objective):
    
    # inputs: model, list of boundary dicts {rxnID:(lb,ub)}, objective function rxnID (string)
    import pandas as pd
    from datetime import datetime
    from cobra.flux_analysis import flux_variability_analysis, pfba
    from src.analysis import FVA_FBA_analysis
    from src.get_info import make_rxn_df, make_compact_rxn_df
    from datetime import datetime

    #  create empty dicts
    bounds = dict() # changed model bounds
    ovs = dict() # objective values
    f = dict() # fba fluxes
    pf = dict() # pfba fluxes
    f_min = dict() # fva min
    f_max = dict() # fva max
    uptake = dict() # uptake fluxes
    secretion = dict() # secretion fluxes
        
    #   set counter, i, to 1
    i = 1
        
    with model.copy() as m:
        # set objective
        m.objective = objective
        # run analysis for every dict in list
        for d in boundary_dicts:
            for k in d:   # set bounds for all keys (rxnIDs in dict)
                m.reactions.get_by_id(k).bounds = d[k] # set bounds defined in dict

            # run analysis
            fba = m.optimize() # fba
            p_fba = pfba(m) # pfba
            # write analysis info / results into dicts
            
            # model bounds
            bounds[i] = pd.DataFrame(d, index = ['lb', 'ub']).T    
            # objective values
            ovs[i] = pd.DataFrame([objective, fba.objective_value], columns = [m.objective.direction], index = ['ID', 'value']).T   
            # uptake / secretion
            uptake_summary = m.summary().uptake_flux
            secretion_summary = m.summary().secretion_flux
            uptake[i] = pd.DataFrame(uptake_summary['flux'])
            secretion[i] = pd.DataFrame(secretion_summary['flux'])
            # all fluxes
            f[i] = fba.to_frame()['fluxes']
            pf[i] = p_fba.to_frame()['fluxes']
            # update counter
            i=i+1
    
    # prepare dfs for excel sheet
    model_info = pd.DataFrame([m.id, m.name,m.compartments,m.annotation],\
             index=['id','name','compartments','annotation'],columns = ['model'])
    rxn_df = make_rxn_df(m)   
    compact_rxn_df = make_compact_rxn_df(m)   
    met_df = pd.DataFrame([[mi.name,mi.compartment,mi.formula,mi.charge,[r.id for r in list(mi.reactions)]] for mi in m.metabolites],\
             index=[m.id for m in m.metabolites],columns=['name','compartment','formula','charge','reactions'])
    bounds_df = pd.concat(bounds)
    ovs_df = pd.concat(ovs)
    
    # uptake / secretion dfs
    uptake_mets = pd.DataFrame([[met,m.metabolites.get_by_id(met).name] for met in uptake_summary['metabolite']],\
             index=uptake_summary.index,columns=['met_id','met_name'])
    secretion_mets = pd.DataFrame([[met,m.metabolites.get_by_id(met).name] for met in secretion_summary['metabolite']],\
             index=secretion_summary['metabolite'].index,columns=['met_id','met_name'])
    uptake_df = pd.merge(uptake_mets, pd.concat(uptake,axis=1), left_index=True, right_index=True)
    secretion_df = pd.merge(secretion_mets, pd.concat(secretion,axis=1), left_index=True, right_index=True)   
    # sort uptake / secretion dfs
    uptake_df = uptake_df.sort_values(by=[c for c in uptake_df.columns if 'flux' in c],ignore_index=True,ascending=False)
    secretion_df = secretion_df.sort_values(by=[c for c in secretion_df.columns if 'flux' in c],ignore_index=True,ascending=True)
    
    # fluxes df
    fluxes_df = pd.concat([pd.DataFrame(pf), pd.DataFrame(f)],\
          keys=["parsimonious flux","flux"],axis=1)
    fluxes_df = fluxes_df =pd.merge(compact_rxn_df, fluxes_df,left_index=True, right_index=True)
    # sort table on absolute flux size 
    fluxes_df= fluxes_df.reindex(fluxes_df[[c for c in fluxes_df.columns if 'flux' in c]].abs().sort_values(by=[c for c in fluxes_df.columns if 'flux' in c],ascending=False).index)
    # select internal fluces only
    fluxes_df[fluxes_df.index.isin([r.id for r in m.reactions if len(r.products)>0])]
    
    # date stamp
    datestr = datetime.strftime(datetime.now(), '%H%M_%d-%m-%Y')   
    
    # write excel file
    with pd.ExcelWriter('results_' + datestr + '.xlsx') as writer:  
        model_info.to_excel(writer, sheet_name = 'model')
        rxn_df.to_excel(writer, sheet_name = 'reactions')
        met_df.to_excel(writer, sheet_name = 'metabolites')
        bounds_df.to_excel(writer, sheet_name = 'altered_bounds')   
        ovs_df.to_excel(writer, sheet_name = 'objective_values')  
        uptake_df.to_excel(writer, sheet_name = 'uptake')    
        secretion_df.to_excel(writer, sheet_name = 'secretion')    
        fluxes_df.to_excel(writer, sheet_name = 'fluxes')   
    return bounds, ovs, f, pf, f_min, f_max, uptake, secretion



def FVA_FBA_analysis(model,DM_atp_c__PR_rxn_id):
    from cobra.flux_analysis import flux_variability_analysis
    from src.get_info import make_rxn_df
    import pandas as pd
        
    # make df with rxns in model
    #rxns = make_rxn_df(model)
    
    df_FVA_full = pd.DataFrame()
    df_FBA_fluxes = pd.DataFrame()
    df_FBA_costs = pd.DataFrame()
    df_pFBA_fluxes = pd.DataFrame()
    df_pFBA_costs = pd.DataFrame()
    
    RPE_ATP = [0, 40, 80]

    for x in RPE_ATP:
        model.reactions.get_by_id(DM_atp_c__PR_rxn_id).bounds=(x,x) # set bounds
        try:
            fba_results = model.optimize() # traditional FBA
            fba = fba_results.to_frame() # df FBA results
            fba.columns = fba.columns + '_' + str(x) # add RPE_ATP to FBA column names
            df_FBA_fluxes = pd.concat([df_FBA_fluxes, fba.iloc[:,0]], axis=1) # fill df flux results
            df_FBA_costs = pd.concat([df_FBA_costs, fba.iloc[:,1]], axis=1) # fill df reduced costs
        except:
            pass
        
        #try:
            #fva_full = flux_variability_analysis(model, loopless=False) # FVA 
            #fva_full.columns = fva_full.columns + '_' + str(x) # add RPE_ATP to FVA column names
            #df_FVA_full = pd.concat([df_FVA_full, fva_full], axis=1) # fill df FVA results
        #except:
           # pass
        df_l = [df_FVA_full, df_FBA_fluxes]
        # df_l = [df_FVA_full, df_FBA_fluxes, df_FBA_costs]
        #df_results = rxns
        df_results = pd.DataFrame()
        df_results = pd.concat(df_l, axis=1) 
    # fix column order
    min_max_flux = [i for i in list(df_results.columns) if 'min' in i] + \
                [i for i in list(df_results.columns) if 'max' in i] + \
                [i for i in list(df_results.columns) if 'flux' in i] 

    df_results = df_results.reindex([i for i in list(df_results.columns) if i not in min_max_flux] + min_max_flux, axis = 1)
    
    return df_results

# mod_RPE_PR_results = FVA_FBA_analysis(mod_RPE_PR,'MAR03964_RPE')
#mod_Human1_Human1_results = FVA_FBA_analysis(mod_Human1_Human1,'MAR03964_RPE')
#mod_Recon3D_Recon3D_results = FVA_FBA_analysis(mod_Recon3D_Recon3D,'DM_atp_c__RPE')

def save_df(df, filename):
        from datetime import datetime
        datestr = datetime.strftime(datetime.now(), '%H%M_%d-%m-%Y')
        filename = '' + datestr + '.xlsx'
        results_to_excel(df,Path().cwd() / 'results' / filename)

# post-processing (sorting) and saving of the results

def results_to_excel(df, path):
    # Write to Multiple Sheets
    with pd.ExcelWriter(path) as writer:
        df.to_excel(writer, sheet_name='full')
        df[df['subsystem'] == 'Exchange/demand reactions'].to_excel(writer, sheet_name='Exchange/demand reactions')

def create_permutation_dicts(rxn_id_bounds_dict):
    import itertools
    keys, values = zip(*rxn_id_bounds_dict.items())
    permutations_dicts = [dict(zip(keys, v)) for v in itertools.product(*values)]
    return permutations_dicts

