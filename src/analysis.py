# perform FVA and FBA on a model setting RPE_ATP demand to fixed values

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

