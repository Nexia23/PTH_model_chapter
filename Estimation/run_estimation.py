import sys
sys.path.append('../../Parameter_Sampler/')

from Estimator import ParameterEstimator
from estimation import *
from collections import OrderedDict
import time
import os
import json

def get_params_bounds(model_name):
    bounds = OrderedDict({
        #'Hkt_init': (0.35, 0.55, False),   
        'k_E_infect': (1e-8 , 4.5e-5, True),    # jetzt paramscan frueher (1e-8 , 1e-4),
        'tropism': (2, 30, False),
        'M':  (1e1, 5e3, True),
        'a_P_d': (15.6e1, 1e7, True),        
        'k_P_d':  (1e-3, 10, True),      
        'r_P_d': (1, 15, False),
        'fac_R_d': (1e-16, 1, True),
        'k_P_art_max': (1e-9, 1e1, True),
        't_mat_P': (5, 9, False),
        'slope_rpi':(1e0, 500, True),           # own idea
        #'k_M_death': (60, 100, False),
        #'t_E_death': (100, 130, False),
        'k_iE_pit_frac': (0, 1, False),         # Anteil der iE die durch ART gepitted werden, 0-1. Rest sterben durch ART
        #'k_iE_art_max': (1, 1000),             # maximale abtötrate von iE durch ART, #parameterscan zu unsensibel 
        #'h_art' :(3, 6, False),                # assume similar to in vitro: [1] R. K. Tyagi u. a., doi: 10.1186/s12916-018-1156-x.
        #'ID50': (1e-1, 1000),                  # ART dosis bei der 50% der iE getötet werden, #parameterscan zu unsensibel       
        #'pre_t': (2,6, True),                  # time of ART addition, 3 and 5 in medians in data for non-pth and pth respectively
        
        #### Pth specific parameteres
        #'Hkt_init_pth': (0.35, 0.55, False),   
        #'tropism_pth': (2, 200, False),
        't_E_death_pth': (100, 130, False),
        's_BH_pth': (1e-8, 1e-4, True),         # slope of linear function defining bystander heamolysis strength
        'LDH_pth': (140, 280, False),           # LDH concentration in blood plasma
        'k_M_death_pth': (30, 100, False),
        #'M_pth':  (1e2, 5e3, True),
        #'k_iE_pit_frac_pth': (0, 1, False),        # Anteil der iE die durch ART gepitted werden, 0-1. Rest sterben durch ART

        #### non-Pth specific parameteres
        #'Hkt_init_non': (0.35, 0.55, False),  
        #'tropism_non': (2, 200, False),
        't_E_death_non': (100, 130, False),
        's_BH_non': (1e-15, 1e-5, True),      # slope of linear function defining bystander heamolysis strength
        'LDH_non': (140, 280, False),       # LDH concentration in blood plasma
        'k_M_death_non': (30, 100, False),
        #'t_E_death_inf_non': (40, 130, False),
        #'M_non':  (1e2, 5e3, True),
        #'k_iE_pit_frac_non': (0, 1, False),        # Anteil der iE die durch ART gepitted werden, 0-1. Rest sterben durch ART
        })
    extra_bounds = {}
    if model_name == 'Hapto':
        #bounds.pop('s_BH_pth') # no pop as switch_fHb is very sensitive even tho just multiplied
        #bounds.pop('s_BH_non')
        extra_bounds = {
            't_halb_HP_decay' : (2, 6, False),         # 2-5 Tage zotero, 1Quellen: https://link.springer.com/chapter/10.1007/978-3-662-48986-4_1389 und weiterverfolgen https://archive.org/stream/WilliamsHematology9thEditionMcGrawHill_201805/Williams%20Hematology%2C%209th%20Edition%20McGraw-Hill_djvu.txt
            't_halb_HCC_decay': (3e-3, 9e-3, False), 
            #### Pth specific parameteres
            'par1_fHb_pth': (4e3, 6e3, False),
            'par2_fHb_pth': (4e-4, 4e-3, False),
            'switch_fHb_pth' : (1e-6, 8e3, True),
            #### non-Pth specific parameteres
            'par1_fHb_non': (1e3, 1e4, False),
            'par2_fHb_non': (1e-4, 1e-2, True),
            'switch_fHb_non' : (1e-6, 6e3, True),
            }
    elif model_name=="immune":
        # delete not used parameters from bounds
        bounds.pop('s_BH_pth')
        bounds.pop('s_BH_non')
        #bounds.pop('fac_R_d')

        extra_bounds = {
            'Treg':(1,1e3, True),
            'k_digest_R': (1e-15, 1e-8, True), 
            'k_digest_iE': (1e-7, 1e-4, True), 
            'delta_Treg':(1e-3, 5e2, True),       # 109
            'mu_tox':(1e-3,1e2, True),            # 22
            #'mu_in_tox':(1e0,1e8, True),         # 1e5
            'V_f':(1e-2, 5e2, True),              # 165
            'K_f':(1e0, 2e2, True),               # 20
            'delta_Ttox':(1e-3, 1e2, True),       # 0.6
            'epsilon':(1e-6, 2e0, True),          # 0.8
        
            #### Pth specific parameteres
            'beta_in_Treg_pth':(2e-4, 2e2, True),
            'mu_in_tox_pth':(1e0,1e8, True),      # 1e5

            #### non-Pth specific parameteres
            'beta_in_Treg_non':(2e-4, 2e2, True),
            'mu_in_tox_non':(1e0,1e8, True),      # 1e5
            }
    bounds.update(extra_bounds)
    return bounds

def calculate_cma_std(bounds):
    """Calculate Stds for params for CMA-ES from parameter space."""
    stds = []
    for p_key in bounds.keys():
        p = bounds[p_key]
        if p[2]:
            p_range = np.log10(p[1]) - np.log10(p[0])
        else:
            p_range = p[1] - p[0]
        stds.append(p_range / 4)
    print(f'Set CMA stds to {dict(zip(bounds.keys(), stds))}.')
    return stds
    

def save_estimation(best_score, best_parameters, update_params, ParamEster, fit_data:str, 
                    bounds: OrderedDict, run_id:int, location: str = 'general', patient: int=0):
    
    keys = update_params.keys()
    for key in update_params:
        best_pars = {k.replace("_"+key,""):v for k,v in best_parameters.items() if key in k or not k.split('_')[-1] in keys}
        best_pars.pop("pre_t")
        save_dict = update_params[key]
        save_dict.update(best_pars)
        best_pars["pre_t"] = best_parameters["pre_t"]
       
        result_dict = {'fitted_data': fit_data, 'parameter_space': bounds, 'run_id': run_id,
                       'best_score': best_score, 'best_parameters': best_pars, 
                       'update_parameters':save_dict,
                       'cost_history': ParamEster.cost_history,
                       'function_calls': ParamEster.function_calls,
                       'full_cost_history': ParamEster.complete_cost_history}

        t = time.strftime("%Y_%m_%d_%H_%M")
        name = t

        save_to = key
        if key == "non": save_to ="non_pth"
        if not patient: folder_loc=f'{location}/{save_to}'
        else: folder_loc=f'{location}/{save_to}/patient{patient}'
        os.makedirs(folder_loc, exist_ok=True)

        with open(folder_loc+f'/{name}_{run_id}.json', "w") as write_file:
            json.dump(result_dict, write_file, indent=4)


def extract_patient_fitting_data(df: pd.DataFrame, features: list=['Hb', 'LDH'], patient_num: int=20):
 
    patient_df = df[df['patientnumber']==patient_num]
    fitting_df = patient_df[['time'] + features]
    # Change naming to fit estimation.py cost function
    features_fit = [x+"_mean" for x in features]
    fitting_df.columns = ['Time'] + features_fit

    return fitting_df


def main():
    model_name = sys.argv[1]
    run_id = sys.argv[2]
    patient_num = 0
    if len(sys.argv) == 4:
        patient_num  = sys.argv[3]
        df = pd.read_csv('../datasets/OIE_data.csv')
        dpatient_fitting = extract_patient_fitting_data(df, features=['Hb', 'LDH', '[R]'], patient_num=patient_num)
        data = {"pth":dpatient_fitting,       
                "non":dpatient_fitting}
        data_used =[f"patient{patient_num}",
                    f"patient{patient_num}"]
    else:
        data = {"pth":pd.read_csv("pth.csv"),       
                "non":pd.read_csv("non_pth.csv")}
        data_used =["pth.csv",
                    "non_pth.csv"]

    model = te.loada(f'../LCT_model/{model_name}_LCT_OIE.ant')
    pop =  30
    print(data)
    if model_name=='immune':
        pop = 500
    est_obj = FitManager(model, data, model_name) 
    ParamEster = ParameterEstimator()
    bounds = get_params_bounds(model_name)
    stds = calculate_cma_std(bounds)
    ParamEster.initialize(est_obj.objective_function, bounds)
    best_score, best_parameters, runtime = ParamEster.run(method='cma', iterations=20,
                                                          run_id=run_id, n_lhs=1,
                                                          optimizer_args={'CMA_stds': stds,
                                                                          'popsize': pop,})
    # Saving of estimation results
    if 'pre_t' in best_parameters:  
        pre_t = best_parameters.pop('pre_t')
    else:
        pre_t = est_obj.default_pre_t
    
    update_parameters={} 
    keys = ['pth','non']
    for key in data.keys():
        usedpars = {k.replace("_"+key,""):v for k,v in best_parameters.items() if key in k or not k.split('_')[-1] in keys}
        dummy_dict = get_steady_state(model, usedpars, model_name)
        update_parameters[key] = {k:v for k,v in dummy_dict.items() }

    best_parameters['pre_t'] = pre_t
    save_estimation(best_score, best_parameters, update_parameters, ParamEster,
                    data_used, bounds, run_id, location=model_name, patient=patient_num)

    print(best_score, best_parameters, runtime)

if __name__=='__main__':
    main()