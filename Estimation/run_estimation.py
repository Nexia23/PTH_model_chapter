import sys
sys.path.append('../../Parameter_Sampler/')

from Estimator import ParameterEstimator
from estimation import *
from collections import OrderedDict
import time
import os
import json

def get_params_bounds():
    bounds = OrderedDict({
    #'Hkt_init': (0.35, 0.55, False),   
    'k_E_infect': (7e-7 , 4.5e-5, True),    # jetzt paramscan früher (1e-8 , 1e-4),
    'tropism': (2, 200, False),
    'M':  (1e2, 5e3, True),
    'a_P_d': (15.6e1, 1e7, True),        
    'k_P_d':  (1e-3, 10, True),      
    'r_P_d': (1, 15, False),
    'fac_R_d': (1e-16, 1, False),
    'k_P_art_max': (1e-9, 1e1, True),
    't_mat_P': (5, 14, False),
    'slope_rpi':(1e0, 500, True),            # own idea
    #'k_M_death': (60, 100, False),
    #'t_E_death': (100, 130, False),
    'k_iE_pit_frac': (0, 1, False),         # Anteil der iE die durch ART gepitted werden, 0-1. Rest sterben durch ART
    #'k_iE_art_max': (1, 1000),             # maximale abtötrate von iE durch ART, #parameterscan zu unsensibel 
    #'h_art' :(3, 6, False),                # assume similar to in vitro: [1] R. K. Tyagi u. a., doi: 10.1186/s12916-018-1156-x.
    #'ID50': (1e-1, 1000),                  # ART dosis bei der 50% der iE getötet werden, #parameterscan zu unsensibel       
    #'pre_t': (2,6, True),                  # time of ART addition, 3 and 5 in medians in data for non-pth and pth respectively
    
    # Pth specific parameteres

    #'Hkt_init_pth': (0.35, 0.55, False),   
    #'t_mat_P_pth': (5,14, False),
    #'tropism_pth': (2, 200, False),
    #'k_P_art_max_pth': (1e-9, 1e1, True),
    't_E_death_pth': (100, 130, False),
    's_BH_pth': (1e-9, 1e-5, True),         # slope of linear function defining bystander heamolysis strength
    'LDH_pth': (140, 280, False),           # LDH concentration in blood plasma
    'k_M_death_pth': (30, 100, False),
    #'M_pth':  (1e2, 5e3, True),
    #'a_P_d_pth': (15.6e1, 1e7, True),        
    #'k_P_d_pth':  (1e-3, 10, True),      
    #'r_P_d_pth': (4, 15, False),
    #'fac_R_d_pth': (1e-16, 1, False),
    #'k_iE_pit_frac_pth': (0, 1, False),        # Anteil der iE die durch ART gepitted werden, 0-1. Rest sterben durch ART

    # non-Pth specific parameteres

    #'Hkt_init_non': (0.35, 0.55, False),  
    #'t_mat_P_non': (5,14, False),
    #'tropism_non': (2, 200, False),
    #'k_P_art_max_non': (1e-9, 1e1, True),
    't_E_death_non': (100, 130, False),
    's_BH_non': (1e-11, 1e-8, True),      # slope of linear function defining bystander heamolysis strength
    'LDH_non': (140, 280, False),       # LDH concentration in blood plasma
    'k_M_death_non': (30, 100, False),
    #'t_E_death_inf_non': (40, 130, False),
    #'M_non':  (1e2, 5e3, True),
    #'a_P_d_non': (15.6e1, 1e7, True),        
    #'k_P_d_non':  (1e-3, 10, True),      
    #'r_P_d_non': (4, 15, False),
    #'fac_R_d_non': (1e-16, 1, False),
    #'k_iE_pit_frac_non': (0, 1, False),        # Anteil der iE die durch ART gepitted werden, 0-1. Rest sterben durch ART
    })
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
                    bounds: OrderedDict, run_id:int):
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
        
        os.makedirs(f'{save_to}', exist_ok=True)

        with open(f'{save_to}/{name}_{run_id}.json', "w") as write_file:
            json.dump(result_dict, write_file, indent=4)


def main():
    # fit_data = sys.argv[1]
    run_id = sys.argv[1]
    model = te.loada('../LCT_model/LCT_OIE.ant')

    #data = pd.read_csv(fit_data)
    data = {"pth":pd.read_csv("pth.csv"),       
            "non":pd.read_csv("non_pth.csv")}
    data_used =["pth.csv",
                "non_pth.csv"]

    est_obj = FitManager(model, data)    
    ParamEster = ParameterEstimator()
    bounds = get_params_bounds()
    stds = calculate_cma_std(bounds)
    ParamEster.initialize(est_obj.objective_function, bounds)
    best_score, best_parameters, runtime = ParamEster.run(method='cma', iterations=20,
                                                          run_id=run_id, n_lhs=1,
                                                          optimizer_args={'CMA_stds': stds})
    # Saving of estimation results
    if 'pre_t' in best_parameters:  
        pre_t = best_parameters.pop('pre_t')
    else:
        pre_t = est_obj.default_pre_t
    
    update_parameters={} 
    keys = ['pth','non']
    for key in data.keys():
        usedpars = {k.replace("_"+key,""):v for k,v in best_parameters.items() if key in k or not k.split('_')[-1] in keys}
        dummy_dict = get_steady_state(model, usedpars)
        update_parameters[key] = {k:v for k,v in dummy_dict.items() }

    best_parameters['pre_t'] = pre_t
    save_estimation(best_score, best_parameters, update_parameters, ParamEster,
                    data_used, bounds, run_id)

    print(best_score, best_parameters, runtime)

if __name__=='__main__':
    main()