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
   # 'Hkt_init': (0.35, 0.55, False),   
    'k_E_infect': (1e-7 , 2.5e-6, True),    #jetzt paramscan früher (1e-8 , 1e-4),

   # 's_P_d': (1e-6, 1e-1, True),     # sigmoid(1e-2, 1e2) ,      #linear #parameterscan zu unsensibel
  #  'k_P_d0':  (1e-6, 1, True),      # sigmoid(1e-2,1e3) ,     #parameterscan zu unsensibel


    #'k_iE_pit_frac': (0, 1, False),            # Anteil der iE die durch ART gepitted werden, 0-1. Rest sterben durch ART
    #'k_iE_art_max': (1, 1000),             # maximale abtötrate von iE durch ART, #parameterscan zu unsensibel 
   # 'h_art' :(3, 6, False),                 # assume similar to in vitro: [1] R. K. Tyagi u. a., doi: 10.1186/s12916-018-1156-x.
    #'ID50': (1e-1, 1000),               #ART dosis bei der 50% der iE getötet werden, #parameterscan zu unsensibel       

    's_BH': (1e-9, 1e-6, True),      # slope of linear function defining bystander heamolysis strength
 
    #'pre_t': (2,6, True),            # time of ART addition, 3 and 5 in medians in data for non-pth and pth respectively
    'LDH': (140, 280, False),       # LDH concentration in blood plasma
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
    

def save_estimation(best_score, best_parameters, update_params,ParamEster, fit_data:str, bounds: OrderedDict, run_id:int):
    result_dict = {'fitted_data': fit_data, 'parameter_space': bounds, 'run_id': run_id,
            'best_score': best_score, 'best_parameters': best_parameters, 'update_parameters': update_params,
            'cost_history': ParamEster.cost_history,
            'function_calls': ParamEster.function_calls,
            'full_cost_history': ParamEster.complete_cost_history}


    t = time.strftime("%Y_%m_%d_%H_%M")
    name = t
    save_to = fit_data[:-4]
    os.makedirs(f'{save_to}', exist_ok=True)

    with open(f'{save_to}/{name}_{run_id}.json', "w") as write_file:
        json.dump(result_dict, write_file, indent=4)


def main():
    fit_data = sys.argv[1]
    run_id = sys.argv[2]
    model = te.loada('../LCT_model/LCT_OIE.ant')
    
    data = pd.read_csv(fit_data)
    est_obj = FitManager(model, data)
    ParamEster = ParameterEstimator()
    bounds = get_params_bounds()
    stds = calculate_cma_std(bounds)
    ParamEster.initialize(est_obj.objective_function, bounds)
    best_score, best_parameters, runtime = ParamEster.run(method='cma', iterations=20, run_id=run_id, n_lhs=1,optimizer_args={'CMA_stds': stds})
    if 'pre_t' in best_parameters:  
        pre_t = best_parameters.pop('pre_t')
    else:
        pre_t = est_obj.default_pre_t
    update_parameters = get_steady_state(model, best_parameters)
    best_parameters['pre_t'] = pre_t
    save_estimation(best_score, best_parameters, update_parameters,ParamEster, fit_data, bounds, run_id)


    print(best_score, best_parameters, runtime)

if __name__=='__main__':
    main()