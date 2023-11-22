#packages

# paths
import sys
sys.path.append('../datasets') 
import dataset_long as dsl


# Model
import numpy as np
import tellurium as te
# Plotting
import matplotlib.pyplot as plt
import seaborn as sns
# Dataframe
import pandas as pd
# Estimation
from scipy.optimize import minimize 
from collections import OrderedDict
from scipy.stats import qmc
from scipy.optimize import differential_evolution
# Analyse Estimation
from scipy.stats import chi2
from scipy.stats import norm

# Save optimization params
import json
import time

#debugger
import logging as log
logger = log.getLogger()
logger.setLevel(log.DEBUG)  #log.debug anstatt print, INFO (printed nicht), DEBUG(printed)


def get_lab_values():
    lab_values = [
        # Define lab values of interest
        #'parasitemia',
        'Hkt', 
        '[E]', 
        'Hb',
        'LDH',
        '[R]',
        #'R_percent',
        #'oiE_percent', 
        '[oiE]',
        #'[iE]'                                 #gerade nihct fitten
        #'Pittingquote_absolut_modified',
        #'Pittingquote_modified',
        ]
    return lab_values

def get_params_guess():
    guess = OrderedDict({
    'Hkt_init': 0.45,
    'M': 40e3,
    #'k_R_death': 0,   # aus Parameter scan kein großen einfluss auf alles
    'k_E_infect': 1e-6 ,
    'k_iE_rupture':1,
    'k_M_death': 48 ,

    'param1_Pdeath': 0.00071535, #sigmoid 0.5 ,       #parameterscan zu unsensibel
    'param2_Pdeath': 0.48924947,  #sigmoid 1  ,    #parameterscan zu unsensibel
    #'param3_Pdeath': 14.5,

   # 'I0_iE':    0, 
    'I0_death_iE': 0,    
    #'Imax_iE': 10,    #parameterscan zu unsensibel
    #'hill' :2.0,      
    #'ID50': 20,       #parameterscan zu unsensibel

    't_oiE_death': 10,   
    #'t_halb_LDH_decay': 10, #4,  rauswerfen, da grenzen zu klein um wirlich verändeurng zu haben


    'k_BH': 5e-4,
    'J_oiEdeath_0': 7500,

    'BH_max': 1,
    't_ART_add': 7,
    #'k_iE_kill_proportion': 0.5,  # wenn angeschlate nur scores über 600

    'param1_TRPaging':  3.53276388,
    'par#am2_TRPaging': 5.99745537,  
    'param3_TRPaging': 0.29658879,
    
    })

    return guess


def get_params_bounds():
    bounds = OrderedDict({
    'Hkt_init': (0.35, 0.55),   
    'M': (1e2,1e5),                     # discard for fitting mean
    #'k_R_death': (1e-5, 1e-4),         # k_R_death < k_R_aging=1.38, drei einheiten kleiner setzen
    'k_E_infect': (1e-6 , 2.5e-6),      # jetzt paramscan früher (1e-8 , 1e-4),
    'k_iE_rupture':(1e-1, 10),          # lebenszeit 2 Tage, Half-life 1 day, rate ln(2)
    'k_M_death': (0,1e3) ,              # k_M_death ~ 48 (lebenszeit 30 min), drei einheiten kleiner setzen

    'param1_Pdeath': (1e-6, 1e-1),      # sigmoid(1e-2, 1e2) ,      #linear #parameterscan zu unsensibel
    'param2_Pdeath':  (1e-6, 1),        # sigmoid(1e-2,1e3) ,     #parameterscan zu unsensibel
    #'param3_Pdeath': (5, 18),          # nur bei sigmoid notwendig

    #'I0_iE':    (0, 1e-3),             # vlt auf 0 setzen ausprobieren, Inhibition ohne ART, zw. 0 und viel kleiner als Imax_iE 
    'I0_death_iE': (0, 1e-3), 
    #'Imax_iE': (1, 1000),              # maximale abtötrate von iE durch ART, #parameterscan zu unsensibel 
    #'hill' :(1, 4),                    # Anstieg, sehr sensibey= 
    #'ID50': (1e-1, 1000),              # ART dosis bei der 50% der iE getötet werden, #parameterscan zu unsensibel       

    # 't_oiE_death': (5,15),            # Zeit nach der oiE sterben, nach ca. 7- 14
    # 't_halb_LDH_decay': (2,14),       # Halbwertzeit von LDH, 3-5 Tage+ Puffer


    'k_BH': (5e-5, 5e-3),               # Anstieg Kurve 
    'J_oiEdeath_0': (5000, 10000),      # in Modell max. 15.000
    

    'BH_max': (2e-1, 20),  
    't_ART_add': (3,20),
    #'k_iE_kill_proportion': (0,1),  # wenn angeschlate nur scores über 600

    'param1_TRPaging':(1, 10),   #  maximal 11             = 3.53276388
    'param2_TRPaging': (1e-1, 20) ,                   #5.99745537
    'param3_TRPaging': (0.15, 0.60) #0.29658879'


    })

    return bounds


def preprocess_experimental_df(df, patient, lab_values):
    #patient = 'patientnumber'
    df_pat = df[df['patientnumber'] == patient]
    df_experiment_data = df_pat.set_index('time')[lab_values]
    return df_experiment_data

def update_model(model, params=None):
    if params is not None:
        if not isinstance(params, OrderedDict):
                raise ValueError("Params must be an ordered dictionary")
            
        model_str = model.getAntimony()        #lade Antimony string
        model_lines = model_str.split('\n')    #splitte in für jede Zeile

        """
        new_model_lines = []
        for name, val in params.items():       # gehe alle zeielen durch, ändere wnn parameter geändert werden soll, füge in neuen strig zu
            for line in model_lines:
                if line.startswith(f'{name}=') or line.startswith(f'{name} '): 
                    line = f'{name} = {10**val}'
                    print(line)
                new_model_lines.append(line)
        new_model_str = '\n'.join(new_model_lines)  
        print(len(new_model_str), len(model_str))  
        """

        new_model_lines = []
        for line in model_lines:
            for name, val in params.items():       # gehe alle zeielen durch, ändere wnn parameter geändert werden soll, füge in neuen strig zu
                if line.startswith(f'{name}=') or line.startswith(f'{name} '): 
                    line = f'{name} = {10**val}'
                    break
            new_model_lines.append(line)
        new_model_str = '\n'.join(new_model_lines)  
        model = te.loada(new_model_str)         #lade überarbeiteten string

        return model


def simulate_model(ant_model, lab_values, params = None):
    # Simulates the model and returns the lab values at different timepoints

    # Load the model
    model= te.loada(ant_model)

    # Create a list of time-points and uncertainties
    time_uncertainty = OrderedDict({0: 0, 2: 0, 6: 2, 13: 3, 27: 4})
    time_experiment = list(time_uncertainty.keys())
    times_experiment_uncertainty =list(time_uncertainty.values())
    #times_model = [int(timepoint + model.t_ART_add) for timepoint in time_experiment]

    if 't_ART_add' in params:
        t_ART_add_new = params['t_ART_add']
        times_model = [int(timepoint + t_ART_add_new) for timepoint in time_experiment]
    else: times_model = [int(timepoint + model.t_ART_add) for timepoint in time_experiment]

    model = update_model(model, params)

    # Simulate model, and store lab values at the time points of medication
    max_time = times_model[-1]+ times_experiment_uncertainty[-1]
    step_size = max_time+1
    result = model.simulate(0, max_time, step_size, selections = ['time']+lab_values )
    simulation_result = pd.DataFrame(result, columns= result.colnames)
    df_simulated_data= pd.DataFrame(columns=lab_values)

    # Calculate mean of lab values at timepoints with uncertainty  
    for i, timepoint in enumerate(times_model):
        uncertain = times_experiment_uncertainty[i]
        min = timepoint - uncertain
        max = timepoint + uncertain
        df_simulated_data.loc[timepoint] = simulation_result.loc[min:max,].mean()
    df_simulated_data.index = time_experiment

    return df_simulated_data


def calculate_objectiv_score(df_simulated_data, df_experimental_data):
    lab_values= get_lab_values()
    objectiv_function = 0

    for lab_value in lab_values:
        model_values = df_simulated_data[lab_value]
        experimental_values = df_experimental_data[lab_value]
        #log.debug(model_values, experimental_values)
        objectiv_function += np.nansum((((experimental_values - model_values)/experimental_values.mean())**2)) #objectiv_function += (((experimental_value - model_value)**2)/experimental_std**2).sum()
        #log.debug(lab_value, objectiv_function)
    return objectiv_function


def objectiv_function (param_values, df_experimental_data, ant_model, param_names):
    #log.debug(param_values)
    params = OrderedDict(dict(zip(param_names, param_values)))
    df_simulated_data = simulate_model(ant_model, get_lab_values(),  params)
    calculated_error = calculate_objectiv_score(df_simulated_data, df_experimental_data)
    log.debug(calculated_error)

    return calculated_error


def simulated_data_for_plot(ant_model, lab_values, params = None):
       # Simulates the model and returns the lab values at different timepoints

    # Load the model
    model= te.loada(ant_model)

    # Create a list of time-points and uncertainties
    time_uncertainty = OrderedDict({0: 0, 2: 0, 6: 2, 13: 3, 27: 4})
    time_experiment = list(time_uncertainty.keys())
    times_experiment_uncertainty =list(time_uncertainty.values())
    times_model = [int(timepoint + model.t_ART_add) for timepoint in time_experiment]

    # Reset all parameters given as arguments
    if params is not None:
        if not isinstance(params, OrderedDict):
            raise ValueError("Params must be an ordered dictionary")
        for name, val in params.items():
            #model.setValue(name, val)
            model.setValue(name, 10**val)
    
    # Simulate model, and store lab values at the time points of medication
    max_time = times_model[-1]+ times_experiment_uncertainty[-1]
    step_size = max_time+1
    result = model.simulate(0, max_time, step_size, selections = ['time']+lab_values )
    df_simulated_data= pd.DataFrame(result, columns= result.colnames)

    df_simulated_data['time'] -= model.t_ART_add  # subtract model.t_ART_add from the 'time' column
    df = df_simulated_data.loc[df_simulated_data['time'] >= 0]  # select only rows where time >= 0
    df = df.drop(columns=['time']).reset_index(drop=True)
    return df


def plot_simulation_and_experiment(df_experimental_data, df_simulated_data, figsize=(10, 6)):
    n_subplots = len(get_lab_values())
    n_cols = 3
    n_rows = int(np.ceil(n_subplots/n_cols))

    fig, ax = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=figsize)

    for i, lab_value in enumerate(get_lab_values()):
        col = i%n_cols
        row = i//n_cols
        ax[row, col].plot(df_simulated_data.time, df_simulated_data[lab_value])
        ax[row, col].plot(df_experimental_data.time, df_experimental_data[lab_value],'o')
        ax[row, col].set_title(lab_value)
    plt.tight_layout()
    return fig
   # plt.show()


def convert_bounds_to_logscale(bounds):
    """
    converts bounds for parameters from linear scale to log10 scale. If a bound is 0, it's replaced with 1e-8 to avoid infinity values in log scale.
    Input and output: Ordered dictionaries
    """
    log_bounds = OrderedDict({})

    for param, bound in bounds.items():
        if bound[0] == 0:
            bound = (bound[1]*1e-3, bound[1])
        log_bounds[param] = (np.log10(bound[0]), np.log10(bound[1]))

    return log_bounds


def generate_LHS_sampling(bounds, n_samples):
    """
    Generate LHS samples in giving scale
    Input:
        bounds = OrderedDict({'param2_EPOprod': (0,300),...})
        n_samples = int, number of samples to generate
    Output: 
        lhs_sample = list of lists, each list is a generated sample [[guess1],[guess 2],..]  bsp [[0,4,6],[1,3,5],..]
    """
    sample_dim = len(bounds)                     # define dimension of sample = number of parameters 
    sampler = qmc.LatinHypercube(d=sample_dim)   # Create the sampler object
    lhs_sample = sampler.random(n=n_samples)     # Generate n_samples(=number) of LHS samples

    lower_bounds = [bounds[param][0] for param in bounds]
    upper_bounds = [bounds[param][1] for param in bounds]
    lhs_sample = qmc.scale(lhs_sample, lower_bounds, upper_bounds)  #scale LHS sample to bounds of params

    return lhs_sample


def create_intial_guess(bounds,sample):
    """Creates an ordered dictionary of initial guess using LHS."""
    intial_guess = OrderedDict({})

    for i, param in enumerate(bounds.keys()):
        intial_guess[param] = sample[i]
    
    return intial_guess


def convert_params_bounds_to_normalscale(log_bounds):
    """ Convert the samples back to the original scale"""
    bounds= OrderedDict({})
    for i, (param, bound) in enumerate(log_bounds.items()):
        bounds[param] = 10 ** np.array(bound)

    return bounds




if __name__ =='__main__':
    # Differential Evolution Algorithm
    start = time.time()
    timestr = time.strftime("%Y%m%d-%H-%M-%S")
    ID = 1
    patient_id = '36'
    ID =  sys.argv[1]
    patient_id = str(sys.argv[2])

    # parameter bounds
    params_bounds = get_params_bounds()
    log_bounds = convert_bounds_to_logscale(params_bounds)

    # LHS
    num_samples = 1                                      # number of LHS samples to generate
    samples = generate_LHS_sampling(log_bounds, num_samples)    # best_samples 

    # Model
    ant_model = '.././model/OIE_model.ant'

    # experiment data
    patient = patient_id
    data_df = pd.read_excel('.././datasets/haemolysismodel_conRetis.xlsx')   
    df = dsl.long_format(data_df)
    df_experimental_data = preprocess_experimental_df(df, patient, get_lab_values())

    optimization_results = []

    # Objectiv function
    for index, sample in enumerate(samples):
        intial_guess = create_intial_guess(params_bounds, sample)

        try:
            start = time.time()
            result = differential_evolution(
                objectiv_function,
                x0=list( intial_guess.values()),
                args=(df_experimental_data, ant_model, intial_guess.keys()),
                bounds=list(log_bounds.values()),
                maxiter=1000
            )
           
        
            results_dict = {
                    'index': index,
                    'initial_guess': intial_guess,
                    'best_score': result.fun,
                    'best_params': dict(zip(list(intial_guess.keys()), result.x.tolist())),
                    'termination_cause': result.message,
                    
                    #'covariance_matrix': result.hess_inv.todense().tolist()
                }
            optimization_results.append(results_dict)

            #with open(timestr + '_paras_optimized_sample_' + str(index) + '.json', 'w') as file:       #einzeln speichern
            #    json.dump( optimization_results, file, indent=4)


        except RuntimeError as err:    #für debugging später
            results_dict = {
                'index': index,
                'initial_guess': intial_guess,
                'error': str(err),
            }
            optimization_results.append(results_dict)

    end = time.time()
    
    output= {
        'metadata': {
            'description': 'k_P_death linear, ohne iE, 1x durchlauf,maxiter 1000, mit fixed update function, ohne hill, paramTRPaging, paramPdeaths',
            'algorithm':'differential_evolution', 
            'patient': patient, 
            'estimation_duration': end - start,
        },
        'results': optimization_results
    }

    with open( timestr + f'paras_optimized_{ID}_{patient_id}.json', 'w') as file:
        json.dump(output, file, indent=4)

