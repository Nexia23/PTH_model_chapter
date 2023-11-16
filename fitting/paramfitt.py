#packages

#paths
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

#Estimation
from scipy.optimize import minimize 
from collections import OrderedDict
from scipy.stats import qmc
#Analyse Estimation
from scipy.stats import chi2
from scipy.stats import norm

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
        '[iE]' 
        #'Pittingquote_absolut_modified',
        #'Pittingquote_modified',
        ]
    return lab_values


def preprocess_experimental_df(df, filter_by_PTH, lab_values):
    #if filter_by_PTH = 1, only cases with PTH,   if 0, only cases without PTH 
    df_filtered = df.loc[df['PTH'] == filter_by_PTH ]                                                 
    df_experiment_data= df_filtered.groupby('time')[lab_values].agg(['mean', 'std', 'sem'])              # group the dataframe df by the 'time' and calculate  mean, standard deviation and SEM of each lab-value
    df_experiment_data.columns = ['{}_{}'.format(col, stat) for col, stat in df_experiment_data.columns]  # rename columns of the resulting df to 'Hb_mean', 'Hb_std', ...

    return df_experiment_data


def simulate_model(ant_model, lab_values, params = None):
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
        model_value = df_simulated_data[lab_value]
        experimental_value = df_experimental_data[f"{lab_value}_mean"]
        experimental_std = df_experimental_data[f"{lab_value}_std"]
        objectiv_function += (((experimental_value - model_value)**2)/experimental_std**2).sum()

    return objectiv_function


def objectiv_function (param_values, df_experimental_data, ant_model, param_names):
    params = OrderedDict(dict(zip(param_names, param_values)))

    if 'k_R_death' in params:
        k_R_death = params['k_R_death'] 
        R0 = 46696.30405453991
        Hkt = 0.45
        param1_TRPaging= 3.53276388
        param2_TRPaging=-5.99745537
        param3_TRPaging= 0.29658879
        t_R_aging = param1_TRPaging/ (1+np.exp(-param2_TRPaging*(Hkt - param3_TRPaging)))
        k_R_aging = np.log(2) / (t_R_aging/2)
        t_P_aging = 11-(param1_TRPaging/ (1+np.exp(-param2_TRPaging*(Hkt - param3_TRPaging))))
        k_P_aging= np.log(2) / (t_P_aging/2)   

        P0= (R0 * k_R_death + R0 * k_R_aging) / (k_P_aging * 2**10)
        params['P'] = P0

    df_simulated_data = simulate_model(ant_model, get_lab_values(),  params)
    calculated_error = calculate_objectiv_score(df_simulated_data, df_experimental_data)
    print(calculated_error)

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
    n_subplots = len(df_simulated_data.columns)
    n_cols = 3
    n_rows = int(np.ceil(n_subplots/n_cols))

    fig, ax = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=figsize)

    for i, lab_value in enumerate(df_simulated_data.columns):
        col = i%n_cols
        row = i//n_cols
        ax[row, col].plot(df_simulated_data.index, df_simulated_data[lab_value])
        ax[row, col].errorbar(df_experimental_data.index, df_experimental_data[f"{lab_value}_mean"], yerr=df_experimental_data[f"{lab_value}_std"], fmt='o')
        ax[row, col].set_title(lab_value)
    plt.tight_layout()
    plt.show()


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
    #Parameter
    params = OrderedDict({'param2_EPOprod': 100, 'param3_EPOprod': 0.04})

    #Model
    ant_model = 'OIE_model.ant'

    #Data von Pinkus 
    data_df = pd.read_excel('./haemolysismodel_conRetis.xlsx')   
    df = dsl.long_format(data_df)
    df_experimental_data = preprocess_experimental_df(df, 1, get_lab_values())

    #Residual Error
    result_objectiv_score = objectiv_function(params.values(), df_experimental_data, ant_model, params.keys())
    print(result_objectiv_score)


"""
def objectiv_function(param_values, df_experimental_data, ant_model, param_names):

    Calculates the residual error ofbetween the experimental data and simulated data from the model, given the parameter values.

    Parameters:
    params: a list of tuples where each tuple contains the name of the parameter to change
            and its new value, e.g. [('param2_EPOprod', 100)], [('param3_EPOprod',0.04)]
    df_experimental_data: a pandas df containing the experimental data for the lab values of interest
    ant_model: Antimony model as a string

    Returns:
    objectiv_function: a float representing the residual error between the experimental and simulated data

    lab_values = [
        'parasitemia',
        #'PfHRP2', 
        'Hkt', 
        '[E]', 
        'Hb',
        'LDH',
        'RPI',
        'oiE_percent', 
        '[oiE]',
        '[iE]' 
        #'Pittingquote_absolut_modified',
        #'Pittingquote_modified',
        ]

    # Load the model
    model= te.loada(ant_model)

    # Create a list of time points from the experimental data
    time_experiment = list(df_experimental_data.index)
    times_uncertainty = [0,0,2,3,4]
    times_model = [int(timepoint + model.t_ART_add) for timepoint in list(df_experimental_data.index)]

    # make params dict
    params = dict(zip(param_names, param_values))

    # Reset all parameters given as arguments
    for name, val in params.items():
        model.setValue(name, val)
        #model.setValue(name, 10**val)
    
    # Simulate model, and store lab values at the time points of medication
    result = model.simulate(0, times_model[-1] , times_model[-1]+1, selections = ['time']+lab_values )
    simulation_result = pd.DataFrame(result, columns= result.colnames)
    df_simulated_data= pd.DataFrame(columns=lab_values)

    #mean von lab_values an zeitpunkten mit unsicherheit zeitpunkt
    for i, timepoint in enumerate(times_model):
        uncertain = times_uncertainty[i]
        min = timepoint - uncertain
        max = timepoint + uncertain
        df_simulated_data.loc[timepoint] = simulation_result.loc[min:max,].mean()

    df_simulated_data.index = time_experiment

    # Compare the experimental data with the simulated data from the model
    objectiv_function = 0
    for lab_value in lab_values:
        model_value = df_simulated_data[lab_value]
        experimental_value = df_experimental_data[f"{lab_value}_mean"]
        experimental_std = df_experimental_data[f"{lab_value}_std"]
        objectiv_function += (((experimental_value - model_value)**2)/experimental_std**2).sum()

    return objectiv_function


#Example Usage
params = OrderedDict({'param2_EPOprod': 100, 'param3_EPOprod': 0.04})
df_experimental_data = df_input_PTH
ant_model = 'OIE_model.ant'
objectiv_function(params.values(), df_experimental_data, ant_model, params.keys())

"""