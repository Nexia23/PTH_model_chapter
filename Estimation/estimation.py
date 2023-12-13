import sys
import numpy as np
import pandas as pd
import tellurium as te


def get_steady_state(model, pars: dict):
    for p in pars:
        model.setValue(p, pars[p])
    
    # setting Hkt defeins aging times of precursors (P) and reticulocytes (R)      
    t_R_a_init = model.t_R_a_max/ (1+np.exp (model.s_R_a*(model.Hkt_init - model.Hkt_0)))     # t_R_a_max, s_R_a, Hkt_0 stem from fitting, see plot_dependencies.ipynb (RPI)
    t_P_a_init = 11-(model.t_R_a_max/ (1+np.exp (model.s_R_a*(model.Hkt_init - model.Hkt_0)))) 

    # E_init determined by Hkt_init, t_R_a_init and t_E_death (fixed)
    E_init = (model.Hkt_init * model.Vol_blood) / (model.Vol_E + (t_R_a_init/model.t_E_death)*model.Vol_R)
    # R_init from E_init , t_R_a_init and t_E_death (fixed)
    R_init = E_init * t_R_a_init/ model.t_E_death
    # P_init from R_init, t_R_a_init and t_P_a_init
    P_init = (R_init * (model.k_R_death + 2*np.log(2)/t_R_a_init))/ (2**10 *2*np.log(2) / t_P_a_init)    

    # Hb_init, needed for J_P_death
    Hb_init = (model.Vol_E * E_init * model.Hb_conc_E + model.Vol_R * R_init * model.Hb_conc_R) / (10*model.Vol_blood)  # conc unit from g/l to g/dl, thus div. by 10

    # equilibrating precursors (P)
    J_P_death_init = P_init * (model.s_P_d * Hb_init + model.k_P_d0 )
    J_P_aging_init = P_init * np.log(2) / (t_P_a_init/2)
    k_P_birth_init   = J_P_death_init + J_P_aging_init

    # equilibrating LDH
    J_LDH_decay_init = model.LDH * (np.log(2) / model.t_halb_LDH_decay)
    J_E_death_init = E_init * 2*np.log(2) / model.t_E_death
    J_R_death_init = R_init * model.k_R_death
    LDH_RBC_init = (J_LDH_decay_init * model.Vol_blood) / (J_E_death_init + J_R_death_init ) 

    # setting initial values/params
    eq_dict = {}
    eq_dict['E'] = E_init
    eq_dict['R'] = R_init
    eq_dict['P'] = P_init
    eq_dict['k_P_birth'] = k_P_birth_init
    eq_dict['LDH_RBC'] = LDH_RBC_init

    # include chosen parameters
    eq_dict.update(pars)

    return eq_dict

def set_model_to_ss(model, pars: dict):
    update_pars = get_steady_state(model, pars)
    for p in update_pars:
        model.setValue(p, update_pars[p])    # CAUTION: if you include a volume which is not 1 in the model, species might get fucked up
   # print(update_pars)
    return model



class FitManager():

    def __init__(self, model, data ) -> None:
        self.model = model
        self.data = data
        self.default_pre_t = 3


    def objective_function(self, **pars):
        procID = pars['process_id']

        del pars['process_id']
        if 'pre_t' in pars:
            pre_t = pars['pre_t']
            del pars['pre_t']
        else:
            pre_t = self.default_pre_t
        
        # from here loop for both models
        error_sum = 0
        for key in self.data:
            usedpars = {k.replace("_"+key,""):v for k,v in pars.items() if key in k}
            usedpars['k_E_infect'] = pars["k_E_infect"]

            # set model to steady state
            self.model = set_model_to_ss(self.model, usedpars)
            
            # simulate 
            t_max = self.data[key]['Time'].max()
            res = self.model.simulate(-pre_t, int(t_max), int(t_max + pre_t + 1), selections=['time', 'Hb', '[LDH]'])
            res_df = pd.DataFrame(res, columns=res.colnames)
            
            # only keep timepoints which are in data
            res_df = res_df[res_df['time'].isin(self.data[key]['Time'])]
            # TODO: loop over data columns?
            Hb_error = ((self.data[key]['Hb_mean'].values - res_df['Hb'].values) / (self.data[key]['Hb_mean'].max() - self.data[key]['Hb_mean'].min()))**2
            LDH_error = ((self.data[key]['LDH_mean'].values - res_df['[LDH]'].values )  / (self.data[key]['LDH_mean'].max() - self.data[key]['LDH_mean'].min()))**2

            self.model.resetToOrigin()
            # print(Hb_error.sum(), LDH_error.sum())
            # TODO: weighted sum?
            error_sum += Hb_error.sum() + LDH_error.sum()
        return error_sum
