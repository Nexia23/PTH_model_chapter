#packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from IPython.display import display


def rename_patientnumber(df):
    """ Convert patient number from the form "S14" to 14000"""
    df.patientnumber.str.replace('S', '').astype(int) * 1000
    return df

def rename_columns(df):
    """ deletes hastags in col-names, renames RPI-berechnet and V4_PTH"""
    df.columns = df.columns.str.replace('RPI-berechnet', 'RPI')
    df.columns =  df.columns.str.replace('#', '')
    df.columns  = df.columns.str.replace('V4_PTH', 'PTH')
    df.columns  = df.columns.str.replace('Parasite_%', 'parasitemia')
    df.columns =  df.columns.str.replace('OIE_per_µl', '[oiE]')
    df.columns =  df.columns.str.replace('Parasite_per_µl', '[iE]')
    df.columns =  df.columns.str.replace('Ery', '[E]')
    df.columns =  df.columns.str.replace('Reti-prozent', 'R_percent')
    df.columns =  df.columns.str.replace('Reti-abs', '[R]')
    df.columns =  df.columns.str.replace('OIE%', 'oiE_percent')
    return df

def create_immunity_column(df):
    df['Immunity'] = df['Groups'].str.extract('(semi-immune|non-immune)')
    return df

def switch_column_name_prefix_suffix(df):
    """ Returns Dataframe, where suffix and prefix of column_names are switched, e.g. 'V1-Hkt' ->'Hkt-V1'"""
    rename_columns(df)
    df.columns = ['-'.join(col.split('-')[::-1]) for col in df.columns]
    return df

def get_column_suffixes():
    """Creates list with Suffixe of col_names of dataframe, e.g. ['V1', 'V2', 'V3', 'V4', 'V5']"""
    suffixes = [f'V{i}' for i in range(1,6)]
    return suffixes

def get_column_stubnames(df):
    """ Creates list with Stubnames of col_names of dataframe, e.g.['Hkt', 'Ery', 'RPI_berechnet']"""
    stubnames = [col.split('-')[0] for col in df.columns if col.split('-')[-1] in get_column_suffixes()]
    stubnames = list(dict.fromkeys(stubnames))
    return stubnames

def get_index_columns(df) -> list:
    """returns list of col-names, who havent multiple timepoint measurements, e.g. patientnr."""
    ic = [col.split('-')[0] for col in df.columns if not col.split('-')[-1] in get_column_suffixes()] 
    return ic

def create_time_column(df):
    time_of_measurment = {'V1':0, 'V2':2, 'V3':6, 'V4':13,'V5':27}
    df['time'] = df['measurement'].map(time_of_measurment)
    return df

def drop_different_drug_rows(df):
    allowed_drug = 'Dihydroartemisinin-Piperaquine'
    df = df[(df['drug1'].isna() | (df['drug1'] ==  allowed_drug)) 
            & (df['drug2'].isna() | (df['drug2'] == allowed_drug))]
    return df

def drop_columns(df, columns_to_drop = [ 'OIE_any', 'measurement', 'drug1', 'drug2']): #,'malclass','age', 'sex', 'Groups', 'Symptombeginn']):
    df = df.drop(columns_to_drop, axis=1)
    return df

def add_patient_id(df):
    # add patient_id
    patient_names = df['patientnumber'].unique()

    for i, p_name in enumerate(patient_names):
        df.loc[df['patientnumber'] == p_name, 'patient_id'] = int(i)
    
    return df

def long_format(df):
    """ reshapes dataframe to long format"""
    #df_clean = rename_patientnumber(df)
    df_clean = rename_columns(df)
    df_clean = create_immunity_column(df_clean)
    df_clean = switch_column_name_prefix_suffix(df_clean)
    suffix = f"(!?{'|'.join(get_column_suffixes())})"
    stubnames = get_column_stubnames(df_clean)
    index = get_index_columns(df_clean)
    df = pd.wide_to_long(df_clean, stubnames= stubnames, i=index, sep='-', j='measurement', suffix=suffix)
    df = df.sort_values(['patientnumber', 'measurement'])
    df = df.reset_index()
    df.index.name = 'index'
    df = create_time_column(df)
    df = drop_different_drug_rows(df)
    df = drop_columns(df)
    df = add_patient_id(df)
    df['[E]'] *= 1e6  #Erys sind in Daten in x*1e6, deswegen Experiment mit 1e6 multiplizieren
    df['[R]'] *= 1e3  #Retis in daten per nl (mein modell in mikroliter)
    return df

def extract_fitting_data(df: pd.DataFrame, features: list=['Hb', 'LDH']):
    """ Returns two dataframes, with mean of features for PTH=1, one for PTH=0."""

    # extract means for features of pth patients
    pth_df = pd.DataFrame(columns=['Time'] + features)
    pth_df['Time'] = df['time'].unique()    
    for feature in features:
        pth_df[f'{feature}_mean'] = get_mean_of_feature(df, feature, pth=1)

    # extract means for features of non-pth patients
    nonpth_df = pd.DataFrame(columns=['Time'] + features)
    nonpth_df['Time'] = df['time'].unique()
    for feature in features:
        nonpth_df[f'{feature}_mean'] = get_mean_of_feature(df, feature, pth=0)

    return pth_df, nonpth_df  


def get_mean_of_feature(data: pd.DataFrame, feature: str, pth: int = 0):
    pth_data = data[data['PTH'] == pth]
    feature_array = pth_data[feature].values
    f_reshaped = feature_array.reshape(len(feature_array)//5 , 5)

    mean_feature = np.nanmean(f_reshaped, axis=0)
    # median_feature = np.nanmedian(f_reshaped, axis=0)  
    return mean_feature 

if __name__=='__main__':
    data_df = pd.read_excel('datasets/haemolysismodel_conRetis.xlsx')   #import Data von Pinkus 
    df =long_format(data_df)
    df.to_csv('datasets/OIE_data.csv', index=False)    
    pth, non_pth = extract_fitting_data(df, ['Hb', 'LDH','[R]']) 
    pth.to_csv('Estimation/pth.csv', index=False)
    non_pth.to_csv('Estimation/non_pth.csv', index=False)

