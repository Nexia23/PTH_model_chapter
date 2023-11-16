#packages
import pandas as pd
import numpy as np

class data(object):

    def __init__(self, csvpath ,imp_type='excel') -> None:
        imp = {
            'excel': pd.read_excel,
            'csv': pd.read_csv
        }
        self.df = imp[imp_type](csvpath)
        self.longdf = long_format(self.df)
        self.patient_count = self.df.size

    def __repr__(self):
        return self.df.to_string    

    def filter_df(self, filter=0):
        pass    

def cols_renaming(df):
    """
    deletes hashtags in col-names, renames RPI-berechnet and V4_PTH
    """
    df.columns = [col.replace('RPI-berechnet', 'RPI_berechnet') for col in df.columns]
    df.columns = [col.replace('#', '') for col in df.columns]
    df.columns = [col.replace('V4_PTH', 'PTH') for col in df.columns]
    #df.columns = [col.replace('Groups', 'Immunity') for col in df.columns]
    return df

def rearrange_immun(df):
    col = 'Groups'
    new_col = 'Immunity'
    immun_type = ['semi-immune', 'non-immune']

    for i in immun_type:
        mask_immun = df[col].str.contains(i, na = False, regex = False)
        df.loc[mask_immun, new_col] = i
    return df

def colnames_fixswitch(df):
    """
    Returns Dataframe, where suffix and prefix of column_names are switched, e.g. 'V1-Hkt' ->'Hkt-V1'
    """
    cols_renaming(df)
    df.columns = ['-'.join(col.split('-')[::-1]) for col in df.columns]
    return df

def colnames_suffixes():
    """
    Creates list with Suffixe of col_names of dataframe, e.g. ['V1', 'V2', 'V3', 'V4', 'V5']
    """
    suffixes = [f'V{i}' for i in range(1,6)]
    return suffixes

def cols_stubnames(df):
    """
    Creates list with Stubnames of col_names of dataframe, e.g.['Hkt', 'Ery', 'RPI_berechnet']
    """
    sn = [col.split('-')[0] for col in df.columns if col.split('-')[-1] in colnames_suffixes()]
    sn = list(dict.fromkeys(sn))
    return sn

def index_col(df) -> list:
    """
    returns list of col-names, who havent multiple timepoint measurements, e.g. patientnr.
    """
    ic = [col.split('-')[0] for col in df.columns if not col.split('-')[-1] in colnames_suffixes()] 
    return ic

def long_format(df):
    """ reshapes dataframe to long format"""
    df_clean = cols_renaming(df)
    df_clean = rearrange_immun(df_clean)
    df_clean = colnames_fixswitch(df_clean)
    suffix = f"(!?{'|'.join(colnames_suffixes())})"
    stubnames = cols_stubnames(df_clean)
    index = index_col(df_clean)

    wide_long_df = pd.wide_to_long(df_clean, stubnames= stubnames, i=index, sep='-', j='Timepoint', suffix=suffix)
    long_df = wide_long_df.reset_index()
    return long_df

if __name__=='__main__':
    D = data('/haemolysismodel_conRetis.xlsx')
    D.df.to_string
    df = D.longdf
    print(df)
    #D.patient_count


""""
data_df2 = pd.read_excel('/haemolysismodel_conRetis.xlsx')
test_columns = ['patientnumber', 'age'] + [col for col in data_df2.columns if col.endswith(('Hkt', 'Ery', 'RPI-berechnet'))]
test_df2 = data_df2[test_columns].head(10)
test_df2

#test_df2.columns

# Column renaming
new_cols = [col.replace('RPI-berechnet', 'RPI_berechnet').replace('#', '') for col in test_df2.columns]

# switch suffix and prefix
new_cols = ['-'.join(col.split('-')[::-1]) for col in new_cols]

new_cols
# make new df with updated columns
new_col_df = test_df2.copy()
new_col_df.columns = new_cols
new_col_df

suffixes = [f'V{i}' for i in range(1,6)]
stubnames = [col.split('-')[0] for col in new_col_df.columns if col.split('-')[-1] in suffixes]
index_col = [col.split('-')[0] for col in new_col_df.columns if not col.split('-')[-1] in suffixes]
stubnames = list(dict.fromkeys(stubnames))

suffixes_sorted = f"(!?{'|'.join(suffixes)})"

long_df = pd.wide_to_long(new_col_df, stubnames=stubnames, i=index_col, sep='-', j='Time', suffix=suffixes_sorted)

final = long_df.reset_index()
"""