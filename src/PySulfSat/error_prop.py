import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.optimize as optimize
from scipy.special import erf
import warnings as w


from PySulfSat.core_calcs import *



def add_noise_1var(var, error_var,
error_type="Abs", error_dist="normal", N_dup=1000, sample_i=0, df_values=None):

    """ Takes 1 row of pd series or a float, returns 1000 duplicates following specified error position
    """
    sample_i=0
    if isinstance(var, pd.Series):
        len_loop=2
    else:
        len_loop=1

    if len_loop==1:
        df_c=pd.DataFrame(data={'var': var}, index=[0])
    else:
        df_c=pd.DataFrame(data={'var': var})


    # Temperature error distribution
    if error_type=='Abs':
        error_var=error_var
    if error_type =='Perc':
        error_var=df_c['var'].iloc[sample_i]*error_var/100
    if error_dist=='normal':
        Noise_to_add_var = np.random.normal(0, error_var, N_dup)
    if error_dist=='uniform':
        Noise_to_add_var = np.random.uniform(- error_var, +
                                                      error_var, N_dup)

    var_with_noise=Noise_to_add_var+df_c['var'].iloc[sample_i]



    return var_with_noise


def add_noise_series(var, error_var,  error_type="Abs", error_dist="normal", N_dup=1000, Sample_ID=None, no_noise=False, df_values=None):



    len_loop=len(var)
    All_outputs=pd.DataFrame([])
    Std_dev_var=np.empty(len_loop)
    Mean_var=np.empty(len_loop)
    Med_var=np.empty(len_loop)



    for i in range(0, len_loop):

        # If user has entered a pandas series for error, takes right one for each loop
        if type(error_var) is pd.Series:
            error_var=error_var.iloc[i]
        else:
            error_var=error_var

        if type(var) is pd.Series:
            var_i=var.iloc[i]
        else:
            var_i=var


        Error_MC=add_noise_1var(var=var_i,
        sample_i=i, error_var=error_var, N_dup=N_dup,
        error_dist=error_dist, error_type=error_type)






        # MC for each FI
        #All_outputs=pd.concat([All_outputs, Error_MC], axis=0)
        if i==0:
            All_outputs=Error_MC
        else:
            All_outputs=np.concatenate((All_outputs, Error_MC), axis=0)
    #     # get av and mean
    #     Std_dev_var[i]=np.nanstd(Error_MC)
    #     Mean_var[i]=np.nanmean(Error_MC)
    #     Med_var[i]=np.nanmedian(Error_MC)
    #
    #
    #
    # Av_outputs=pd.DataFrame(data={'Sample_ID': Sample,
    #                                   'Mean_var': Mean_var,
    #                                   'Std_var': Std_var,
    #                                    'Med_var': Med_var})



    return All_outputs

def duplicate_dataframe(df, N_dup=1000):
    df_values=df
    Dupdf = pd.DataFrame(np.repeat(df_values.values,
    N_dup, axis=0))
    Dupdf.columns = df_values.columns

    return Dupdf

def add_noise_2_dataframes(df_values, df_noise,
        error_type="Abs", error_dist="normal", N_dups=10):

    df1=df_values.copy()
    df2=df_noise.copy()

    #def create_noise_2df(df1, df2, error_type="Abs", error_dist="normal", N_dup=3):

    df2.columns = [col[:-4] if col.endswith('_Err') else col for col in df2.columns]
    df2 = df2[df1.columns]  # Reorder the columns to match df1

    # Now check if any columns were missing and return a warning if so
    if set(df1.columns) == set(df2.columns):
        print('columns match in 2 dataframes')
    else:
        # Find columns that are in df1 but not in df2
        print("Columns in df_values but not in df_noise:", set(df1.columns) - set(df2.columns))
        # Find columns that are in df2 but not in df1
        print("Columns in df_noise but not in df_values:", set(df2.columns) - set(df1.columns))
        raise TypeError('Columns in df_noise dont match those in df_values, cant do calculation Ammend inputs')

    # get only columns which are not objects
    num_cols = df1.select_dtypes(exclude=['object']).columns.tolist()

    #noisy_array=np.empty((N_dups*len(df1), len(num_cols)), dtype=float)
    col_i=0
    row_i=0
    for cols in num_cols:
        noise=add_noise_series(var=df1[cols], error_var=df2[cols], error_type=error_type, error_dist=error_dist, N_dup=N_dups)
        if col_i==0:
            noisy_array=noise
        else:
            noisy_array=np.vstack((noisy_array, noise))
        col_i=col_i+1
    noisy_array_T=noisy_array.T

    df_noisy=pd.DataFrame(noisy_array_T.tolist(), columns=num_cols)
    df_noisy

    # Get sample names
    if 'Sample_ID' in df1.columns:
        s=df1['Sample_ID']
    else:
        s=df1.index
    s_repeated = pd.Series(np.repeat(s.values, N_dups))

    # Tile the repeated values to get the end-on-end repetition
    s_end_on_end = pd.Series(np.tile(s_repeated, N_dups))


    df_noisy['Sample_ID']=s_end_on_end



    return df_noisy