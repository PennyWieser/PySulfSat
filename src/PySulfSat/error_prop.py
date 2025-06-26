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

    """ Takes 1 row of pd series, or a float, returns N_dup duplicates following specified error.


    Parameters
    ---------------
    var: pd.Series, float, int
        value to add noise to.

        if pd.Series specify sample_i, the row to use.

    error_var: int, float
        error to add (e.g. enter 5 to add 1sigma=5 Kelvin error to temperature)

    N_dup: int
        Number of duplicates to make (e.g., 1000 synthetic values)

    error_type: str, 'Abs' or 'Perc'
        Whether the inputted errors are absolute errors, or percentage errors


    err_dist: str, "normal" (default) or "uniform"
        whether the added error is normally or uniformly distributed

    Returns
    ----------
    panda series with a length of N_dup (e.g. var+noise, var+noise, var+noise... N_Dup)

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


def add_noise_series(var, error_var,  error_type="Abs", error_dist="normal", N_dup=1000):

    """  This function takes a panda series 'var' and adds noise from 'error_var', with attributions
    depending on the following parameters. It returns a panda series having duplicated each input row N_dup times, with noise added.
    (it relies on the function add_noise_1var to duplicate each individual row N times).

    Parameters
    ---------------
    var: pd.Series
        column of data you wish to add noise to.

    error_var:pd.Series, int, float
        error to add. Can be a column of errors with the same length as var, or a fixed error for all rows in var.

    N_dup: int
        Number of duplicates to make (e.g., 1000 synthetic values)

    error_type: str, 'Abs' or 'Perc'
        Whether the inputted errors are absolute errors, or percentage errors

    err_dist: str, "normal" (default) or "uniform"
        whether the added error is normally or uniformly distributed

    Returns
    ----------
    panda series with a length of N_dup (e.g. var+noise, var+noise, var+noise... N_Dup)

    """
    if isinstance(error_var, pd.Series):
        if len(var)!=len(error_var):
            raise TypeError('Input variable and error variable panda series are not the same length')



    len_loop=len(var)
    All_outputs=pd.DataFrame([])
    Std_dev_var=np.empty(len_loop)
    Mean_var=np.empty(len_loop)
    Med_var=np.empty(len_loop)



    for i in range(0, len_loop):

        # If user has entered a pandas series for error, takes right one for each loop
        if isinstance(error_var, pd.Series):
            this_error = error_var.iloc[i]
        else:
            this_error = error_var

        if isinstance(var, pd.Series):
            var_i = var.iloc[i]
        else:
            var_i = var


        Error_MC = add_noise_1var(
            var=var_i,
            sample_i=0,
            error_var=this_error,
            N_dup=N_dup,
            error_dist=error_dist,
            error_type=error_type
        )





        if i==0:
            All_outputs=Error_MC
        else:
            All_outputs=np.concatenate((All_outputs, Error_MC), axis=0)

    if isinstance(All_outputs, np.ndarray):
        All_outputs=pd.Series(All_outputs)


    return All_outputs

def av_noise_samples_series(calc, sampleID):
    '''
    This function calculates the mean, median, standard devation, maximum and
    minimum value of the values in the pd.Series "calc", averaging rows with the same value of SampleID.

    Parameters
    -------
    calc: Series
        Panda series of inputs you want to average.
    SampleID: str
        column heading for the thing you want to average by (e.g., Sample_ID_Cpx)

    Returns
    -------

    Dataframe with headings "Sample", "Mean_calc", "Median_calc",
    "St_dev_calc", "Max_calc", "Min_calc", with a number of rows corresponding to the
    number of unique values in SampleID

    '''


    if isinstance(calc, pd.Series):
        N = sampleID.unique()
        Av_mean = np.empty(len(N), dtype=float)
        Av_median = np.empty(len(N), dtype=float)
        Max = np.empty(len(N), dtype=float)
        Min = np.empty(len(N), dtype=float)
        Std = np.empty(len(N), dtype=float)
        i=0
        for ID in sampleID.unique():
            sam=ID
            # print(sam)
            # print(i)
            # print(np.nanmean(calc[sampleID == sam]))



            Av_mean[i] = np.nanmean(calc[sampleID == sam])
            Av_median[i] = np.nanmedian(calc[sampleID == sam])
            Std[i] = np.nanstd(calc[sampleID == sam])
            Min[i] = np.nanmin(calc[sampleID == sam])
            Max[i] = np.nanmax(calc[sampleID == sam])

            i=i+1
    len1=len(calc[sampleID == sam])
    Err_out = pd.DataFrame(data={'Sample': N, '# averaged': len1, 'Mean_calc': Av_mean,
    'Median_calc': Av_median, 'St_dev_calc': Std,
    'Max_calc': Max, 'Min_calc': Min})

    return Err_out





def duplicate_dataframe(df, N_dup=1000):
    """ This function takes a dataframe, and for each row, makes N_dup duplicates, end on end
    (e.g. row1-row2-row3 goes to row1-row1-row1-rowN_dup, row2-row2-row2-rowN_dup)

    Parameters
    --------------
    df: pd.DataFrame
        dataframe that needs duplicating. Can have whatever rows and columns you want

    N_dup: int
        Number of times to duplicate each row.


    Returns
    --------------
    pd.DataFrame
        dataframe with each row from original dataframe duplicated N times (Row1-row1-row1.... row2-row2-row2, row3-row3-row3)

    """

    df_values=df
    Dupdf = pd.DataFrame(np.repeat(df_values.values,
    N_dup, axis=0))
    Dupdf.columns = df_values.columns

    return Dupdf

def add_noise_2_dataframes(df_values, df_err,
        error_type="Abs", error_dist="normal", N_dups=10, sample_name_col='Sample_ID'):
    """ This function takes each value in df_values and adds noise in the corresponding row and column from df_err, and duplicates the row N times.
    e.g. for N_dups=10, 10 duplicate rows are made for each row in df_values, with the values pertubed based on noise in df_err

    Parameters
    --------------
    sample_name_col: str
        Name of column with sample name in each dataframe.

    df_values: pd.DataFrame
        dataframe of preferred values

    df_err: pd.DataFrame
        dataframe of error values. Needs to have the same number of rows, and exact same columns as df_values, but with the string '_Err' after each column.
        The only exception is the column 'sample_name_col' doesnt have to have _Err after. Note, the order has to be the same,
        the function doesnt automatically match rows based on sample names (it could, but wit would be slower than turning it all into numpy arrays and doing matrix math.

    N_dup: int
        Number of duplicates to make (e.g., 1000 synthetic values)

    error_type: str, 'Abs' or 'Perc'
        Whether the inputted errors are absolute errors, or percentage errors

    err_dist: str, "normal" (default) or "uniform"
        whether the added error is normally or uniformly distributed

    Returns
    --------------
    pd.DataFrame with the length of the original dataframe*  N_Dups

    """
    if len(df_values)!=len(df_err):
        raise InputError('Length of df1 and df2 need to be the same')
    df1=df_values.copy()
    df2=df_err.copy()

    # Dealing with sample names

    if sample_name_col in df1.columns:
        s=df1[sample_name_col]
        df1=df1.drop(sample_name_col, axis=1)

    else:
        s=df1.index

    if sample_name_col in df2.columns:

        df2=df2.drop(sample_name_col, axis=1)


    s_repeated = pd.Series(np.repeat(s.values, N_dups))


    #def create_noise_2df(df1, df2, error_type="Abs", error_dist="normal", N_dup=3):

    df2.columns = [col[:-4] if col.endswith('_Err') else col for col in df2.columns]
    df2 = df2[df1.columns]  # Reorder the columns to match df1

    # Now check if any columns were missing and return a warning if so
    if set(df1.columns) == set(df2.columns):
        print('Yay. columns match in 2 dataframes')
    else:
        # Find columns that are in df1 but not in df2
        print("Columns in df_values but not in df_err:", set(df1.columns) - set(df2.columns))
        # Find columns that are in df2 but not in df1
        print("Columns in df_err but not in df_values:", set(df2.columns) - set(df1.columns))
        raise TypeError('Columns in df_err dont match those in df_values, cant do calculation Ammend inputs')

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


    # Tile the repeated values to get the end-on-end repetition
    s_end_on_end = pd.Series(np.tile(s_repeated, N_dups))


    df_noisy[sample_name_col]=s_end_on_end



    return df_noisy