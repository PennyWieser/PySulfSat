import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



def propagate_s_in_sulfide_ind(sample_i=0,  N_dup=1000, S_Sulf=None, Vol=None,
sulf_dens=None, melt_dens=None,
error_S_Sulf=0, error_type_S_Sulf='Abs', error_dist_S_Sulf='normal',
error_Vol=0, error_type_Vol='Abs', error_dist_Vol='normal',
error_sulf_dens=0, error_type_sulf_dens='Abs', error_dist_sulf_dens='normal',
error_melt_dens=0, error_type_melt_dens='Abs', error_dist_melt_dens='normal', len_loop=1):

    """ This function propagates uncertainty in reconstruction of melt inclusion -sulfide volumes for a single row in a dataframe
    and returns a dataframe

    Parameters
    ----------------
    sample_i: int
        if your inputs are panda series, says which row to take

    N_dup: int
        Number of duplicates when generating errors for Monte Carlo simulations

    S_Sulf: int
        S in sulifde in wt%

    Vol: int, float, pd.series
        Volume percent of sulfide in melt inclusion

    sulf_dens: int, float, pd.series
        Density of the sulfide in kg/m3

    melt_dens:int, float, pd.series
        Density of the melt in kg/m3

    error_Vol, error_sulf_dens, error_melt_dens: int, float, pd.Series
        Error for each variable, can be absolute or %

    error_type_Vol, error_type_sulf_dens, error_type_melt_dens: 'Abs' or 'Perc'
        whether given error is perc or absolute

    error_dist_Vol, error_dist_sulf_dens, error_dist_melt_dens: 'normal' or 'uniform'
        Distribution of simulated error

    Returns
    ------------------
    pd.DataFrame:
        Input variable duplicated N_dup times with noise added.




    """


    if len_loop==1:
        df_c=pd.DataFrame(data={'S_Sulf': S_Sulf,
                            'Vol': Vol,
                            'sulf_dens': sulf_dens,
                            'melt_dens': melt_dens }, index=[0])
    else:
        df_c=pd.DataFrame(data={'S_Sulf': S_Sulf,
                            'Vol': Vol,
                            'sulf_dens': sulf_dens,
                            'melt_dens': melt_dens})


    # S in Sulf error distribution

    if error_type_S_Sulf=='Abs':
        error_S_Sulf=error_S_Sulf
    if error_type_S_Sulf =='Perc':
        error_S_Sulf=df_c['S_Sulf'].iloc[sample_i]*error_S_Sulf/100
    if error_dist_S_Sulf=='normal':
        Noise_to_add_S_Sulf = np.random.normal(0, error_S_Sulf, N_dup)
    if error_dist_S_Sulf=='uniform':
        Noise_to_add_S_Sulf = np.random.uniform(- error_S_Sulf, +
                                                      error_S_Sulf, N_dup)

    S_Sulf_with_noise=Noise_to_add_S_Sulf+df_c['S_Sulf'].iloc[sample_i]
    S_Sulf_with_noise[S_Sulf_with_noise < 0.000000000000001] = 0.000000000000001


    # Volume error distribution

    if error_type_Vol=='Abs':
        error_Vol=error_Vol
    if error_type_Vol =='Perc':
        error_Vol=df_c['Vol'].iloc[sample_i]*error_Vol/100
    if error_dist_Vol=='normal':
        Noise_to_add_Vol = np.random.normal(0, error_Vol, N_dup)
    if error_dist_Vol=='uniform':
        Noise_to_add_Vol = np.random.uniform(- error_Vol, +
                                                      error_Vol, N_dup)

    Vol_with_noise=Noise_to_add_Vol+df_c['Vol'].iloc[sample_i]
    Vol_with_noise[Vol_with_noise < 0.000000000000001] = 0.000000000000001

    #

    # Volume error distribution

    if error_type_sulf_dens=='Abs':
        error_sulf_dens=error_sulf_dens
    if error_type_sulf_dens =='Perc':
        error_sulf_dens=df_c['Vol'].iloc[sample_i]*error_sulf_dens/100
    if error_dist_sulf_dens=='normal':
        Noise_to_add_sulf_dens = np.random.normal(0, error_sulf_dens, N_dup)
    if error_dist_sulf_dens=='uniform':
        Noise_to_add_sulf_dens = np.random.uniform(- error_sulf_dens, +
                                                      error_sulf_dens, N_dup)

    sulf_dens_with_noise=Noise_to_add_sulf_dens+df_c['sulf_dens'].iloc[sample_i]
    sulf_dens_with_noise[sulf_dens_with_noise < 0.000000000000001] = 0.000000000000001

    # Volume error distribution

    if error_type_melt_dens=='Abs':
        error_melt_dens=error_melt_dens
    if error_type_melt_dens =='Perc':
        error_melt_dens=df_c['Vol'].iloc[sample_i]*error_melt_dens/100
    if error_dist_melt_dens=='normal':
        Noise_to_add_melt_dens = np.random.normal(0, error_melt_dens, N_dup)
    if error_dist_melt_dens=='uniform':
        Noise_to_add_melt_dens = np.random.uniform(- error_melt_dens, +
                                                      error_melt_dens, N_dup)

    melt_dens_with_noise=Noise_to_add_melt_dens+df_c['melt_dens'].iloc[sample_i]
    melt_dens_with_noise[melt_dens_with_noise < 0.000000000000001] = 0.000000000000001
    S_eq_melt_ind=10**4 * (0.01*df_c['Vol']*df_c['S_Sulf']*df_c['sulf_dens'])/df_c['melt_dens']
    df_out=pd.DataFrame(data={
        'S_Sulf_with_noise': S_Sulf_with_noise,
                                'Vol_with_noise': Vol_with_noise,
                                'sulf_dens_with_noise': sulf_dens_with_noise,
                                'melt_dens_with_noise': melt_dens_with_noise,
                                'S_Sulf': df_c['S_Sulf'].iloc[sample_i],
                                'Vol': df_c['Vol'].iloc[sample_i],
                                'Crustal Density_kg_m3': sulf_dens,
                                'error_S_Sulf': error_S_Sulf,
                                'error_type_S_Sulf': error_type_S_Sulf,
                                'error_dist_S_Sulf': error_dist_S_Sulf,
                                'error_Vol': error_Vol,
                                'error_type_Vol': error_type_Vol,
                                'error_dist_Vol': error_dist_Vol,
                                'error_sulf_dens': error_sulf_dens,
                                'error_type_sulf_dens': error_type_sulf_dens,
                                'error_dist_sulf_dens': error_dist_sulf_dens,
                                })

    S_eq_melt=10**4*((df_out['S_Sulf_with_noise']*0.01*df_out['Vol_with_noise']*df_out['sulf_dens_with_noise']))/(df_out['melt_dens_with_noise'])

    df_out['S_eq_melt']=S_eq_melt
    df_out['S_eq_melt_ind']=S_eq_melt_ind




    return df_out


## Now do it for all samples

def propagate_s_in_sulfide(sample_ID, Vol, S_Sulf, N_dup=1000,
sulf_dens=None, melt_dens=None,
error_S_Sulf=0, error_type_S_Sulf='Abs', error_dist_S_Sulf='normal',
error_Vol=0, error_type_Vol='Abs', error_dist_Vol='normal',
error_sulf_dens=0, error_type_sulf_dens='Abs', error_dist_sulf_dens='normal',
error_melt_dens=0, error_type_melt_dens='Abs', error_dist_melt_dens='normal'):



    """ This function propagates uncertainty in reconstruction of melt inclusion -sulfide volumes
    by feeding each row into propagate_s_in_sulfide_ind

    Parameters
    ----------------

    N_dup: int
        Number of duplicates when generating errors for Monte Carlo simulations

    S_Sulf: int
        S in sulifde in wt%

    Vol: int, float, pd.series
        Volume proportion of sulfide in melt inclusion

    sulf_dens: int, float, pd.series
        Density of the sulfide in kg/m3

    melt_dens:int, float, pd.series
        Density of the melt in kg/m3

    error_Vol, error_sulf_dens, error_melt_dens: int, float, pd.Series
        Error for each variable, can be absolute or %

    error_type_Vol, error_type_sulf_dens, error_type_melt_dens: 'Abs' or 'Perc'
        whether given error is perc or absolute

    error_dist_Vol, error_dist_sulf_dens, error_dist_melt_dens: 'normal' or 'uniform'
        Distribution of simulated error

    Returns
    ------------------
    pd.DataFrame: df_step, All_outputs
        All outputs has calculations for every simulaiton
        df_step has average and standard deviation for each sample

    """



    # Set up empty things to fill up.

    if type(Vol) is pd.Series:
        len_loop=len(Vol)
    else:
        len_loop=1



    mean_S_eq_melt = np.empty(len_loop, dtype=float)
    mean_S_eq_melt_ind = np.empty(len_loop, dtype=float)
    med_S_eq_melt  = np.empty(len_loop, dtype=float)
    std_S_eq_melt  = np.empty(len_loop, dtype=float)
    Sample=np.empty(len_loop,  dtype=np.dtype('U100') )


    All_outputs=pd.DataFrame([])


    #This loops through each sample
    for i in range(0, len_loop):
        if i % 20 == 0:
            print('working on sample number '+str(i))


        # If user has entered a pandas series for error, takes right one for each loop

        # S in the sulfide and error
        if type(error_S_Sulf) is pd.Series:
            error_S_Sulf_i=error_S_Sulf.iloc[i]
        else:
            error_S_Sulf_i=error_S_Sulf

        if type(S_Sulf) is pd.Series:
            S_Sulf_i=S_Sulf.iloc[i]
        else:
            S_Sulf_i=S_Sulf

        # Vol % and error


        if type(Vol) is pd.Series:
            Vol_i=Vol.iloc[i]
        else:
            Vol_i=Vol


        if type(error_Vol) is pd.Series:
            error_Vol_i=error_Vol.iloc[i]
        else:
            error_Vol_i=error_Vol

        # Sulf dens and error


        if type(sulf_dens) is pd.Series:
            sulf_dens_i=sulf_dens.iloc[i]
        else:
            sulf_dens_i=sulf_dens

        if type(error_sulf_dens) is pd.Series:
            error_sulf_dens_i=error_sulf_dens.iloc[i]
        else:
            error_sulf_dens_i=error_sulf_dens

        # Melt dens and error

        if type(melt_dens) is pd.Series:
            melt_dens_i=melt_dens.iloc[i]
        else:
            melt_dens_i=melt_dens


        if type(error_melt_dens) is pd.Series:
            error_melt_dens_i=error_melt_dens.iloc[i]
        else:
            error_melt_dens_i=error_melt_dens






        # This is the function doing the work to actually make the simulations for each variable.
        df_synthetic=propagate_s_in_sulfide_ind(N_dup=N_dup,
        S_Sulf=S_Sulf_i, Vol=Vol_i,
        sulf_dens=sulf_dens_i,
        melt_dens=melt_dens_i,
        error_S_Sulf=error_S_Sulf_i,
        error_type_S_Sulf=error_type_S_Sulf,
        error_dist_S_Sulf=error_dist_S_Sulf,
        error_Vol=error_Vol_i,
        error_type_Vol=error_type_Vol,
        error_dist_Vol=error_dist_Vol,
        error_sulf_dens=error_sulf_dens_i,
         error_type_sulf_dens=error_type_sulf_dens,
        error_melt_dens=error_melt_dens_i,
        error_type_melt_dens=error_type_melt_dens,
        error_dist_melt_dens=error_dist_melt_dens,
        error_dist_sulf_dens=error_dist_sulf_dens,
        len_loop=1)



        # Convert to densities for MC



        df=pd.DataFrame(data={'S_eq_melt': df_synthetic['S_eq_melt'],
                                'S_ind_eq_melt': df_synthetic['S_eq_melt_ind'] })






        # Check of
        if sample_ID is None:
            Sample[i]=i

        elif isinstance(sample_ID, str):
            Sample[i]=sample_ID
        else:
            Sample[i]=sample_ID.iloc[i]

        df.insert(0, 'Filename', Sample[i])







        All_outputs=pd.concat([All_outputs, df], axis=0)

        mean_S_eq_melt[i]=np.nanmean(df['S_eq_melt'])
        med_S_eq_melt[i]=np.nanmedian(df['S_eq_melt'])
        std_S_eq_melt[i]=np.nanstd(df['S_eq_melt'])

        mean_S_eq_melt_ind[i]=df['S_ind_eq_melt'].iloc[0]





    df_step=pd.DataFrame(data={'Filename': Sample,
                        'std_S_eq_melt': std_S_eq_melt,
                        'med_S_eq_melt': med_S_eq_melt,
                        'mean_S_eq_melt': mean_S_eq_melt,
                        'ind_s_eq_melt':mean_S_eq_melt_ind
                         })



    return df_step, All_outputs