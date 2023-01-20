import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.optimize as optimize
from scipy.special import erf



from PySulfSat.core_calcs import *
#SCAS (Li and Ripley, 2009;
#Baker and Moretti, 2011; Masotta and Keppler, 2015;
#Chowdhury and Dasgupta, 2019; Zajacz and Tsay, 2019)
df_ideal_liq_Chow = pd.DataFrame(columns=['SiO2_Liq', 'TiO2_Liq', 'Al2O3_Liq',
'FeOt_Liq', 'MnO_Liq', 'MgO_Liq', 'CaO_Liq', 'Na2O_Liq', 'K2O_Liq', 'H2O_Liq'])

def norm_liqs_with_H2O(Liqs):
    Liqs_c=Liqs.copy()
    Liqs_c=Liqs_c.fillna(0)
    Liqs2=Liqs_c.reindex(
        df_ideal_liq_Chow.columns, axis=1).fillna(0)

    sum_rows=Liqs2.sum(axis=1)

    Liqs_norm=Liqs2.divide(sum_rows/100, axis='rows')
    return Liqs_norm


def calculate_CD2019_SCAS(*, df, T_K, H2O_Liq=None):
    """" Calculates SCAS using the model of Chowdhury and Dasgupta, 2019


    Parameters
    -------
    df: pandas.DataFrame
        Dataframe of liquid compositions. Needs to have the headings "SiO2_Liq", "TiO2_Liq" etc, with
        compositions in oxide wt%. all FeO should be entered as FeOt_Liq

    T_K: int, float, pandas.Series
        Temperature in Kelvin.

    H2O_Liq: int, float, pandas.Series
        Option input, overwrites H2O_Liq in input dataframe




    Returns
    -------
    df of input dataframe, with new column heading "Calc SCAS (ppm)"

    """


    df_c=df.copy()
    df_norm=norm_liqs_with_H2O(Liqs=df_c)
    hyd_prop=calculate_hydrous_mol_proportions_liquid(liq_comps=df_norm)
    sum_hyd_prop=hyd_prop.sum(axis=1)
    hyd_frac=calculate_hydrous_mol_fractions_liquid(liq_comps=df_norm)
    sum_hydr_frac=hyd_frac.sum(axis=1)

    # Coefficients
    a=-13.23
    b=-0.5
    dsi=3.02
    dca=36.7
    dmg=2.84
    dfe=10.14
    dal=44.28
    dna=26.27
    dk=-25.77
    e=0.09
    f=0.54

    lnXS=(a+(b*(10**4/T_K))+dsi*hyd_frac['SiO2_Liq_mol_frac']
        +dca*hyd_frac['CaO_Liq_mol_frac']+dmg*hyd_frac['MgO_Liq_mol_frac']
        +dfe*hyd_frac['FeOt_Liq_mol_frac']+dal*hyd_frac['Al2O3_Liq_mol_frac']
        +dna*hyd_frac['Na2O_Liq_mol_frac']+dk*hyd_frac['K2O_Liq_mol_frac']
        +e*df_c['H2O_Liq']-f*np.log(hyd_frac['CaO_Liq_mol_frac']))
    Xs=np.exp(lnXS)
    molesS=Xs*(sum_hyd_prop+Xs)
    SCAS=molesS*32.065
    SCASppm=SCAS*10000
    # print(SCASppm)
    # print(len(SCASppm))
    out=pd.concat([df_c, hyd_prop, hyd_frac], axis=1)

    out.insert(0, 'Calc SCAS (ppm)', SCASppm)
    out.insert(1, 'lnXS', lnXS)
    out.insert(2, 'Xs', Xs)
    out.insert(3, 'molesS', molesS)
    return out

def calculate_ZT2022_SCAS(*, df, T_K, H2O_Liq=None):
    """" Calculates SCAS using the model of Zajacz and Tsay, 2022

    Parameters
    -------
    df: pandas.DataFrame
        Dataframe of liquid compositions. Needs to have the headings "SiO2_Liq", "TiO2_Liq" etc, with
        compositions in oxide wt%. all FeO should be entered as FeOt_Liq

    T_K: int, float, pandas.Series
        Temperature in Kelvin.

    H2O_Liq: int, float, pandas.Series
        Option input, overwrites H2O_Liq in input dataframe




    Returns
    -------
    df of input dataframe, with new column heading "Calc SCAS (ppm)"

    """


    df_c=df.copy()
    df_c['MnO_Liq']=0
    if H2O_Liq is not None:
        df_c['H2O_Liq']=H2O_Liq

    hyd_prop=calculate_hydrous_mol_proportions_liquid(liq_comps=df_c)
    sum_hyd_prop=hyd_prop.sum(axis=1)
    hyd_frac=calculate_hydrous_mol_fractions_liquid(liq_comps=df_c)
    sum_hydr_frac=hyd_frac.sum(axis=1)


    #NBO
    param1=(hyd_frac['Na2O_Liq_mol_frac']*2+hyd_frac['K2O_Liq_mol_frac']*2+2*(hyd_frac['CaO_Liq_mol_frac']
    +hyd_frac['MgO_Liq_mol_frac']+hyd_frac['FeOt_Liq_mol_frac'])-hyd_frac['Al2O3_Liq_mol_frac']*2)/(hyd_frac['Al2O3_Liq_mol_frac']*2+hyd_frac['SiO2_Liq_mol_frac'])
    param1_log=param1>0
    # If param1>0
    NBOT=(hyd_frac['Na2O_Liq_mol_frac']*2+hyd_frac['K2O_Liq_mol_frac']*2+2*(hyd_frac['CaO_Liq_mol_frac']+hyd_frac['MgO_Liq_mol_frac']+hyd_frac['FeOt_Liq_mol_frac'])-hyd_frac['Al2O3_Liq_mol_frac']*2)/(hyd_frac['Al2O3_Liq_mol_frac']*2+hyd_frac['SiO2_Liq_mol_frac'])
    NBOT[~param1_log]=0

    # P_Rhyolite
    param2=hyd_frac['K2O_Liq_mol_frac']+hyd_frac['Na2O_Liq_mol_frac']+hyd_frac['CaO_Liq_mol_frac']
    param2_log=param2>hyd_frac['Al2O3_Liq_mol_frac']
    P_Rhyo=3.11*(hyd_frac['K2O_Liq_mol_frac']+hyd_frac['Na2O_Liq_mol_frac']+hyd_frac['CaO_Liq_mol_frac']-hyd_frac['Al2O3_Liq_mol_frac'])
    P_Rhyo[~param2_log]=1.54*(hyd_frac['Al2O3_Liq_mol_frac']-(hyd_frac['CaO_Liq_mol_frac']+hyd_frac['Na2O_Liq_mol_frac']+hyd_frac['K2O_Liq_mol_frac']))

    # P_c
    P_c=(P_Rhyo+251*hyd_frac['CaO_Liq_mol_frac']**2+57*hyd_frac['MgO_Liq_mol_frac']**2+154*hyd_frac['FeOt_Liq_mol_frac']**2)/(2*hyd_frac['Al2O3_Liq_mol_frac']+hyd_frac['SiO2_Liq_mol_frac'])/(1+4.8*NBOT)
    P_T=np.exp(-7890/(T_K))
    P_H2O=hyd_frac['H2O_Liq_mol_frac']*(2.09-1.65*NBOT)+0.42*NBOT+0.23

    # Ksp
    Ksp=np.exp(1.226*np.log(P_c*P_T*P_H2O)+0.079)
    XSmelt=Ksp/hyd_frac['CaO_Liq_mol_frac']
    Smelt=XSmelt*sum_hyd_prop*32.07*10000

    df_c.insert(0, 'Calc SCAS (ppm)', Smelt)

    return df_c


## Masotta et al. (2013) calculations
def calculate_MK2015_SCAS(liq_comps, T_K):
    """ Calculates S6+ dissolved in the melt using Masotta and Kepler (2015)
    doi: https://doi.org/10.1016/j.gca.2015.02.033

    Parameters
    -------

    liq_comps: pandas.DataFrame
                Panda DataFrame of liquid compositions with column headings SiO2_Liq, TiO2_Liq etc.

    T_K: int, float, pandas.Series
        Temperature in Kelvin


    Returns
    -------
    pandas dataframe with calculated S in wt% and ppm, other intermediate
    calculations, and input dataframe
    """

    df=liq_comps.copy()
    mol_prop=calculate_anhydrous_mol_proportions_liquid(liq_comps=df)
    mol_frac=calculate_anhydrous_mol_fractions_liquid(liq_comps=df)
    cat_frac=calculate_anhydrous_cat_fractions_liquid(liq_comps=df)

    mol_prop_sum = mol_prop.sum(axis='columns')

    NBOT=((2*(2*cat_frac['Si_Liq_cat_frac']
    +2*cat_frac['Ti_Liq_cat_frac']+1.5*cat_frac['Al_Liq_cat_frac']
        +cat_frac['Fet_Liq_cat_frac']+cat_frac['Mg_Liq_cat_frac']
    +cat_frac['Ca_Liq_cat_frac']+0.5*cat_frac['Na_Liq_cat_frac']
    +0.5*cat_frac['K_Liq_cat_frac'])-4*(cat_frac['Si_Liq_cat_frac']
    +cat_frac['Ti_Liq_cat_frac']+cat_frac['Al_Liq_cat_frac']))
    /(cat_frac['Si_Liq_cat_frac']+cat_frac['Ti_Liq_cat_frac']
    +cat_frac['Al_Liq_cat_frac']))

    LnKps=(8.95+-146.5*NBOT+(-26960/(T_K))+197200*NBOT*(1/(T_K))
           +0.409*liq_comps['H2O_Liq'])

    # Wont match exactly as they include SO3 in their sum, which we dont know
    CaO_MF=((liq_comps['CaO_Liq']/56.0774)/(mol_prop_sum-mol_prop['CaO_Liq_mol_prop']
            +liq_comps['CaO_Liq']/56.0774))
    SO3_MF=np.exp(LnKps)/CaO_MF
    C_SO3_melt=SO3_MF*(mol_prop_sum-mol_prop['CaO_Liq_mol_prop']
                +liq_comps['CaO_Liq']/56.0774)*80.066
    S_melt_wt=C_SO3_melt*0.40048
    S_melt_ppm=S_melt_wt*10000

    new=pd.DataFrame(data={'S_melt_ppm': S_melt_ppm,
                              'S_melt_wt': S_melt_wt,
                              'C_SO3_melt': C_SO3_melt,
                              'SO3_MF': SO3_MF,
                              'CaO_MF': CaO_MF,
                              'LnKps': LnKps,
                              'NBOT': NBOT})
    df_out=pd.concat([new, df], axis=1)
    return df_out