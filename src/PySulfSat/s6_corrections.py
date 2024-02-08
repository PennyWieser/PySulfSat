import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import scipy
import scipy.optimize as optimize
from scipy.special import erf

# import warnings
# warnings.simplefilter('error')



from PySulfSat.core_calcs import *

def calculate_SCSS_Total(SCSS, S6St_Liq):
    '''
    Calculates SCSS Total from the SCSS (2-)
    following Jugo et al. (2009, 2010)

    Parameters
    -----------
    SCSS: int, float, pandas.Series
        SCSS from different calculators (only uses S2-)

    S6St_Liq: int, float, pandas.Series
        Proportion of S6+ in the liquid

    Returns
    -----------
    float, pandas.Series
        SCSS total, i.e inputted SCSS (2-) scaled up to account for S present in the S6+ form.

    '''
    SCSS_Total=SCSS/(1-S6St_Liq)
    return SCSS_Total


def calculate_SCAS_Total(SCAS, S2St_Liq=None, S6St_Liq=None):
    '''
    Calculates SCAS Total from the SCAS (6+)
    following Jugo et al. (2009, 2010)

    Parameters
    -----------
    SCAS: int, float, pandas.Series, np.array
        SCAS from different calculators (only uses S6+)

    S2St_Liq: int, float, pandas.Series, np.array
        Proportion of S2- in the liquid

    or
    S6St_Liq: proportion of S6+ in the liquid

    Returns
    -----------
    float, pandas.Series, np.array
        SCAStotal, i.e inputted SCAS (6+) scaled up to account for S present in the S2- form.

    '''
    if S6St_Liq is not None:
        S2St_Liq=1-S6St_Liq
    SCAS_Total=SCAS/(1-S2St_Liq)
    return SCAS_Total

def calculate_S6St_Jugo2010_eq10(deltaQFM):
    '''
    Calculates The S6/St ratio from deltaQFM using
    Jugo et al. (2010) eq10

    Parameters
    -----------
    deltaQFM: int, float, pandas.Series, np.array
        Offset from the QFM buffer in log units. Make sure relative to frost buffer

    Returns
    -----------
    float, pandas.Series, np.array
        Proportion of S6+ in the liquid

    '''
    S6St_Liq=1/(1+10**(2.1-2*deltaQFM))

    return S6St_Liq

def calculate_S2St_Jugo2009_eq8(deltaQFM):
    '''
    Calculates The S2/St ratio from deltaQFM using
    Jugo (2009) eq8 (re-arranged).

    Parameters
    -----------
    deltaQFM: int, float, pandas.Series, np.array
        Offset from the QFM buffer in log units.

    Returns
    -----------
    float, pandas.Series, np.array
        Proportion of S6+ in the liquid

    '''
    S2St_Liq=1/(1+np.exp(2.23*deltaQFM-2.89))

    return S2St_Liq


def calculate_S6St_Nash2019(T_K, Fe3Fet_Liq):
    '''
    Calculates The S6/St ratio from the temperature and the proportion
    of Fe3+ in the liquid using Nash et al. (2019).

    Parameters
    -----------
    T_K: int, float, pandas.Series, np.array
        Temperature in Kelvin.

    Fe3Fet_Liq: int, float, pandas.Series, np.array
        Proportion of Fe3 in the liquid
    Returns
    -----------
    float, pandas.Series,  np.array
        Proportion of S6+ in the liquid

    '''
    # First, calculate Fe3/Fe2 from Fe3/Fet
    Fe3Fe2=Fe3Fet_Liq/(1-Fe3Fet_Liq)
    #print(Fe3Fe2)
    if isinstance(Fe3Fe2, float) or isinstance(Fe3Fe2, int):
        log_Fe3Fe2=np.log10(Fe3Fe2)

    else:
        log_Fe3Fe2=np.log10(Fe3Fe2.astype(float))


    log_S6S2=8*log_Fe3Fe2 + 8.7436*10**6/(T_K**2) - 27703/T_K + 20.273
    S6S2=10**(log_S6S2)

    S6=S6S2/(1+S6S2)

    S2=1/(1+S6S2)
    S6St=S6/(S6+S2)

    return S6St


def calculate_S_Tot_Kleinsasser2022_dacite(*, SCSS2=None,  SCAS=None, deltaQFM=None):
    """ Calculates SCSS total as a function of S2- and S6+ in dacitic melts,
    following Kleinsasser et al. 2022

    Parameters
    -----------
    SCSS2: int, float, pd.Series, np.array
        SCSS in ppm

    SCAS: int, float, pd.Series, np.array
        SCAS in ppm

    deltaQFM: int, float, pd.Series, np.array

    Returns
    --------
    same format as inputs
        Proportion of S6+ in the liquid

    """

    if SCSS2 is not None:
        SCSS_Tot=SCSS2*(1+10**(2*deltaQFM-3.05))
        return SCSS_Tot
    if SCAS is not None:
        SCAS_Tot=SCAS*(1+np.exp(1.26-2*deltaQFM))
        return SCAS_Tot



def calculate_S_Tot_Jugo2010(*, SCSS2=None,  SCAS=None,deltaQFM):
    '''Calculates S total as a function of S2- and S6+
    following Jugo et al. 2010 for given SCSS or SCAS

    Parameters
    -----------
    deltaQFM: int, float, pandas.Series
        Offset from the QFM buffer in log units.

    Returns
    -----------
    float, pandas.Series
        Proportion of S6+ in the liquid

    '''
    if SCSS2 is not None:
        SCSS_Tot=SCSS2*(1+10**(2*deltaFMQ-2.1))
        return SCSS_Tot

        SCAS_Tot=SCAS*(1+np.exp(2.89-2.23*deltaFMQ))
        return SCAS_Tot


def calculate_S_Total_SCSS_SCAS(*, SCSS, SCAS, deltaQFM=None,  model=None, S6St_Liq=None,
                                T_K=None, Fe3Fet_Liq=None, df=None, logfo2=None):
    """
    Calculates the total amount of S accounting for the SCSS and SCAS

    Parameters
    -----------
    SCSS: int, float, pandas.Series
        SCSS (2-) you have calculated from whatever model you want.

    SCAS: int, float, pandas.Series
        SCAS (6+) you have calculated from whatever model you want.

    Choose either
    S6St_Liq: int, float, str
        S6/ST in the liquid

    Or if you dont know this, specify a model to calculate it.

    model: str
        Model used to calculate S6St_Liq, choice of:
        'Nash': Uses Nash et al. (2019) based on Fe3Fet_Liq and T_K
        'Jugo': Uses Jugo et al. (2010) Eq 10 based on DeltaQFM
        'Kleinsasser': Uses Kleinsasser et al. 2022 based on deltaQFM
        'OM2022': Need to also enter your dataframe of liquid compositions.

    deltaQFM: int, float, panda.Series (for Jugo et al. 2010, or Kleinsasser et al. 2022)

    logfo2: int, float, panda.Series (works for OM2022)

    Fe3Fet_Liq: int, float, panda.Series
        Proportion of Fe3Fet in the liquid, needed if model='Nash'

    T_K: int, float, panda.Series
        Temperature in Kelvin, needed if model='Nash'

    df: panda.DataFrame
        Option, needed for OM2022




    """
    if S6St_Liq is not None and model is not None:
        raise TypeError('Please specify either S6St_Liq or a model to calculate this parameter, not both')

    if model =="Kleinsasser":
        SCSS_Tot=calculate_S_Tot_Kleinsasser2022_dacite(SCSS2=SCSS,
                                                           deltaQFM=deltaQFM)
        SCAS_Tot=calculate_S_Tot_Kleinsasser2022_dacite(SCAS=SCAS,
                                                           deltaQFM=deltaQFM)
        S6St_Liq=np.nan




    else:

        if model=='OM2022':
            if T_K is None:
                raise TypeError('You need to specify a temp')
            if df is None:
                raise TypeError('You need to input a liquid dataframe')
            if Fe3Fet_Liq is None and logfo2 is None:
                raise TypeError('you need to input either Fe3Fet_Liq or logfo2')
            if Fe3Fet_Liq is not None:
                calcS_OM2022_GivenFe3=calculate_OM2022_S6St(df=df, T_K=T_K,
                            Fe3Fet_Liq=Fe3Fet_Liq)
            if logfo2 is not None:
                calcS_OM2022_GivenFe3=calculate_OM2022_S6St(df=df,
                        T_K=T_K,
                        logfo2=logfo2)
            S6St_Liq=calcS_OM2022_GivenFe3['S6St_Liq']

        if model =="Jugo":
            S6St_Liq=calculate_S6St_Jugo2010_eq10(deltaQFM=deltaQFM)

        if model =="Nash":
            if T_K is None:
                raise Exception('Need a T_K input to use Nash')
            if Fe3Fet_Liq is None:
                raise Exception('Need a Fe3Fet_Liq input to use Nash')
            S6St_Liq=calculate_S6St_Nash2019(T_K=T_K, Fe3Fet_Liq=Fe3Fet_Liq)

        if S6St_Liq is not None:
            S6St_Liq=S6St_Liq


        SCSS_Tot=calculate_SCSS_Total(SCSS=SCSS, S6St_Liq=S6St_Liq)
        SCAS_Tot=calculate_SCAS_Total(SCAS=SCAS, S2St_Liq=1-S6St_Liq)



    SCSS_S6_cont=SCSS_Tot-SCSS
    SCAS_S2_cont=SCAS_Tot-SCAS

    if isinstance(SCSS_Tot, float):
        df_Species=pd.DataFrame(data={'deltaQFM': deltaQFM,
                                    'S6St_Liq': S6St_Liq,
                                    'SCSS_2_ppm': SCSS,
                                'SCAS_6_ppm': SCAS,
                                'SCSS_Tot': SCSS_Tot,
                                'SCAS_Tot': SCAS_Tot,
                                    'S6 in SCSS_Tot': SCSS_S6_cont,
                                'S2 in SCAS_Tot': SCAS_S2_cont
                                }, index=[0])
    else:

        df_Species=pd.DataFrame(data={'deltaQFM': deltaQFM,
                                    'S6St_Liq': S6St_Liq,
                                    'SCSS_2_ppm': SCSS,
                                'SCAS_6_ppm': SCAS,
                                'SCSS_Tot': SCSS_Tot,
                                'SCAS_Tot': SCAS_Tot,
                                    'S6 in SCSS_Tot': SCSS_S6_cont,
                                'S2 in SCAS_Tot': SCAS_S2_cont
                                })

    # Cant have more S6 than the SCAS


    toomuchS6=df_Species['S6 in SCSS_Tot']>df_Species['SCAS_6_ppm']
    df_Species['SCSS_Tot_check']=df_Species['SCSS_Tot']
    df_Species.loc[toomuchS6, 'SCSS_Tot_check']=df_Species['SCAS_6_ppm']+df_Species['SCSS_2_ppm']

    # Cant have more S2- than the SCSS
    toomuchS2=df_Species['S2 in SCAS_Tot']>df_Species['SCSS_2_ppm']
    df_Species['SCAS_Tot_check']=df_Species['SCAS_Tot']
    # Ignore this error for now, seems a panda issue. wasted a lot of time on it.

    df_Species['SCAS_6_ppm'] = pd.to_numeric(df_Species['SCAS_6_ppm'], errors='coerce')
    df_Species['SCSS_2_ppm'] = pd.to_numeric(df_Species['SCSS_2_ppm'], errors='coerce')
    df_Species.loc[toomuchS2, 'SCAS_Tot_check']=df_Species['SCAS_6_ppm']+df_Species['SCSS_2_ppm']

    #Set as the minimum one
    df_Species.insert(0, 'Total_S_ppm',df_Species['SCAS_Tot_check'])
    SCAS_higher=df_Species['SCAS_Tot_check']>df_Species['SCSS_Tot_check']
    df_Species.loc[SCAS_higher, 'Total_S_ppm']=df_Species['SCSS_Tot_check']

    df_Species.insert(1, 'S2_Tot_ppm', df_Species['Total_S_ppm']*(1-df_Species['S6St_Liq']))
    df_Species.insert(2, 'S6_Tot_ppm', df_Species['Total_S_ppm']*(df_Species['S6St_Liq']))

    df_Species2=df_Species.copy()
    df_Species2.drop(['SCAS_Tot_check', 'SCSS_Tot_check', 'SCSS_Tot_check'], axis=1, inplace=True)
    return df_Species2


def calculate_BW2022_CS6(*, df, T_K):
    """ Calculates logCs6 and Cs6 using the expression of
    Boulliung and Wood, 2022. Also converts into the same form as ONeill and Mavrogenes

    Parameters
    -----------
    df: pandas Dataframe
        input dataframe of liquid compositions with _Liq after each oxide


    T_K: int, float, pd.series
        Temperature in Kelvin

    Returns
    --------------
    pd.DataFrame:
        includes columns for CS6 in BW22 and O22 format, and cat fractions


    """

    liqs=df.copy()
    liqs=calculate_anhydrous_cat_fractions_liquid(liq_comps=liqs)

    logCs6=(-12.659+(3692*liqs['Ca_Liq_cat_frac'] - 7592*liqs['Si_Liq_cat_frac']
        -13736*liqs['Ti_Liq_cat_frac']+3762*liqs['Al_Liq_cat_frac']+34483)/T_K)
    lnK=55921/T_K-25.07+0.6465*np.log(T_K)

    lnCs6_ONeill_Format=logCs6*np.log(10) + np.log10(10**4)*np.log(10)-lnK
    liqs.insert(0, 'LogCS6_calc_BW22_format', logCs6)
    liqs.insert(1, 'LnCS6_calc_OM22_format', lnCs6_ONeill_Format)


    return liqs



def calculate_OM2022_CS6(df, T_K, Fe3Fet_Liq=None, logfo2=None):
    """ Calculates ln Cs6+ using ONeill and Mavrogenes (2022).
    Also converts to logCs6 in the form of Boulling and Wood (2022)

    Parameters
    -----------
    df: pandas Dataframe
        input dataframe of liquid compositions with _Liq after each oxide


    T_K: int, float, pd.series
        Temperature in Kelvin

    Fe3Fet_Liq: float, pd.series
        Fe3/Fet ratio in the liquid

    logfo2: float, pd.series
        logfo2 value


    Returns
    -----------
    pd.DataFrame
        input dataframe, with Cs6 calculated both Oneill and Bouilling way.


    """
    liqs=df.copy()



    liqs=calculate_anhydrous_cat_fractions_liquid(liq_comps=liqs)


    if Fe3Fet_Liq is not None:
        C5=1
    if logfo2 is not None:
        C5=2
        logfo2_calc=logfo2

    if Fe3Fet_Liq is not None and logfo2 is not None:
        raise Exception('enter one of Fe3FetLIq or logfo2, not both as this is ambigous')

    if C5==1: # specify Fe3Fet not logfo2
        deltaQFM=(4*(np.log10(Fe3Fet_Liq/(1-Fe3Fet_Liq))+1.36-2*liqs['Na_Liq_cat_frac']
                    -3.7*liqs['K_Liq_cat_frac']-2.4*liqs['Ca_Liq_cat_frac']))
        logfo2_calc=deltaQFM-25050/T_K+8.58
        liqs['deltaQFM_calc']=deltaQFM
        liqs['logfo2_calc']=logfo2_calc
    if C5==2: #If specify Fe3Fet
        deltaQFM=logfo2-25050/T_K+8.58
        liqs['deltaQFM_calc']=deltaQFM

    if C5==2: # e.g. if Fe3Fet not given
        Fe2Fet_Liq=1/(1+10**(0.25*deltaQFM-1.36+2*liqs['Na_Liq_cat_frac']+2.4*liqs['Ca_Liq_cat_frac']+3.7*liqs['K_Liq_cat_frac']))
        liqs['Fe2Fet_Liq_calc']=Fe2Fet_Liq
    if C5==1: # IF Fe3Fet is given
        Fe2Fet_Liq=1-Fe3Fet_Liq

    liqs['Fe2_Liq_cat_frac']=liqs['Fet_Liq_cat_frac']*Fe2Fet_Liq

    LnCS6_calc=(-8.02+(21100+44000*liqs['Na_Liq_cat_frac']+18700*liqs['Mg_Liq_cat_frac']
    +4300*liqs['Al_Liq_cat_frac']+35600*liqs['Ca_Liq_cat_frac']
    +44200*liqs['K_Liq_cat_frac']+16500*liqs['Fe2_Liq_cat_frac']+12600*liqs['Mn_Liq_cat_frac'])/T_K)
    lnK=55921/T_K-25.07+0.6465*np.log(T_K)
    logCS6_BW22=(LnCS6_calc+lnK)/np.log(10) - np.log10(10**4)

    liqs.insert(0, 'LnCS6_calc_OM22_format', LnCS6_calc)
    liqs.insert(1, 'LogCS6_calc_BW22_format', logCS6_BW22 )

    return liqs

def calculate_OM2022_S6St(df, T_K, logfo2=None,
                    Fe3Fet_Liq=None):
    """
    Calculates S6/ST (as well as ln Cs2- and ln Cs6+) Using ONeill and Mavrogenes (2022)

    Parameters
    -----------
    df: pandas Dataframe
        input dataframe of liquid compositions with _Liq after each oxide


    T_K: int, float, pd.series
        Temperature in Kelvin


    Fe3Fet_Liq: int, float, pd.Series
        Fe3Fet ratio in the liquid

    OR

    logfo2: int, float, pd.Series
        logfo2 value

    Returns
    -----------

    pd.DataFrame: Contains S6/ST, Cs2- and Cs6+, and all intermediate calculations + user inputs


    """

    liqs=df.copy()
    liqs=calculate_anhydrous_cat_fractions_liquid(liq_comps=liqs)
    T_K=T_K-0.15
    if Fe3Fet_Liq is not None:
        C5=1
    if logfo2 is not None:
        C5=2
        logfo2_calc=logfo2

    if Fe3Fet_Liq is not None and logfo2 is not None:
        raise Exception('enter one of Fe3FetLIq or logfo2, not both as this is ambigous')

    if C5==1: # specify Fe3Fet not logfo2
        deltaQFM=(4*(np.log10(Fe3Fet_Liq/(1-Fe3Fet_Liq))+1.36-2*liqs['Na_Liq_cat_frac']
                    -3.7*liqs['K_Liq_cat_frac']-2.4*liqs['Ca_Liq_cat_frac']))
        logfo2_calc=deltaQFM-25050/T_K+8.58
        liqs['deltaQFM_calc']=deltaQFM
        liqs['logfo2_calc']=logfo2_calc
    if C5==2: #If specify logfo2
        deltaQFM=logfo2-8.58+25050/T_K
        liqs['deltaQFM_calc']=deltaQFM

    if C5==2: # e.g. if Fe3Fet not given
        Fe2Fet_Liq=1/(1+10**(0.25*deltaQFM-1.36+2*liqs['Na_Liq_cat_frac']+2.4*liqs['Ca_Liq_cat_frac']+3.7*liqs['K_Liq_cat_frac']))
        liqs['Fe2Fet_Liq_calc']=Fe2Fet_Liq
    if C5==1: # IF Fe3Fet is given
        Fe2Fet_Liq=1-Fe3Fet_Liq

    liqs['Fe2_Liq_cat_frac']=liqs['Fet_Liq_cat_frac']*Fe2Fet_Liq

    liqs['LnCS2_calc']=(8.77-23590/T_K+(1673/T_K)*(6.7*(liqs['Na_Liq_cat_frac']+liqs['K_Liq_cat_frac'])
        +4.9*liqs['Mg_Liq_cat_frac']+8.1*liqs['Ca_Liq_cat_frac']+8.9*(liqs['Fet_Liq_cat_frac']+liqs['Mn_Liq_cat_frac'])
        +5*liqs['Ti_Liq_cat_frac']+1.8*liqs['Al_Liq_cat_frac']
        -22.2*liqs['Ti_Liq_cat_frac']*(liqs['Fet_Liq_cat_frac']+liqs['Mn_Liq_cat_frac'])
            +7.2*((liqs['Fet_Liq_cat_frac']+liqs['Mn_Liq_cat_frac'])*liqs['Si_Liq_cat_frac']))-2.06*erf(-7.2*(liqs['Fet_Liq_cat_frac']+liqs['Mn_Liq_cat_frac'])))

    liqs['LnCS6_calc']=(-8.02+(21100+44000*liqs['Na_Liq_cat_frac']+18700*liqs['Mg_Liq_cat_frac']
    +4300*liqs['Al_Liq_cat_frac']+35600*liqs['Ca_Liq_cat_frac']
    +44200*liqs['K_Liq_cat_frac']+16500*liqs['Fe2_Liq_cat_frac']+12600*liqs['Mn_Liq_cat_frac'])/T_K)

    liqs['LnKSO2S2']=-55921/T_K+25.07-0.6465*np.log(T_K)

    liqs['LnS6S2']=(liqs['LnCS6_calc']-liqs['LnKSO2S2']-liqs['LnCS2_calc'])+2*np.log(10)*logfo2_calc
    liqs['S6St_Liq']=1-1/(1+np.exp(liqs['LnS6S2']))

    cols_to_move = ['S6St_Liq', 'LnCS2_calc', 'LnCS6_calc', 'LnKSO2S2', 'LnS6S2',
                    'deltaQFM_calc']
    liqs = liqs[cols_to_move +
                                    [col for col in liqs.columns if col not in cols_to_move]]
    return liqs


def calculate_BW2022_OM2022_S6St(df, T_K, logfo2=None,
                    Fe3Fet_Liq=None):
    """
    Calculates S6/ST (as well as ln Cs2- and ln Cs6+) Using ONeill and Mavrogenes (2022)
    for CS2- and Boulliung and Wood (2022) for CS6+

    Parameters
    -----------
    df: pandas Dataframe
        input dataframe of liquid compositions with _Liq after each oxide


    T_K: int, float, pd.series
        Temperature in Kelvin


    Fe3Fet_Liq: int, float, pd.Series
        Fe3Fet ratio in the liquid

    OR

    logfo2: int, float, pd.Series
        logfo2 value

    Returns
    -----------

    pd.DataFrame: Contains S6/ST, Cs2- and Cs6+, and all intermediate calculations + user inputs


    """

    liqs=df.copy()
    liqs=calculate_anhydrous_cat_fractions_liquid(liq_comps=liqs)
    T_K=T_K-0.15
    if Fe3Fet_Liq is not None:
        C5=1
    if logfo2 is not None:
        C5=2
        logfo2_calc=logfo2

    if Fe3Fet_Liq is not None and logfo2 is not None:
        raise Exception('enter one of Fe3FetLIq or logfo2, not both as this is ambigous')

    if C5==1: # specify Fe3Fet not logfo2
        deltaQFM=(4*(np.log10(Fe3Fet_Liq/(1-Fe3Fet_Liq))+1.36-2*liqs['Na_Liq_cat_frac']
                    -3.7*liqs['K_Liq_cat_frac']-2.4*liqs['Ca_Liq_cat_frac']))
        logfo2_calc=deltaQFM-25050/T_K+8.58
        liqs['deltaQFM_calc']=deltaQFM
        liqs['logfo2_calc']=logfo2_calc
    if C5==2: #If specify Fe3Fet
        deltaQFM=logfo2-25050/T_K+8.58
        liqs['deltaQFM_calc']=deltaQFM

    if C5==2: # e.g. if Fe3Fet not given
        Fe2Fet_Liq=1/(1+10**(0.25*deltaQFM-1.36+2*liqs['Na_Liq_cat_frac']+2.4*liqs['Ca_Liq_cat_frac']+3.7*liqs['K_Liq_cat_frac']))
        liqs['Fe2Fet_Liq_calc']=Fe2Fet_Liq
    if C5==1: # IF Fe3Fet is given
        Fe2Fet_Liq=1-Fe3Fet_Liq

    liqs['Fe2_Liq_cat_frac']=liqs['Fet_Liq_cat_frac']*Fe2Fet_Liq

    liqs['LnCS2_calc']=(8.77-23590/T_K+(1673/T_K)*(6.7*(liqs['Na_Liq_cat_frac']+liqs['K_Liq_cat_frac'])
        +4.9*liqs['Mg_Liq_cat_frac']+8.1*liqs['Ca_Liq_cat_frac']+8.9*(liqs['Fet_Liq_cat_frac']+liqs['Mn_Liq_cat_frac'])
        +5*liqs['Ti_Liq_cat_frac']+1.8*liqs['Al_Liq_cat_frac']
        -22.2*liqs['Ti_Liq_cat_frac']*(liqs['Fet_Liq_cat_frac']+liqs['Mn_Liq_cat_frac'])
            +7.2*((liqs['Fet_Liq_cat_frac']+liqs['Mn_Liq_cat_frac'])*liqs['Si_Liq_cat_frac']))-2.06*erf(-7.2*(liqs['Fet_Liq_cat_frac']+liqs['Mn_Liq_cat_frac'])))

    liqs['LnCS6_calc']=calculate_BW2022_CS6(df=df, T_K=T_K).LnCS6_calc_OM22_format

    liqs['LnKSO2S2']=-55921/T_K+25.07-0.6465*np.log(T_K)

    liqs['LnS6S2']=(liqs['LnCS6_calc']-liqs['LnKSO2S2']-liqs['LnCS2_calc'])+2*np.log(10)*logfo2_calc
    liqs['S6St_Liq']=1-1/(1+np.exp(liqs['LnS6S2']))

    cols_to_move = ['S6St_Liq', 'LnCS2_calc', 'LnCS6_calc', 'LnKSO2S2', 'LnS6S2',
                    'deltaQFM_calc']
    liqs = liqs[cols_to_move +
                                    [col for col in liqs.columns if col not in cols_to_move]]
    return liqs


def calculate_fo2_QFM_buffers(logfo2=-8.52, T_K=1200, P_kbar=1):
    """ Calculates fo2 values for Frost, Oneill, MELTS and Petrolog3/Myers 'QFM' buffer positions, and
    the offset from each QFM bufer

    Parameters
    ----------------
    logfo2: float, int, pd.Series, np.array
        logfo2 value

    T_K: float, int, pd.Series, np.array
        Temperature in Kelvin

    P_kbar: float, int, pd.Series, np.array
        Presure in kbar

    Returns
    -------------------
    pd.DataFrame
        logfo2 positions for each buffer, and the delta QFM values for each reference point.
    """
    logfo2_QFM_Oneill=8.58-25050/T_K
    # Frost depends on temperature

    logfo2_QFM_highT=(-25096.3/T_K) + 8.735 + 0.11 * ((P_kbar*1000)-1)/T_K
    T_Choice='HighT Beta Qtz'

    logfo2_QFM_lowT=(-26455.3/T_K) +10.344 + 0.092 * ((P_kbar*1000)-1)/T_K
    T_Choice='Low T alpha Qtz'

    Cut_off_T=573+273.15+0.025*(P_kbar*1000)
    if T_K<Cut_off_T:
        logfo2_QFM_Frost= logfo2_QFM_lowT
    if T_K>=Cut_off_T:
        logfo2_QFM_Frost=logfo2_QFM_highT


    logfo2_QFM_Petrolog=-24442/T_K+8.29 #This is actually Myers 1983
    logfo2_QFM_MeltsExcel=-24442/T_K + 8.29 + 0.111*(P_kbar*1000-1)/T_K

    DeltaQFM_Oneill=logfo2-logfo2_QFM_Oneill
    DeltaQFM_Frost=logfo2-logfo2_QFM_Frost
    DeltaQFM_Petrolog=logfo2-logfo2_QFM_Petrolog
    DeltaQFM_MeltsExcel=logfo2-logfo2_QFM_MeltsExcel



    if isinstance(DeltaQFM_Oneill, float):

        df=pd.DataFrame(data={'logfo2_QFM_ONeill': logfo2_QFM_Oneill,
    'logfo2_QFM_Frost': logfo2_QFM_Frost,
    'logfo2_QFM_Petrolog3': logfo2_QFM_Petrolog,
'DeltaQFM_ONeill': DeltaQFM_Oneill,
    'DeltaQFM_Frost': DeltaQFM_Frost,
    'DeltaQFM_Petrolog3': DeltaQFM_Petrolog,
    'DeltaQFM_MeltsExcel': DeltaQFM_MeltsExcel,
    },
    index=[0])
    else:
        df=pd.DataFrame(data={'logfo2_QFM_ONeill': logfo2_QFM_Oneill,
    'logfo2_QFM_Frost': logfo2_QFM_Frost,
    'logfo2_QFM_Petrolog3': logfo2_QFM_Petrolog,
'DeltaQFM_ONeill': DeltaQFM_Oneill,
    'DeltaQFM_Frost': DeltaQFM_Frost,
    'DeltaQFM_Petrolog3': DeltaQFM_Petrolog,
    'DeltaQFM_MeltsExcel': DeltaQFM_MeltsExcel,
    })

    return df




