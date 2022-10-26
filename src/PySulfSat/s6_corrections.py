import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def calculate_SCSS_Total(SCSS, S6St_Liq):
    '''
    Calculates SCSS Total from the SCSS (2-)
    See Wieser et al. (2019) for details

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


def calculate_SCAS_Total(SCAS, S2St_Liq):
    '''
    Calculates SCAS Total from the SCAS (6+)
    See Wieser et al. (2019) for details

    Parameters
    -----------
    SCAS: int, float, pandas.Series
        SCAS from different calculators (only uses S6+)

    S2St_Liq: int, float, pandas.Series
        Proportion of S2- in the liquid

    Returns
    -----------
    float, pandas.Series
        SCAStotal, i.e inputted SCAS (6+) scaled up to account for S present in the S2- form.

    '''
    SCAS_Total=SCAS/(1-S2St_Liq)
    return SCAS_Total

def calculate_S6St_Jugo2010_eq10(deltaQFM):
    '''
    Calculates The S6/St ratio from deltaQFM using
    Jugo et al. (2010) eq10

    Parameters
    -----------
    deltaQFM: int, float, pandas.Series
        Offset from the QFM buffer in log units.

    Returns
    -----------
    float, pandas.Series
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
    deltaQFM: int, float, pandas.Series
        Offset from the QFM buffer in log units.

    Returns
    -----------
    float, pandas.Series
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
    T_K: int, float, pandas.Series
        Temperature in Kelvin.

    Fe3Fet_Liq: int, float, pandas.Series
        Proportion of Fe3 in the liquid
    Returns
    -----------
    float, pandas.Series
        Proportion of S6+ in the liquid

    '''
    # First, calculate Fe3/Fe2 from Fe3/Fet
    Fe3Fe2=Fe3Fet_Liq/(1-Fe3Fet_Liq)
    log_Fe3Fe2=np.log10(Fe3Fe2)

    log_S6S2=8*log_Fe3Fe2 + 8.7436*10**6/(T_K**2) - 27703/T_K + 20.273
    S6S2=10**(log_S6S2)

    S6=S6S2/(1+S6S2)

    S2=1/(1+S6S2)
    S6St=S6/(S6+S2)

    return S6St


def calculate_S_Tot_Kleinsasser2022_dacite(*, SCSS2=None,  SCAS=None, deltaQFM=None):
    """ Calculates SCSS total as a function of S2- and S6+ in dacitic melts,
    following Kleinsasser et al. 2022"""

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
                                T_K=None, Fe3Fet_Liq=None):
    """
    Calculates the total amount of S accounting for the SCSS and SCAS

    Parameters
    -----------
    SCSS: int, float, pandas.Series
        SCSS (2-) you have calculated from whatever model you want.

    SCAS: int, float, pandas.Series
        SCAS (6+) you have calculated from whatever model you want.

    model: str
        Model used to calculate S6St_Liq, choice of:
        'Nash': Uses Nash et al. (2019) based on Fe3Fet_Liq and T_K
        'Jugo': Uses Jugo et al. (2010) Eq 10 based on DeltaQFM
        'Kleinsasser': Uses Kleinsasser et al. 2022 based on deltaQFM

    deltaQFM: int, float, panda.Series (for Jugo et al. 2010, or Kleinsasser et al. 2022)

    Fe3Fet_Liq: int, float, panda.Series
        Proportion of Fe3Fet in the liquid, needed if model='Nash'

    T_K: int, float, panda.Series
        Temperature in Kelvin, needed if model='Nash'


    """
    if model =="Kleinsasser":
        SCSS_Tot=calculate_S_Tot_Kleinsasser2022_dacite(SCSS2=SCSS,
                                                           deltaQFM=deltaQFM)
        SCAS_Tot=calculate_S_Tot_Kleinsasser2022_dacite(SCAS=SCAS,
                                                           deltaQFM=deltaQFM)
        S6_St=np.nan
    else:
        if model =="Jugo":
            S6_St=calculate_S6St_Jugo2010_eq10(deltaQFM=deltaQFM)

        if model =="Nash":
            if T_K is None:
                raise Exception('Need a T_K input to use Nash')
            if Fe3Fet_Liq is None:
                raise Exception('Need a Fe3Fet_Liq input to use Nash')
            S6_St=calculate_S6St_Nash2019(T_K=T_K, Fe3Fet_Liq=Fe3Fet_Liq)

        if S6St_Liq is not None:
            S6_St=S6St_Liq


        SCSS_Tot=calculate_SCSS_Total(SCSS=SCSS, S6St_Liq=S6_St)
        SCAS_Tot=calculate_SCAS_Total(SCAS=SCAS, S2St_Liq=1-S6_St)



    SCSS_S6_cont=SCSS_Tot-SCSS
    SCAS_S2_cont=SCAS_Tot-SCAS

    df_Species=pd.DataFrame(data={'deltaQFM': deltaQFM,
                                  'S6_St': S6_St,
                                  'SCSS_2': SCSS,
                              'SCAS_6': SCAS,
                              'SCSS_Tot': SCSS_Tot,
                              'SCAS_Tot': SCAS_Tot,
                                 'S6 in SCSS_Tot': SCSS_S6_cont,
                               'S2 in SCAS_Tot': SCAS_S2_cont
                              })

    # Cant have more S6 than the SCAS
    toomuchS6=df_Species['S6 in SCSS_Tot']>df_Species['SCAS_6']
    df_Species['SCSS_Tot_check']=df_Species['SCSS_Tot']
    df_Species.loc[toomuchS6, 'SCSS_Tot_check']=df_Species['SCAS_6']+df_Species['SCSS_2']

    # Cant have more S2- than the SCSS
    toomuchS2=df_Species['S2 in SCAS_Tot']>df_Species['SCSS_2']
    df_Species['SCAS_Tot_check']=df_Species['SCAS_Tot']
    df_Species.loc[toomuchS2, 'SCAS_Tot_check']=df_Species['SCAS_6']+df_Species['SCSS_2']

    #Set as the minimum one
    df_Species.insert(0, 'Total_S',df_Species['SCAS_Tot_check'])
    SCAS_higher=df_Species['SCAS_Tot_check']>df_Species['SCSS_Tot_check']
    df_Species.loc[SCAS_higher, 'Total_S']=df_Species['SCSS_Tot_check']

    df_Species.insert(1, 'S2_Tot', df_Species['Total_S']*(1-df_Species['S6_St']))
    df_Species.insert(2, 'S6_Tot', df_Species['Total_S']*(df_Species['S6_St']))

    df_Species2=df_Species.copy()
    df_Species2.drop(['SCAS_Tot_check', 'SCSS_Tot_check', 'SCSS_Tot_check'], axis=1, inplace=True)
    return df_Species2

