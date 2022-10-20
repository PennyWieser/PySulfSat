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

def calculate_S6St_Jugo2010(DeltaQFM):
    '''
    Calculates The S6/St ratio from DeltaQFM using
    Jugo et al. (2010)

    Parameters
    -----------
    DeltaQFM: int, float, pandas.Series
        Offset from the QFM buffer in log units.

    Returns
    -----------
    float, pandas.Series
        Proportion of S6+ in the liquid

    '''
    S6St_Liq=1/(1+10**(2.1-2*DeltaQFM))

    return S6St_Liq

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

    log_S6S2=8*log_Fe3Fe2 + 8.7436*10**6/T_K**2 - 27703/T_K + 20.273
    S6S2=10**(log_S6S2)

    S6=S6S2/(1+S6S2)

    S2=1/(1+S6S2)
    S6St=S6/(S6+S2)

    return S6St


def calculate_S_Tot_Kleinsasser2022_dacite(*, SCSS2=None,  SCAS=None, deltaFMQ=None):
    """ Calculates SCSS total as a function of S2- and S6+ in dacitic melts,
    following Kleinsasser et al. 2022"""

    if SCSS2 is not None:
        SCSS_Tot=SCSS2*(1+10**(2*deltaFMQ-3.05))
        return SCSS_Tot
    if SCAS is not None:
        SCAS_Tot=SCAS*(1+np.exp(1.26-2*deltaFMQ))
        return SCAS_Tot



def calculate_S_Tot_Jugo2010(*, SCSS2=None,  SCAS=None,DeltaQFM):
    '''Calculates S total as a function of S2- and S6+
    following Jugo et al. 2010 for given SCSS or SCAS

    Parameters
    -----------
    DeltaQFM: int, float, pandas.Series
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




