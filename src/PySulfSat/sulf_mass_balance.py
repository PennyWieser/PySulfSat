import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


st_ratio=1/22.6436 # From https://doi.org/10.1016/S0016-7037(01)00611-1


def convert_d34_to_3432S(d34S, st_ratio=1/22.6436):
    """ Converts d34S to 3432S using a st ratio (34/32S value)
    """
    S3432=((d34S/1000)+1)*(st_ratio)
    return S3432

def convert_3432S_to_d34(S3432, st_ratio=1/22.6436):
    """ Converts d34S to 3432S using a st ratio (34/32S value)
    """
    d34S=(((S3432)/(st_ratio)) -1)*1000
    return d34S

def calculate_std_ratio_used(d34S, S3432=1/22.6436):
    """ When you have both a d34S value and a S3432 value,
    calculates what ratio was used by the study
    """
    S3432_standard=(1000*S3432)/(d34S+1000)

    return S3432_standard


def crystallize_S_incomp(S_init=1200, F_melt=None):
    """
    Calculates amount of S remaining in the melt for different melt fractoins (F_melt)
    assuming S is completely incompatible in silicate minerals


    Parameters
    --------------
    S_init: int, float, pd.Series
        initial S content in ppm
    F_melt: int, float, pd.Series
        melt fraction (between 0 and 1)


    Returns
    --------------
    S in the melt for each F

    """


    S_melt_FC=S_init/F_melt
    return S_melt_FC


def calculate_mass_frac_sulf(S_model=None, S_init=None, F_melt=None, S_sulf=None):
    """
    Calculates mass fraction of sulfide removed for a given SCSS value, a S content in the melt,
    and S content in the sulfide

    Parameters
    --------------
    S_model: int, float, pd.Series
        modelled amount of S present in the melt in ppm (could be SCSS2-, SCSStot, or SCAS, or STot from another method)
    S_init: int, float, pd.Series
        initial S content of the system where F_melt=1 in ppm
    F_melt:  int, float, pd.Series
        melt fraction (between 0 and 1)
    S_sulf:  int, float, pd.Series
        S content of the sulfide in ppm

    Returns
    --------------
    Mass fraction of sulfide

    """


    X_sulf=(S_init-F_melt*S_model)/S_sulf
    return X_sulf
