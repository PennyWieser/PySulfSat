import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


st_ratio=1/22.6436 # From https://doi.org/10.1016/S0016-7037(01)00611-1


def convert_d34_to_3432S(d34S, st_ratio=1/22.6436):
    """ Converts d34S to 3432S using a st ratio (34/32S value)

    Parameters
    --------------
    d34S: int, float, pd.Series, np.array
        d34S

    st_ratio: S isotope ratio of standard to reference too

    Returns
    -------------
    34/32S ratio.
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


def calculate_S_isotope_factors(*, T_K, S6St_Liq=None):
    """ Calculates S fractionation factor as a function of temperature. If S6St_Liq is entered,
    the fractionaton factor for S2 and SO4 are combined to give one for the overall melt.

    All calcs use FeS_H2S fractionatoin factor from Ohmoto and Rye 1979
    Other fractionation factors from:
        Miyoshi et al. (1984), _M, https://doi.org/10.2343/geochemj.18.75
        or
        Fiege et al. (2015), _F, /10.1016/j.chemgeo.2014.11.012.

    Parameters
    -------------------
    T_K: int, float, pd.Series, np.array
        Temperature in Kelvin
    S6st_Liq: optional, int, float, pd.Series, np.array
        S6/St ratio in liquid, used to calculate a total FeS-S in melt fractionation factor.

    Returns
    -----------------
    df: pd.DataFrame of different fractionatoin factors
    '..._OR79': fractionation factor from Ohmoto and Rye 1979
    '..._M84' Fractionaton factor from Miyoshi et al. (1984)
    '..._F15' Fracionation factor from Fiege et al. (2015)
    if S6ST_Liq is not None:
    'a_FeS_ST_F15_M84' Sulfide melt fractionatoin factor using a_FeS_S2_F15 and a_FeS_SO4_M84 (Fiege doesnt have SO4-FeS)
    'a_FeS_ST_M84' Sulfide melt fractionatoin factor using a_FeS_S2_M84 and a_FeS_SO4_M84




    """
    # if (isinstance(S6St_Liq, float) or instance(S6St_Liq, int)) and (isinstance(T_K, float) or isinstance(T_K, int)):
    #     type='all floats/integers'
    # elif (isinstance(S6St_Liq, float) or instance(S6St_Liq, int)) and  ~(isinstance(T_K, float) or isinstance(T_K, int)):
    #     # T isnt a float or integer, but S is, so need to scale up S6St to be the same length
    #
    #
    if isinstance(T_K, (int, float)) and isinstance(S6St_Liq, (pd.Series, np.ndarray)):
        #print('replacing T with series')
        T_K = pd.Series(np.full(S6St_Liq.shape, T_K))
    elif isinstance(S6St_Liq, (int, float)) and isinstance(T_K, (pd.Series, np.ndarray)):
        #print('replacing S6 with pd.Series')
        S6St_Liq = pd.Series(np.full(T_K.shape, S6St_Liq))


    T=T_K
    T_C=T_K-273.15

    lna_FeS_H2S_1000=0.1*(1000/T)**3 # From Ohmoto and Rye 1979
    lna_S2_SO4_1000_M=-7.4*(1000/T)**2+0.19 # from Miyoshi 1984
    lna_H2S_SO4_1000_M=-6.5*(1000/T)**2+0.19 # from Miyoshi 1984
    lna_H2S_S2_1000_M=lna_H2S_SO4_1000_M-lna_S2_SO4_1000_M
    lna_FeS_S2_1000_M=lna_FeS_H2S_1000+lna_H2S_S2_1000_M
    lna_FeS_SO4_1000_M=lna_FeS_H2S_1000+lna_H2S_SO4_1000_M

    # From Fiege
    #lna_H2S_SO4_1000_F=-6.5*(1000/T)**2+0.19# Fiege doesnt have this term - its from Miyoshi

    lna_H2S_S2_1000_F=10.84*(1000/T)**2-2.5 # From Fiege
    lna_S2_SO4_1000_F=lna_H2S_SO4_1000_M-lna_H2S_S2_1000_F # from Fiege

    lna_FeS_S2_1000_F=lna_FeS_H2S_1000+lna_H2S_S2_1000_F
    #lna_FeS_SO4_1000_F=lna_FeS_H2S_1000+lna_H2S_SO4_1000_M

    if isinstance(T_K, int) or isinstance(T_K, float):

        df=pd.DataFrame(data={'T_K': T_K,
                              'T_C': T_C,
                              'lna_FeS_H2S_1000_OR79': lna_FeS_H2S_1000,
                              'lna_S2_SO4_1000_M84':  lna_S2_SO4_1000_M,
                              'lna_H2S_SO4_1000_M84': lna_H2S_SO4_1000_M,
                              'lna_H2S_S2_1000_M84': lna_H2S_S2_1000_M,
                              'lna_FeS_S2_1000_M84': lna_FeS_S2_1000_M,
                              'lna_FeS_SO4_1000_M84': lna_FeS_SO4_1000_M,

                              'lna_H2S_S2_1000_F15': lna_H2S_S2_1000_F,
                              'lna_S2_SO4_1000_F15': lna_S2_SO4_1000_F,
                              'lna_FeS_S2_1000_F15': lna_FeS_S2_1000_F,
                              }, index=[0])
    else:
            df=pd.DataFrame(data={'T_K': T_K,
                              'T_C': T_C,
                              'lna_FeS_H2S_1000_OR79': lna_FeS_H2S_1000,
                              'lna_S2_SO4_1000_M84':  lna_S2_SO4_1000_M,
                              'lna_H2S_SO4_1000_M84': lna_H2S_SO4_1000_M,
                              'lna_H2S_S2_1000_M84': lna_H2S_S2_1000_M,
                              'lna_FeS_S2_1000_M84': lna_FeS_S2_1000_M,
                              'lna_FeS_SO4_1000_M84': lna_FeS_SO4_1000_M,

                              'lna_H2S_S2_1000_F15': lna_H2S_S2_1000_F,
                              'lna_S2_SO4_1000_F15': lna_S2_SO4_1000_F,
                              'lna_FeS_S2_1000_F15': lna_FeS_S2_1000_F,
                              })

    df['a_FeS_S2_F15']=np.exp(df['lna_FeS_S2_1000_F15']/1000)
    #df['a_FeS_SO4_F']=np.exp(df['lna_FeS_SO4_1000_F']/1000)
    df['a_FeS_S2_M84']=np.exp(df['lna_FeS_S2_1000_M84']/1000)
    df['a_FeS_SO4_M84']=np.exp(df['lna_FeS_SO4_1000_M84']/1000)
    df['a_S2_SO4_M84']=np.exp(df['lna_S2_SO4_1000_M84']/1000)
    df['a_S2_SO4_F15']=np.exp(df['lna_S2_SO4_1000_F15']/1000)
    df['a_FeS_H2S_OR79']=np.exp(df['lna_FeS_H2S_1000_OR79']/1000)
    df['a_H2S_S2_M84']=np.exp(df['lna_H2S_S2_1000_M84']/1000)
    df['a_H2S_S2_F15']=np.exp(df['lna_H2S_S2_1000_F15']/1000)
    df['a_H2S_SO4_M84']=np.exp(df['lna_H2S_SO4_1000_M84']/1000)


    if S6St_Liq is not None:
        Sprop=S6St_Liq
        if len(df)==1:
            df['a_FeS_ST_F15_M84']=Sprop*df['a_FeS_SO4_M84'].iloc[0]+ (1-Sprop)*df['a_FeS_S2_F15'].iloc[0]
            df['a_FeS_ST_M84']=Sprop*df['a_FeS_SO4_M84'].iloc[0]+ (1-Sprop)*df['a_FeS_S2_M84'].iloc[0]

        else:
            df['a_FeS_ST_F15_M84']=Sprop*df['a_FeS_SO4_M84']+ (1-Sprop)*df['a_FeS_S2_F15']
            df['a_FeS_ST_M84']=Sprop*df['a_FeS_SO4_M84']+ (1-Sprop)*df['a_FeS_S2_M84']


    return df


