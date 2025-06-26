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

def calculate_total_fractionation(lna_S2_SO4_1000, lna_sulf_S2_1000, S6ST):
    lna_sulf_silicate=S6ST*lna_S2_SO4_1000 + lna_sulf_S2_1000

    a_sulf_sil=np.exp(lna_sulf_silicate / 1000)


    return a_sulf_sil


def calculate_S_isotope_factors(*, T_K, S6ST=0):
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
    S6St: optional, int, float, pd.Series, np.array
        S6/St ratio in liquid, used to calculate a total FeS-S in melt fractionation factor.

    Returns
    -----------------
    df: pd.DataFrame of different fractionatoin factors
    '..._OR79': fractionation factor from Ohmoto and Rye 1979
    '..._M84' Fractionaton factor from Miyoshi et al. (1984)
    '..._F15' Fracionation factor from Fiege et al. (2015)





    """

    if isinstance(T_K, (int, float)) and isinstance(S6ST, (pd.Series, np.ndarray)):
        #print('replacing T with series')
        T_K = pd.Series(np.full(S6ST.shape, T_K))
    elif isinstance(S6ST, (int, float)) and isinstance(T_K, (pd.Series, np.ndarray)):
        #print('replacing S6 with pd.Series')
        S6St = pd.Series(np.full(T_K.shape, S6ST))


    T=T_K
    T_C=T_K-273.15

    # ----------- Sulfide - H2S fractionation factors
    lna_FeS_H2S_1000_OR79=0.1*(1000/T)**2 # Original Ohmoto and Rye (1979) expression, taken from Marini et al. (2011) as cant find original paper.

    lna_FeS_H2S_1000_LL06=0.25*(1000/T)**2 # From Liu and Li (2006), update based on their models

    lna_FeCuS2_H2S_1000_LL06=0.05*(1000/T)**2 # From Liu and Li (2006) based on chalcopyrite

    lna_Cu5FeS4_H2S_1000_LL06=-0.07*(1000/T)**2 # From Liu and Li (2006) based on bornite

    lna_CuFe2S3_H2S_1000_LL06=0.04*(1000/T)**2 # From Liu and Li (2006) based on cubanite

    lna_FeNi2S4_H2S_1000_LL06=-0.27*(1000/T)**2 # From Liu and Li (2006) based on cubanite

    #------------------Now lets get the H2S - S2 fractionation factor from Fiege et al (2015)----------------------

    lna_H2S_S2_1000_F15=10.84*(1000/T)**2-2.5 # From Fiege et al. 2015, page 52, Eq 8



    # This allows us to calculate sulfide-S2 fractionation factors using all of the equations above.

    lna_FeS_S2_1000_OR79_F15=lna_FeS_H2S_1000_OR79+lna_H2S_S2_1000_F15  # using ohmoto and Rye

    lna_FeS_S2_1000_LL06_F15=lna_FeS_H2S_1000_LL06+lna_H2S_S2_1000_F15  # using new Liu and Li

    lna_FeCuS2_S2_1000_LL06_F15=lna_FeCuS2_H2S_1000_LL06+lna_H2S_S2_1000_F15  # using new Liu and Li

    lna_Cu5FeS4_S2_1000_LL06_F15=lna_Cu5FeS4_H2S_1000_LL06+lna_H2S_S2_1000_F15  # using new Liu and Li

    lna_CuFe2S3_S2_1000_LL06_F15=lna_CuFe2S3_H2S_1000_LL06+lna_H2S_S2_1000_F15  # using new Liu and Li

    lna_FeNi2S4_S2_1000_LL06_F15=lna_FeNi2S4_H2S_1000_LL06+lna_H2S_S2_1000_F15  # using new Liu and Li

    # Alternative approach uses Miyoshi's two expressions to get the S2-H2S fractionation factor.

    lna_S2_SO4_1000_M84=-7.4 * 10**6/(T**2)+0.19 # this is from the Miyoshi 1984 paper  - inverted
    lna_H2S_SO4_1000_M84=-6.5*10**6/T**2 # This is from Miyoshi 1984 directly, inverted to get other way up with a minus sign.
    lna_H2S_S2_1000_M84=lna_H2S_SO4_1000_M84-lna_S2_SO4_1000_M84

    # Using some of Fiege, some of Miyoshi - but the big difference is the H2S-S2 term from Fiege.
    lna_S2_SO4_1000_M84_F15=lna_H2S_SO4_1000_M84-lna_H2S_S2_1000_F15



    lna_FeS_S2_1000_OR79_M84=lna_FeS_H2S_1000_OR79+lna_H2S_S2_1000_M84  # using ohmoto and Rye

    lna_FeS_S2_1000_LL06_M84=lna_FeS_H2S_1000_LL06+lna_H2S_S2_1000_M84  # using new Liu and Li

    lna_FeCuS2_S2_1000_LL06_M84=lna_FeCuS2_H2S_1000_LL06+lna_H2S_S2_1000_M84  # using new Liu and Li

    lna_Cu5FeS4_S2_1000_LL06_M84=lna_Cu5FeS4_H2S_1000_LL06+lna_H2S_S2_1000_M84  # using new Liu and Li

    lna_CuFe2S3_S2_1000_LL06_M84=lna_CuFe2S3_H2S_1000_LL06+lna_H2S_S2_1000_M84  # using new Liu and Li

    lna_FeNi2S4_S2_1000_LL06_M84=lna_FeNi2S4_H2S_1000_LL06+lna_H2S_S2_1000_M84  # using new Liu and Li

    # An alternative way still would use the Ohmoto and Lsaga (1982) SO4-H2S approach




    var_names = [
        'T_K', 'T_C','S6ST',
        'lna_FeS_H2S_1000_OR79', 'lna_FeS_H2S_1000_LL06',
        'lna_FeCuS2_H2S_1000_LL06', 'lna_Cu5FeS4_H2S_1000_LL06',
        'lna_CuFe2S3_H2S_1000_LL06', 'lna_FeNi2S4_H2S_1000_LL06',
        'lna_H2S_S2_1000_F15',

        'lna_FeS_S2_1000_OR79_F15', 'lna_FeS_S2_1000_LL06_F15',
        'lna_FeCuS2_S2_1000_LL06_F15', 'lna_Cu5FeS4_S2_1000_LL06_F15',
        'lna_CuFe2S3_S2_1000_LL06_F15', 'lna_FeNi2S4_S2_1000_LL06_F15',


        'lna_H2S_S2_1000_M84', 'lna_H2S_SO4_1000_M84', 'lna_S2_SO4_1000_M84','lna_S2_SO4_1000_M84_F15',


        'lna_FeS_S2_1000_OR79_M84', 'lna_FeS_S2_1000_LL06_M84',
        'lna_FeCuS2_S2_1000_LL06_M84', 'lna_Cu5FeS4_S2_1000_LL06_M84',
        'lna_CuFe2S3_S2_1000_LL06_M84', 'lna_FeNi2S4_S2_1000_LL06_M84'
    ]

    # Build the dictionary using locals()
    local_vars = locals()
    missing = [var for var in var_names if var not in local_vars]
    if missing:
        raise KeyError(f"The following variables were not defined: {missing}")
    data_dict = {var: local_vars[var] for var in var_names}


    # If T_K is scalar, wrap in a DataFrame with index=0
    if isinstance(T_K, (int, float)):
        df = pd.DataFrame(data_dict, index=[0])
    else:
        # For array-style T_K, assume all variables are vectors of the same length
        df = pd.DataFrame(data_dict)

    # Now lets get them all into true fractionation factors.

    for col in df.columns:
        if col.startswith('lna_'):
            new_col = 'a_' + col[4:]  # Remove 'lna_' prefix, add 'a_'
            df[new_col] = np.exp(df[col] / 1000)

    # Now lets calculate the sulfide-total fractionation

    #---------------------------------------------------All ones using Fiege for the H2S-S2 equation-------------------------------------
    # This uses Miyoshi alone for S2_SO4, and then the OR79 and F15 for sulf-sulfide
    df['a_FeS_melt_M84_OR79_F15'] = calculate_total_fractionation(lna_S2_SO4_1000_M84, lna_FeS_S2_1000_OR79_F15, S6ST)

    # Miyoshi for S2_S04, LL06 and F15 for sulf-sulfide
    df['a_FeS_melt_M84_LL06_F15'] = calculate_total_fractionation(lna_S2_SO4_1000_M84, lna_FeS_S2_1000_LL06_F15, S6ST)

    # All the other sulfide compositions from Li and Liu
    df['a_FeCuS2_melt_M84_LL06_F15'] = calculate_total_fractionation(lna_S2_SO4_1000_M84, lna_FeCuS2_S2_1000_LL06_F15, S6ST)
    df['a_FeCuS2_melt_M84_LL06_F15'] = calculate_total_fractionation(lna_S2_SO4_1000_M84, lna_FeCuS2_S2_1000_LL06_F15, S6ST)
    df['a_Cu5FeS4_melt_M84_LL06_F15'] = calculate_total_fractionation(lna_S2_SO4_1000_M84, lna_Cu5FeS4_S2_1000_LL06_F15, S6ST)
    df['a_CuFe2S3_melt_M84_LL06_F15'] = calculate_total_fractionation(lna_S2_SO4_1000_M84, lna_CuFe2S3_S2_1000_LL06_F15, S6ST)
    df['a_FeNi2S4_melt_M84_LL06_F15'] = calculate_total_fractionation(lna_S2_SO4_1000_M84, lna_FeNi2S4_S2_1000_LL06_F15, S6ST)

    #--------------------------All ones using Miyoshi for the H2S-S2 equation-------------------------------------

    df['a_FeS_melt_M84_F15_OR79_F15'] = calculate_total_fractionation(lna_S2_SO4_1000_M84_F15, lna_FeS_S2_1000_OR79_F15, S6ST)

    # Miyoshi for S2_S04, LL06 and F15 for sulf-sulfide
    df['a_FeS_melt_M84_F15_LL06_F15'] = calculate_total_fractionation(lna_S2_SO4_1000_M84_F15, lna_FeS_S2_1000_LL06_F15, S6ST)


    df['a_FeCuS2_melt_M84_F15_LL06_F15'] = calculate_total_fractionation(lna_S2_SO4_1000_M84_F15, lna_FeCuS2_S2_1000_LL06_F15, S6ST)
    df['a_FeCuS2_melt_M84_F15_LL06_F15'] = calculate_total_fractionation(lna_S2_SO4_1000_M84_F15, lna_FeCuS2_S2_1000_LL06_F15, S6ST)
    df['a_Cu5FeS4_melt_M84_F15_LL06_F15'] = calculate_total_fractionation(lna_S2_SO4_1000_M84_F15, lna_Cu5FeS4_S2_1000_LL06_F15, S6ST)
    df['a_CuFe2S3_melt_M84_F15_LL06_F15'] = calculate_total_fractionation(lna_S2_SO4_1000_M84_F15, lna_CuFe2S3_S2_1000_LL06_F15, S6ST)
    df['a_FeNi2S4_melt_M84_F15_LL06_F15'] = calculate_total_fractionation(lna_S2_SO4_1000_M84_F15, lna_FeNi2S4_S2_1000_LL06_F15, S6ST)




    df['a_FeS_melt_M84_OR79_M84'] = calculate_total_fractionation(lna_S2_SO4_1000_M84, lna_FeS_S2_1000_OR79_M84, S6ST)

    # Miyoshi for S2_S04, LL06 and F15 for sulf-sulfide
    df['a_FeS_melt_M84_LL06_M84'] = calculate_total_fractionation(lna_S2_SO4_1000_M84, lna_FeS_S2_1000_LL06_M84, S6ST)

    # All the other sulfide compositions from Li and Liu
    df['a_FeCuS2_melt_M84_LL06_M84'] = calculate_total_fractionation(lna_S2_SO4_1000_M84, lna_FeCuS2_S2_1000_LL06_M84, S6ST)
    df['a_FeCuS2_melt_M84_LL06_M84'] = calculate_total_fractionation(lna_S2_SO4_1000_M84, lna_FeCuS2_S2_1000_LL06_M84, S6ST)
    df['a_Cu5FeS4_melt_M84_LL06_M84'] = calculate_total_fractionation(lna_S2_SO4_1000_M84, lna_Cu5FeS4_S2_1000_LL06_M84, S6ST)
    df['a_CuFe2S3_melt_M84_LL06_M84'] = calculate_total_fractionation(lna_S2_SO4_1000_M84, lna_CuFe2S3_S2_1000_LL06_M84, S6ST)
    df['a_FeNi2S4_melt_M84_LL06_M84'] = calculate_total_fractionation(lna_S2_SO4_1000_M84, lna_FeNi2S4_S2_1000_LL06_M84, S6ST)

    # --------------- All ones using combination of Miyoshi and Fiege for the S2-SO4 term

    df['a_FeS_melt_M84_F15_OR79_M84'] = calculate_total_fractionation(lna_S2_SO4_1000_M84_F15, lna_FeS_S2_1000_OR79_M84, S6ST)

    # Miyoshi for S2_S04, LL06 and F15 for sulf-sulfide
    df['a_FeS_melt_M84_F15_LL06_M84'] = calculate_total_fractionation(lna_S2_SO4_1000_M84_F15, lna_FeS_S2_1000_LL06_M84, S6ST)

    # All the other sulfide compositions from Li and Liu
    df['a_FeCuS2_melt_M84_F15_LL06_M84'] = calculate_total_fractionation(lna_S2_SO4_1000_M84_F15, lna_FeCuS2_S2_1000_LL06_M84, S6ST)
    df['a_FeCuS2_melt_M84_F15_LL06_M84'] = calculate_total_fractionation(lna_S2_SO4_1000_M84_F15, lna_FeCuS2_S2_1000_LL06_M84, S6ST)
    df['a_Cu5FeS4_melt_M84_F15_LL06_M84'] = calculate_total_fractionation(lna_S2_SO4_1000_M84_F15, lna_Cu5FeS4_S2_1000_LL06_M84, S6ST)
    df['a_CuFe2S3_melt_M84_F15_LL06_M84'] = calculate_total_fractionation(lna_S2_SO4_1000_M84_F15, lna_CuFe2S3_S2_1000_LL06_M84, S6ST)
    df['a_FeNi2S4_melt_M84_F15_LL06_M84'] = calculate_total_fractionation(lna_S2_SO4_1000_M84_F15, lna_FeNi2S4_S2_1000_LL06_M84, S6ST)





    return df


