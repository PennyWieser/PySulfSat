
import numpy as np
import matplotlib.pyplot as plt
from functools import partial
import inspect
import warnings as w
import numbers
import pandas as pd
from PySulfSat.core_calcs import *




def return_cali_dataset(model=None):
    """
    This function returns the calibration dataset for different models.
    This allows you to make your own plots rather than using the generic_cali_plot() option.


    Parameters
    -------


    model: str
        SCSS models:
        LiZhang2022_SCSS: Li and Zhang (2022)
        B2021_SCSS: Blanchard et al (2021) -
        O2021_SCSS: ONeill (2021) - No dataset yet (said zoom in february)
        Liu2021_SCSS: Liu et al. (2021)
        S2017_SCSS:  Smythe et al. (2017)
        F2015_SCSS: Fontin et al. (2015)

        SCAS models:
        CD2019_SCAS: Chowdhury and Dasgupta (2019)
        ZT2022_SCAS: Zajacz and Tsay (2019)
        MK_2015: Masotta and Keppler(2015)


        CS6:
        O2022_CS6:  ONeill and Mavrogenes (2022)
        BW2022_CS6: Boulliung and Wood (2022)

    """

    #SCSS models
    if model=="S2017_SCSS":
        with open(PySulfSat_dir/'Cali_Smythe17.pkl', 'rb') as f:
            Cali_input=load(f)

    if model=="LiZhang2022_SCSS":
        with open(PySulfSat_dir/'Cali_LiZhang22.pkl', 'rb') as f:
            Cali_input=load(f)

    if model=="Liu2021_SCSS":
        with open(PySulfSat_dir/'Cali_Liu2021.pkl', 'rb') as f:
            Cali_input=load(f)

    if model== "B2021_SCSS":
        with open(PySulfSat_dir/'Cali_Blanchard2021.pkl', 'rb') as f:
            Cali_input=load(f)

    if model=="O2021_SCSS":
        print('waiting on a dataset from Hugh!')
        # with open(PySulfSat_dir/'", 'rb') as f:
        #     Cali_input=load(f)

    if model=="F2015_SCSS":
        with open(PySulfSat_dir/'Cali_Fortin2015.pkl', 'rb') as f:
            Cali_input=load(f)

    #SCAS models
    if model=="CD2019_SCAS":
        with open(PySulfSat_dir/'Cali_ChowDas22.pkl', 'rb') as f:
            Cali_input=load(f)

    if model=="ZT2022_SCAS":
        with open(PySulfSat_dir/'Cali_ZT2019.pkl', 'rb') as f:
            Cali_input=load(f)

    if model=="MK_2015":
        with open(PySulfSat_dir/'Cali_MK2015.pkl', 'rb') as f:
            Cali_input=load(f)

    #CS6

    if model=="O2022_CS6":
        with open(PySulfSat_dir/'Cali_OM2022.pkl', 'rb') as f:
            Cali_input=load(f)

    if model=="BW2022_CS6":
        with open(PySulfSat_dir/'Cali_BW2022.pkl', 'rb') as f:
            Cali_input=load(f)

    return Cali_input