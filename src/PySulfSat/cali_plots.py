
import numpy as np
import matplotlib.pyplot as plt
from functools import partial
import inspect
import warnings as w
import numbers
import pandas as pd
from PySulfSat.core_calcs import *
import os




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

    Returns
    -------
    pd.DataFrame
        dataframe of the calibration dataset, with column headings MgO_Liq, SiO2_Liq etc.


    """
    subfolder_name='Cali_datasets_CSV'
    sup_models = ['S2017_SCSS', 'LiZhang2022_SCSS', 'Liu2021_SCSS', 'B2021_SCSS', 'O2021_SCSS', 'F2015_SCSS',
    'CD2019_SCAS', 'ZT2022_SCAS', 'MK_2015',
    'O2022_CS6', 'BW2022_CS6']



    # Check if the provided model is in the list of supported models
    if model not in sup_models:
        raise ValueError(f"Invalid model: {model}. Supported models with inbuilt calibration datasets are {', '.join(sup_models)}")






    #SCSS models
    if model=="S2017_SCSS":
        csv_file_path = os.path.join(PySulfSat_dir, subfolder_name, 'Cali_Smythe17.csv')
        Cali_input=pd.read_csv(csv_file_path)



    if model=="LiZhang2022_SCSS":
        csv_file_path = os.path.join(PySulfSat_dir, subfolder_name, 'Cali_LiZhang22.csv')
        Cali_input=pd.read_csv(csv_file_path)

    if model=="Liu2021_SCSS":

        csv_file_path = os.path.join(PySulfSat_dir, subfolder_name, 'Cali_Liu2021.csv')
        Cali_input=pd.read_csv(csv_file_path)

    if model== "B2021_SCSS":
        csv_file_path = os.path.join(PySulfSat_dir, subfolder_name, 'Cali_Blanchard2021.csv')
        Cali_input=pd.read_csv(csv_file_path)


    if model=="O2021_SCSS":
        raise TypeError('waiting on a dataset from Hugh!')
        # with open(PySulfSat_dir/'", 'rb') as f:
        #     Cali_input=load(f)

    if model=="F2015_SCSS":
        csv_file_path = os.path.join(PySulfSat_dir, subfolder_name, 'Cali_Fortin2015.csv')
        Cali_input=pd.read_csv(csv_file_path)

    #SCAS models
    if model=="CD2019_SCAS":
        csv_file_path = os.path.join(PySulfSat_dir, subfolder_name, 'Cali_ChowDas22.csv')
        Cali_input=pd.read_csv(csv_file_path)

    if model=="ZT2022_SCAS":
        csv_file_path = os.path.join(PySulfSat_dir, subfolder_name, 'Cali_ZT2019.csv')
        Cali_input=pd.read_csv(csv_file_path)

    if model=="MK_2015":
        csv_file_path = os.path.join(PySulfSat_dir, subfolder_name, 'Cali_MK2015.csv')
        Cali_input=pd.read_csv(csv_file_path)

    #CS6

    if model=="O2022_CS6":
        csv_file_path = os.path.join(PySulfSat_dir, subfolder_name, 'Cali_OM2022.csv')
        Cali_input=pd.read_csv(csv_file_path)

    if model=="BW2022_CS6":
        csv_file_path = os.path.join(PySulfSat_dir, subfolder_name, 'Cali_BW2022.csv')
        Cali_input=pd.read_csv(csv_file_path)

    return Cali_input