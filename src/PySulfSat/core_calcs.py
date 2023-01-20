import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.optimize as optimize
from scipy.special import erf
from pathlib import Path
from pickle import load
import pickle
PySulfSat_dir=Path(__file__).parent





## Anhydrous Molar masses, used for ONeill (2021). Doesnt include P2O5, Cr2O3 or H2O in the calculation.
oxide_mass_liq_anhyd = {'SiO2_Liq': 60.0843, 'MgO_Liq': 40.3044,
'MnO_Liq': 70.9375, 'FeOt_Liq': 71.8464, 'CaO_Liq': 56.0774,
'Al2O3_Liq': 101.961, 'Na2O_Liq': 61.9789, 'K2O_Liq': 94.196,
'TiO2_Liq': 79.8788, 'P2O5_Liq': 141.944522 }


cation_num_liq_anhyd = {'SiO2_Liq': 1, 'MgO_Liq': 1, 'MnO_Liq': 1,
'FeOt_Liq': 1, 'CaO_Liq': 1, 'Al2O3_Liq': 2, 'Na2O_Liq': 2,
'K2O_Liq': 2, 'TiO2_Liq': 1, 'P2O5_Liq': 2}


# Turns dictionary into a dataframe so pandas matrix math functions can be used
oxide_mass_liq_anhyd_df = pd.DataFrame.from_dict(
    oxide_mass_liq_anhyd, orient='index').T
oxide_mass_liq_anhyd_df['Sample_ID_Liq'] = 'MolWt'
oxide_mass_liq_anhyd_df.set_index('Sample_ID_Liq', inplace=True)



#
cation_num_liq_anhyd_df = pd.DataFrame.from_dict(
    cation_num_liq_anhyd, orient='index').T
cation_num_liq_anhyd_df['Sample_ID_Liq'] = 'CatNum'
cation_num_liq_anhyd_df.set_index('Sample_ID_Liq', inplace=True)

## Hydrous Molar masses - Used for Fortin et al. (2015)

oxide_mass_liq_hyd = {'SiO2_Liq': 60.0843, 'MgO_Liq': 40.3044,
'MnO_Liq': 70.9375, 'FeOt_Liq': 71.8464, 'CaO_Liq': 56.0774,
'Al2O3_Liq': 101.961, 'Na2O_Liq': 61.9789, 'K2O_Liq': 94.196,
'TiO2_Liq': 79.8788, 'P2O5_Liq': 141.944, 'H2O_Liq': 18,
'Cr2O3_Liq': 151.99}


cation_num_liq_hyd = {'SiO2_Liq': 1, 'MgO_Liq': 1, 'MnO_Liq': 1,
'FeOt_Liq': 1, 'CaO_Liq': 1, 'Al2O3_Liq': 2, 'Na2O_Liq': 2,
'K2O_Liq': 2, 'TiO2_Liq': 1, 'P2O5_Liq': 2, 'H2O_Liq':2,
'Cr2O3_Liq':2}


# Turns dictionary into a dataframe so pandas matrix math functions can be used
oxide_mass_liq_hyd_df = pd.DataFrame.from_dict(
    oxide_mass_liq_hyd, orient='index').T
oxide_mass_liq_hyd_df['Sample_ID_Liq'] = 'MolWt'
oxide_mass_liq_hyd_df.set_index('Sample_ID_Liq', inplace=True)



#
cation_num_liq_hyd_df = pd.DataFrame.from_dict(
    cation_num_liq_hyd, orient='index').T
cation_num_liq_hyd_df['Sample_ID_Liq'] = 'CatNum'
cation_num_liq_hyd_df.set_index('Sample_ID_Liq', inplace=True)


## With




## Anhydrous proportions, fractions, calculations from Thermobar (Wieser et al. in prep)


def calculate_anhydrous_mol_proportions_liquid(liq_comps):
    '''Import Liq compositions using liq_comps=My_Liquids, returns anhydrous mole proportions

   Parameters
    -------

    liq_comps: pandas.DataFrame
        Panda DataFrame of liquid compositions with column headings SiO2_Liq, TiO2_Liq etc.

    Returns
    -------
    pandas DataFrame
        anhydrous mole proportions for the liquid with column headings of the form SiO2_Liq_mol_prop

    '''
    # This makes the input match the columns in the oxide mass dataframe
    liq_wt = liq_comps.reindex(
        oxide_mass_liq_anhyd_df.columns, axis=1).fillna(0)
    # Combine the molecular weight and weight percent dataframes
    liq_wt_combo = pd.concat([oxide_mass_liq_anhyd_df, liq_wt],)
    # Drop the calculation column
    mol_prop_anhyd = liq_wt_combo.div(
        liq_wt_combo.loc['MolWt', :], axis='columns').drop(['MolWt'])
    mol_prop_anhyd.columns = [
        str(col) + '_mol_prop' for col in mol_prop_anhyd.columns]
    return mol_prop_anhyd


def calculate_anhydrous_mol_fractions_liquid(liq_comps):
    '''Import Liq compositions using liq_comps=My_Liquids, returns anhydrous mole fractions

   Parameters
    -------

    liq_comps: pandas.DataFrame
                Panda DataFrame of liquid compositions with column headings SiO2_Liq, TiO2_Liq etc.

    Returns
    -------
    pandas DataFrame
        anhydrous mole fractions for the liquid with column headings of the form SiO2_Liq_mol_frac

    '''
    mol_prop = calculate_anhydrous_mol_proportions_liquid(liq_comps)
    mol_prop['sum'] = mol_prop.sum(axis='columns')
    mol_frac_anhyd = mol_prop.div(mol_prop['sum'], axis='rows')
    mol_frac_anhyd.drop(['sum'], axis='columns', inplace=True)
    mol_frac_anhyd.columns = [str(col).replace('prop', 'frac')
                              for col in mol_frac_anhyd.columns]
    return mol_frac_anhyd


def calculate_anhydrous_cat_proportions_liquid(liq_comps):
    '''Import Liq compositions using liq_comps=My_Liquids, returns anhydrous cation proportions (e.g., mole proportions * no of cations)

   Parameters
    -------

    liq_comps: pandas.DataFrame
                Panda DataFrame of liquid compositions with column headings SiO2_Liq, TiO2_Liq etc.


    Returns
    -------
    pandas DataFrame
        anhydrous cation proportions for the liquid with column headings of the form S_Liq_cat_prop

    '''
    mol_prop_no_cat_num = calculate_anhydrous_mol_proportions_liquid(liq_comps)
    mol_prop_no_cat_num.columns = [str(col).replace(
        '_mol_prop', '') for col in mol_prop_no_cat_num.columns]
    ox_num_reindex = cation_num_liq_anhyd_df.reindex(
        oxide_mass_liq_anhyd_df.columns, axis=1).fillna(0)
    df_calc_comb = pd.concat([ox_num_reindex, mol_prop_no_cat_num])
    cation_prop_anhyd = df_calc_comb.multiply(
        df_calc_comb.loc['CatNum', :], axis='columns').drop(['CatNum'])
    cation_prop_anhyd.columns = [
        str(col) + '_cat_prop' for col in cation_prop_anhyd.columns]



    return cation_prop_anhyd



def calculate_anhydrous_cat_fractions_liquid(liq_comps):
    '''Import Liq compositions using liq_comps=My_Liquids, returns anhydrous cation fractions

   Parameters
    -------

    liq_comps: pandas.DataFrame
        liquid compositions with column headings SiO2_Liq, TiO2_Liq etc.


    Returns
    -------
    pandas DataFrame
        anhydrous cation fractions for the liquid with column headings of the form _Liq_cat_frac,
        as well as the initial dataframe of liquid compositions.


    '''
    cat_prop = calculate_anhydrous_cat_proportions_liquid(liq_comps=liq_comps)
    mol_prop = calculate_anhydrous_mol_fractions_liquid(liq_comps=liq_comps)
    cat_prop['sum'] = cat_prop.sum(axis='columns')
    cat_frac_anhyd = cat_prop.div(cat_prop['sum'], axis='rows')
    cat_frac_anhyd.drop(['sum'], axis='columns', inplace=True)
    cat_frac_anhyd.columns = [str(col).replace('prop', 'frac')
                              for col in cat_frac_anhyd.columns]
    cat_frac_anhyd = pd.concat([liq_comps, mol_prop, cat_frac_anhyd], axis=1)

    if "Fe3Fet_Liq" in cat_frac_anhyd:
        cat_frac_anhyd['Mg_Number_Liq_NoFe3'] = (cat_frac_anhyd['MgO_Liq'] / 40.3044) / (
            (cat_frac_anhyd['MgO_Liq'] / 40.3044) + (cat_frac_anhyd['FeOt_Liq'] / 71.844))
        cat_frac_anhyd['Mg_Number_Liq_Fe3'] = (cat_frac_anhyd['MgO_Liq'] / 40.3044) / (
            (cat_frac_anhyd['MgO_Liq'] / 40.3044) + (cat_frac_anhyd['FeOt_Liq'] * (1 - cat_frac_anhyd['Fe3Fet_Liq']) / 71.844))
    if "Fe3Fet_Liq" not in cat_frac_anhyd:
        cat_frac_anhyd['Mg_Number_Liq_Fe3'] = (cat_frac_anhyd['MgO_Liq'] / 40.3044) / (
            (cat_frac_anhyd['MgO_Liq'] / 40.3044) + (cat_frac_anhyd['FeOt_Liq'] / 71.844))
        cat_frac_anhyd['Mg_Number_Liq_NoFe3'] = (cat_frac_anhyd['MgO_Liq'] / 40.3044) / (
            (cat_frac_anhyd['MgO_Liq'] / 40.3044) + (cat_frac_anhyd['FeOt_Liq'] / 71.844))


    cation_frac_anhyd2=cat_frac_anhyd.rename(columns={
                        'SiO2_Liq_cat_frac': 'Si_Liq_cat_frac',
                        'TiO2_Liq_cat_frac': 'Ti_Liq_cat_frac',
                        'Al2O3_Liq_cat_frac': 'Al_Liq_cat_frac',
                        'FeOt_Liq_cat_frac': 'Fet_Liq_cat_frac',
                        'MnO_Liq_cat_frac': 'Mn_Liq_cat_frac',
                        'MgO_Liq_cat_frac': 'Mg_Liq_cat_frac',
                        'CaO_Liq_cat_frac': 'Ca_Liq_cat_frac',
                        'Na2O_Liq_cat_frac': 'Na_Liq_cat_frac',
                        'K2O_Liq_cat_frac': 'K_Liq_cat_frac',
                        'Cr2O3_Liq_cat_frac': 'Cr_Liq_cat_frac',
                        'P2O5_Liq_cat_frac': 'P_Liq_cat_frac',

                        })

    return cation_frac_anhyd2


## hydrous proportions, fractions, calculations from Thermobar (Wieser et al. in prep)


def calculate_hydrous_mol_proportions_liquid(liq_comps):
    '''Import Liq compositions using liq_comps=My_Liquids, returns hydrous mole proportions

   Parameters
    -------

    liq_comps: pandas.DataFrame
        Panda DataFrame of liquid compositions with column headings SiO2_Liq, TiO2_Liq etc.

    Returns
    -------
    pandas DataFrame
        hydrous mole proportions for the liquid with column headings of the form SiO2_Liq_mol_prop

    '''
    # This makes the input match the columns in the oxide mass dataframe
    liq_wt = liq_comps.reindex(
        oxide_mass_liq_hyd_df.columns, axis=1).fillna(0)
    # Combine the molecular weight and weight percent dataframes
    liq_wt_combo = pd.concat([oxide_mass_liq_hyd_df, liq_wt],)
    # Drop the calculation column
    mol_prop_hyd = liq_wt_combo.div(
        liq_wt_combo.loc['MolWt', :], axis='columns').drop(['MolWt'])
    mol_prop_hyd.columns = [
        str(col) + '_mol_prop' for col in mol_prop_hyd.columns]
    return mol_prop_hyd


def calculate_hydrous_mol_fractions_liquid(liq_comps):
    '''Import Liq compositions using liq_comps=My_Liquids, returns hydrous mole fractions

   Parameters
    -------

    liq_comps: pandas.DataFrame
                Panda DataFrame of liquid compositions with column headings SiO2_Liq, TiO2_Liq etc.

    Returns
    -------
    pandas DataFrame
        hydrous mole fractions for the liquid with column headings of the form SiO2_Liq_mol_frac

    '''
    mol_prop = calculate_hydrous_mol_proportions_liquid(liq_comps)
    mol_prop['sum'] = mol_prop.sum(axis='columns')
    mol_frac = mol_prop.div(mol_prop['sum'], axis='rows')
    mol_frac.drop(['sum'], axis='columns', inplace=True)
    mol_frac.columns = [str(col).replace('prop', 'frac')
                              for col in mol_frac.columns]
    return mol_frac


def calculate_hydrous_cat_proportions_liquid(liq_comps):
    '''Import Liq compositions using liq_comps=My_Liquids, returns hydrous cation proportions (e.g., mole proportions * no of cations)

   Parameters
    -------

    liq_comps: pandas.DataFrame
                Panda DataFrame of liquid compositions with column headings SiO2_Liq, TiO2_Liq etc.


    Returns
    -------
    pandas DataFrame
        hydrous cation proportions for the liquid with column headings of the form S_Liq_cat_prop

    '''
    mol_prop_no_cat_num = calculate_hydrous_mol_proportions_liquid(liq_comps)
    mol_prop_no_cat_num.columns = [str(col).replace(
        '_mol_prop', '') for col in mol_prop_no_cat_num.columns]
    ox_num_reindex = cation_num_liq_hyd_df.reindex(
        oxide_mass_liq_hyd_df.columns, axis=1).fillna(0)
    df_calc_comb = pd.concat([ox_num_reindex, mol_prop_no_cat_num])
    cation_prop_hyd = df_calc_comb.multiply(
        df_calc_comb.loc['CatNum', :], axis='columns').drop(['CatNum'])
    cation_prop_hyd.columns = [
        str(col) + '_cat_prop' for col in cation_prop_hyd.columns]



    return cation_prop_hyd



def calculate_hydrous_cat_fractions_liquid(liq_comps):
    '''Import Liq compositions using liq_comps=My_Liquids, returns hydrous cation fractions

   Parameters
    -------

    liq_comps: pandas.DataFrame
        liquid compositions with column headings SiO2_Liq, TiO2_Liq etc.


    Returns
    -------
    pandas DataFrame
        hydrous cation fractions for the liquid with column headings of the form _Liq_cat_frac,
        as well as the initial dataframe of liquid compositions.


    '''
    cat_prop = calculate_hydrous_cat_proportions_liquid(liq_comps=liq_comps)
    mol_prop = calculate_hydrous_mol_fractions_liquid(liq_comps=liq_comps)
    cat_prop['sum'] = cat_prop.sum(axis='columns')
    cat_frac_hyd = cat_prop.div(cat_prop['sum'], axis='rows')
    cat_frac_hyd.drop(['sum'], axis='columns', inplace=True)
    cat_frac_hyd.columns = [str(col).replace('prop', 'frac')
                              for col in cat_frac_hyd.columns]
    cat_frac_hyd = pd.concat([liq_comps, mol_prop, cat_frac_hyd], axis=1)

    if "Fe3Fet_Liq" in cat_frac_hyd:
        cat_frac_hyd['Mg_Number_Liq_NoFe3'] = (cat_frac_hyd['MgO_Liq'] / 40.3044) / (
            (cat_frac_hyd['MgO_Liq'] / 40.3044) + (cat_frac_hyd['FeOt_Liq'] / 71.844))
        cat_frac_hyd['Mg_Number_Liq_Fe3'] = (cat_frac_hyd['MgO_Liq'] / 40.3044) / (
            (cat_frac_hyd['MgO_Liq'] / 40.3044) + (cat_frac_hyd['FeOt_Liq'] * (1 - cat_frac_hyd['Fe3Fet_Liq']) / 71.844))
    if "Fe3Fet_Liq" not in cat_frac_hyd:
        cat_frac_hyd['Mg_Number_Liq_Fe3'] = (cat_frac_hyd['MgO_Liq'] / 40.3044) / (
            (cat_frac_hyd['MgO_Liq'] / 40.3044) + (cat_frac_hyd['FeOt_Liq'] / 71.844))
        cat_frac_hyd['Mg_Number_Liq_NoFe3'] = (cat_frac_hyd['MgO_Liq'] / 40.3044) / (
            (cat_frac_hyd['MgO_Liq'] / 40.3044) + (cat_frac_hyd['FeOt_Liq'] / 71.844))


    cation_frac_hyd2=cat_frac_hyd.rename(columns={
                        'SiO2_Liq_cat_frac': 'Si_Liq_cat_frac',
                        'TiO2_Liq_cat_frac': 'Ti_Liq_cat_frac',
                        'Al2O3_Liq_cat_frac': 'Al_Liq_cat_frac',
                        'FeOt_Liq_cat_frac': 'Fet_Liq_cat_frac',
                        'MnO_Liq_cat_frac': 'Mn_Liq_cat_frac',
                        'MgO_Liq_cat_frac': 'Mg_Liq_cat_frac',
                        'CaO_Liq_cat_frac': 'Ca_Liq_cat_frac',
                        'Na2O_Liq_cat_frac': 'Na_Liq_cat_frac',
                        'K2O_Liq_cat_frac': 'K_Liq_cat_frac',
                        'Cr2O3_Liq_cat_frac': 'Cr_Liq_cat_frac',
                        'P2O5_Liq_cat_frac': 'P_Liq_cat_frac',

                        })

    return cation_frac_hyd2

## Hydrous Calculations - Used for Smythe et al. (2017)
oxide_mass_liq = {'SiO2_Liq': 60.0843, 'MgO_Liq': 40.3044,
'MnO_Liq': 70.9375, 'FeO_Liq': 71.8464, 'Fe2O3_Liq': 159.69, 'CaO_Liq': 56.0774,
'Al2O3_Liq': 101.961, 'Na2O_Liq': 61.9789, 'K2O_Liq': 94.196,
'TiO2_Liq': 79.8788, 'P2O5_Liq': 141.944}

oxide_mass_liq_df = pd.DataFrame.from_dict(
    oxide_mass_liq, orient='index').T
oxide_mass_liq_df['Sample_ID_Liq'] = 'MolWt'
oxide_mass_liq_df.set_index('Sample_ID_Liq', inplace=True)




def calculate_mol_proportions_liquid(liq_comps):
    '''Import Liq compositions using liq_comps=My_Liquids, returns anhydrous mole proportions

   Parameters
    -------

    liq_comps: pandas.DataFrame
        Panda DataFrame of liquid compositions with column headings SiO2_Liq, TiO2_Liq etc.

    Returns
    -------
    pandas DataFrame
        anhydrous mole proportions for the liquid with column headings of the form SiO2_Liq_mol_prop

    '''
    # This makes the input match the columns in the oxide mass dataframe
    liq_wt = liq_comps.reindex(
        oxide_mass_liq_df.columns, axis=1).fillna(0)
    # Combine the molecular weight and weight percent dataframes
    liq_wt_combo = pd.concat([oxide_mass_liq_df, liq_wt],)
    # Drop the calculation column
    mol_prop_anhyd = liq_wt_combo.div(
        liq_wt_combo.loc['MolWt', :], axis='columns').drop(['MolWt'])
    mol_prop_anhyd.columns = [
        str(col) + '_mol_prop' for col in mol_prop_anhyd.columns]
    return mol_prop_anhyd




# For oxide to wt% function
oxide_mass_all = {'SiO2_Liq': 60.084, 'MgO_Liq': 40.304, 'FeO_Liq': 71.846, 'Fe2O3_Liq':159.69,
    'CaO_Liq': 56.079, 'Al2O3_Liq': 101.961, 'Na2O_Liq': 61.979, 'K2O_Liq': 94.195, 'MnO_Liq': 70.937, 'TiO2_Liq': 79.898,
                          'P2O5_Liq':141.944, 'H2O_Liq': 18.01528}
oxide_mass_all = pd.DataFrame.from_dict(
    oxide_mass_all, orient='index').T
oxide_mass_all['Sample_ID'] = 'MolWt'
oxide_mass_all.set_index('Sample_ID', inplace=True)

elemental_mass_mult_all = {'SiO2_Liq': 28.0855, 'MgO_Liq': 24.305, 'FeO_Liq': 55.845,
'Fe2O3_Liq': 55.845*2, 'CaO_Liq': 40.078, 'Al2O3_Liq': 26.981539*2,
                          'Na2O_Liq': 22.989769*2, 'K2O_Liq': 39.0983*2, 'MnO_Liq': 54.938044, 'TiO2_Liq': 47.867,
                       'P2O5_Liq': 2*30.973762, 'H2O_Liq':1.00794*2}
elemental_mass_mult_all = pd.DataFrame.from_dict(
    elemental_mass_mult_all, orient='index').T
elemental_mass_mult_all['Sample_ID'] = 'ElWt'
elemental_mass_mult_all.set_index('Sample_ID', inplace=True)



df_ideal_all2 = pd.DataFrame(columns=['SiO2', 'TiO2', 'Al2O3',
'FeO', 'Fe2O3', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O',
'Cr2O3', 'P2O5', 'F', 'Cl', 'H2O'])

def convert_oxide_percent_to_element_weight_percent(df, Fe3Fet_Liq=None):


    df_c=df.copy()
    if Fe3Fet_Liq is not None:
        df_c['Fe3Fet_Liq']=Fe3Fet_Liq

    df_c['FeO_Liq']=df_c['FeOt_Liq']*(1-df_c['Fe3Fet_Liq'])
    df_c['Fe2O3_Liq']=df_c['FeOt_Liq']*(df_c['Fe3Fet_Liq'])*1.11111
    df_c.drop(columns=['FeOt_Liq'], inplace=True)

    df_oxides=df_c.reindex(oxide_mass_all.columns, axis=1).fillna(0)

    liq_wt_combo = pd.concat([oxide_mass_all, df_oxides],)


    mol_prop_anhyd = liq_wt_combo.div(
        liq_wt_combo.loc['MolWt', :], axis='columns').drop(['MolWt'])

    el_combo=pd.concat([elemental_mass_mult_all, mol_prop_anhyd ],)
    wt_perc = el_combo.multiply(
        el_combo.loc['ElWt', :], axis='columns').drop(['ElWt'])


    wt_perc2=pd.DataFrame(data={'Si_wt': wt_perc['SiO2_Liq'],
                                'Mg_wt': wt_perc['MgO_Liq'],
                                'Fe2_wt':wt_perc['FeO_Liq'],
                                'Fe3_wt':wt_perc['Fe2O3_Liq'],
                                'Ca_wt':wt_perc['CaO_Liq'],
                                'Al_wt':wt_perc['Al2O3_Liq'],
                                'Na_wt':wt_perc['Na2O_Liq'],
                                'K_wt':wt_perc['K2O_Liq'],
                                'Mn_wt':wt_perc['MnO_Liq'],
                                'Ti_wt':wt_perc['TiO2_Liq'],
                                'P_wt':wt_perc['P2O5_Liq'],
                                'H_wt': wt_perc['H2O_Liq'],



                                })
    sum_element=wt_perc2.sum(axis=1)

    mol_prop=calculate_mol_proportions_liquid(liq_comps=df_c)



    O_by_Charge=(15.9994/2*(mol_prop['SiO2_Liq_mol_prop']*4+mol_prop['TiO2_Liq_mol_prop']*4
    +2*mol_prop['Al2O3_Liq_mol_prop']*3+mol_prop['MgO_Liq_mol_prop']*2+mol_prop['MnO_Liq_mol_prop']*2+mol_prop['FeO_Liq_mol_prop']*2
    +2*mol_prop['Fe2O3_Liq_mol_prop']*3+mol_prop['CaO_Liq_mol_prop']*2+2*mol_prop['Na2O_Liq_mol_prop']*1+2*mol_prop['K2O_Liq_mol_prop']*1
    +2*mol_prop['P2O5_Liq_mol_prop']*5+df['H2O_Liq']/9.00764))

    wt_perc2['O_wt']=O_by_Charge

    return wt_perc2

atomic_wt=pd.DataFrame(data={'Sample_ID_Liq': "Atom_wt",
                                'Si_wt': 28.0855,
                                'Ti_wt': 47.867,
                                'Al_wt': 26.981538,
                                'Mg_wt':24.305,
                                'Mn_wt':54.938049,
                                'Fe2_wt':55.8457,
                                'Fe3_wt': 55.8457,
                                'Ca_wt': 40.078,
                                'Na_wt': 22.98977,
                                'K_wt': 39.0983,
                                'P_wt': 30.973761,
                                'H_wt': 1.00794,
                                'O_wt': 15.9994
                                    }, index=[0])

atomic_wt.set_index('Sample_ID_Liq', inplace=True)

def calculate_Smythe_silicate_mole_fractions(df, Fe3Fet_Liq=None):
    '''Smythe et al. call the results of these calculations "Silicate Mole Fractions"


   Parameters
    -------

    liq_comps: pandas.DataFrame
        Panda DataFrame of liquid compositions with column headings SiO2_Liq, TiO2_Liq etc.

    Returns
    -------
    pandas DataFrame
        anhydrous mole proportions for the liquid with column headings of the form SiO2_Liq_mol_prop

    '''
    if Fe3Fet_Liq is not None:
        el_wt=convert_oxide_percent_to_element_weight_percent(df, Fe3Fet_Liq)
    else:
        el_wt=convert_oxide_percent_to_element_weight_percent(df)


    # This makes the input match the columns in the oxide mass dataframe
    liq_wt = el_wt.reindex(
        atomic_wt.columns, axis=1).fillna(0)
    # Combine the molecular weight and weight percent dataframes
    liq_wt_combo = pd.concat([atomic_wt, liq_wt],)
    # Drop the calculation column
    atom_prop = liq_wt_combo.div(
        liq_wt_combo.loc['Atom_wt', :], axis='columns').drop(['Atom_wt'])
    atom_prop.columns = [
        str(col) + '_atom' for col in atom_prop.columns]

    atom_prop['sum'] = atom_prop.sum(axis='columns')
    atom_prop2 = atom_prop.div(atom_prop['sum'], axis='rows')
    atom_prop2.drop(['sum'], axis='columns', inplace=True)
    # atom_prop2.columns = [str(col).replace('Atom_wt', 'Atom_perc')
    #                           for col in Atom_wt.columns]
    a=atom_prop2*100

    a_drop=a.drop(columns=['O_wt_atom'])
    a_drop['sum'] = a_drop.sum(axis='columns')
    Final_Calc_step = a_drop.div(a_drop['sum'], axis='rows')
    Final_Calc_step.drop(['sum'], axis='columns', inplace=True)


    return Final_Calc_step


## Redox calculations
def convert_fo2_to_fe_partition(*, liq_comps, T_K, P_kbar,
                                model="Kress1991", fo2, renorm=False, fo2_offset=0):
    '''
    Calculates Fe3Fet_Liq, FeO and Fe2O3 based on user-specified buffer

   Parameters
    -------

    liq_comps: pandas.DataFrame
        Liquid compositions with column headings SiO2_Liq, MgO_Liq etc.

    T_K:  int, flt, pandas.Series
        Temperature in Kelvin (buffer positions are very T-sensitive)

    P_kbar: int, flt, pandas.Series
        Pressure in Kbar (Buffer positions are slightly sensitive to pressure)

    fo2:  str ("QFM", "NNO") or int, flt, pandas.Series
        Either a value of fo2 (enter 10*logfo2), or buffer position.
        So far, includes QFM or NNO

    fo2_offset: int, flt, pandas.Series
        log units offset from buffer, e.g., could specify fo2=QFM, fo2_offset=1
        to perform calculations at QFM+1

    model: str
        "Kress1991" - Uses Kress and Carmichael 1991 to calculate XFe2Fe3 from fo2
        "Put2016_eq6b" - Uses Putirka (2016) expression to calculate XFe2Fe3 from fo2

    renorm: bool
        Following excel code of K. Iacovino.
        If True, renormalizes other oxide concentrations
        to account for change in total following partitioning of Fe into FeO and Fe2O3.

    Returns
    -------

    liquid compositions with calculated Fe3Fet_Liq, FeO_Liq, Fe2O3_Liq, and XFe3Fe2.

    '''
    if isinstance(fo2, str):
        if fo2=="NNO":
        # Buffer position from frost (1991)
            logfo2=(-24930/T_K) + 9.36 + 0.046 * ((P_kbar*1000)-1)/T_K+fo2_offset
            fo2=10**logfo2

        if fo2=="QFM":
        # Buffer position from frost (1991)
            logfo2=(-25096.3/T_K) + 8.735 + 0.11 * ((P_kbar*1000)-1)/T_K+fo2_offset
            fo2=10**logfo2



    liq_comps_c=liq_comps.copy()
    mol_frac_short=calculate_hydrous_mol_fractions_liquid(liq_comps_c)
    mol_frac=pd.concat([mol_frac_short, liq_comps_c], axis=1)
    To=1673.15

    if model=="Kress1991":
        ln_XFe2FeO3_XFeO=((0.196*np.log(fo2))+(11492/T_K)-6.675+((-2.243*mol_frac['Al2O3_Liq_mol_frac'])+(-1.828*mol_frac['FeOt_Liq_mol_frac'])
        +(3.201*mol_frac['CaO_Liq_mol_frac'])+(5.854*mol_frac['Na2O_Liq_mol_frac'])+(6.215*mol_frac['K2O_Liq_mol_frac']))
        -3.36*(1-(To/T_K) - np.log(T_K/To)) -0.000000701*((P_kbar*100000000)/T_K)
         + -0.000000000154*(((T_K-1673)*(P_kbar*100000000))/T_K) + 0.0000000000000000385*((P_kbar*100000000)**2/T_K))
        #print(ln_XFe2FeO3_XFeO)
    if model=="Put2016_eq6b":
        ln_XFe2FeO3_XFeO=(-6.35+10813.8/T_K + 0.19*np.log(fo2)+ 12.4*(mol_frac['Na2O_Liq_mol_frac']
         +mol_frac['K2O_Liq_mol_frac'])
        -3.44*(mol_frac['Al2O3_Liq_mol_frac']/(mol_frac['Al2O3_Liq_mol_frac']+mol_frac['SiO2_Liq_mol_frac']))
        +4.15*mol_frac['CaO_Liq_mol_frac'])

    X_Fe2O3_X_FeO=np.exp(ln_XFe2FeO3_XFeO)
    X_Fe2O3=X_Fe2O3_X_FeO*mol_frac['FeOt_Liq_mol_frac']/(2*X_Fe2O3_X_FeO+1)

    #X_FeO=mol_frac['FeOt_Liq_mol_frac']/(2*X_Fe2O3_X_FeO+1) Kayla's way
    X_FeO=mol_frac['FeOt_Liq_mol_frac']-2*X_Fe2O3
    Sum_all_mol_frac=(mol_frac['SiO2_Liq_mol_frac']+mol_frac['TiO2_Liq_mol_frac']+mol_frac['Al2O3_Liq_mol_frac']+mol_frac['MnO_Liq_mol_frac']
                      +mol_frac['MgO_Liq_mol_frac']+mol_frac['CaO_Liq_mol_frac']+mol_frac['Na2O_Liq_mol_frac']+mol_frac['K2O_Liq_mol_frac']
                      +mol_frac['P2O5_Liq_mol_frac']+X_FeO+X_Fe2O3)

    Fe2O3_unnorm=X_Fe2O3*159.6
    FeO_unnorm=X_FeO*71.844
    Sum_All_mol=(mol_frac['SiO2_Liq_mol_frac']*60.0843+mol_frac['TiO2_Liq_mol_frac']*79.8788
    +mol_frac['Al2O3_Liq_mol_frac']*101.961+mol_frac['MnO_Liq_mol_frac']*70.9375
    +mol_frac['MgO_Liq_mol_frac']*40.3044+mol_frac['CaO_Liq_mol_frac']*56.0774+mol_frac['Na2O_Liq_mol_frac']*61.9789+mol_frac['K2O_Liq_mol_frac']*94.196
    +mol_frac['P2O5_Liq_mol_frac']*141.937+X_Fe2O3*159.6+X_FeO*71.844)
    New_Fe2O3_wt=(100*X_Fe2O3*159.6)/Sum_All_mol
    New_FeO_wt=(100*X_FeO*71.844)/Sum_All_mol

    New_Oxide_out_nonorm=liq_comps.copy()
    New_Oxide_out_nonorm['FeO_Liq']=New_FeO_wt
    New_Oxide_out_nonorm['Fe2O3_Liq']=New_Fe2O3_wt
    New_Oxide_out_nonorm['XFe3Fe2']=X_Fe2O3_X_FeO
    New_Oxide_out_nonorm['Fe3Fet_Liq']=New_Fe2O3_wt*0.8998/(New_FeO_wt+New_Fe2O3_wt*0.8998)


    New_Oxide_out_norm=pd.DataFrame(data={'SiO2_Liq': 100*mol_frac['SiO2_Liq_mol_frac']*60.084/Sum_All_mol,
                                         'TiO2_Liq': 100*mol_frac['TiO2_Liq_mol_frac']*79.8788/Sum_All_mol,
                                         'Al2O3_Liq':100*mol_frac['Al2O3_Liq_mol_frac']*101.961/Sum_All_mol,
                                          'Fe2O3_Liq': (100*X_Fe2O3*159.6)/Sum_All_mol,
                                          'FeO_Liq': (100*X_FeO*71.844)/Sum_All_mol,
                                          'MnO_Liq': 100*mol_frac['MnO_Liq_mol_frac']*70.9375/Sum_All_mol,
                                          'MgO_Liq': 100*mol_frac['MgO_Liq_mol_frac']*40.3044/Sum_All_mol,
                                         'CaO_Liq': 100*mol_frac['CaO_Liq_mol_frac']*56.0774/Sum_All_mol,
                                          'Na2O_Liq': 100*mol_frac['Na2O_Liq_mol_frac']*61.9789/Sum_All_mol,
                                          'K2O_Liq': 100*mol_frac['K2O_Liq_mol_frac']*94.196/Sum_All_mol,
                                         'P2O5_Liq':  100*mol_frac['P2O5_Liq_mol_frac']*141.937/Sum_All_mol,
                                         })
    Old_Sum=(100/liq_comps_c.drop(['Sample_ID_Liq'], axis=1).sum(axis=1))
    New_Oxide_out_New_old_total=New_Oxide_out_norm.div(Old_Sum, axis=0)
    New_Oxide_out_New_old_total['Fe3Fet_Liq']=(New_Oxide_out_norm['Fe2O3_Liq']*0.8998/(New_Oxide_out_norm['FeO_Liq']+New_Oxide_out_norm['Fe2O3_Liq']*0.8998)).fillna(0)



    if renorm==False:
        return New_Oxide_out_nonorm
    else:
        return New_Oxide_out_New_old_total


## Need some functions for calculating mole proportions with Fe partition
oxide_mass_liq_hyd_redox = {'SiO2_Liq': 60.0843, 'MgO_Liq': 40.3044,
'MnO_Liq': 70.9375, 'FeO_Liq': 71.844, 'Fe2O3_Liq': 159.69, 'CaO_Liq': 56.0774,
'Al2O3_Liq': 101.961,'Na2O_Liq': 61.9789, 'K2O_Liq': 94.196,
 'TiO2_Liq': 79.8788, 'P2O5_Liq': 141.937, 'Cr2O3_Liq': 151.9982,
  'H2O_Liq': 18.01528}
# Turns dictionary into a dataframe so pandas matrix math functions can be used
oxide_mass_liq_hyd_df_redox = pd.DataFrame.from_dict(
    oxide_mass_liq_hyd_redox, orient='index').T
oxide_mass_liq_hyd_df_redox['Sample_ID_Liq'] = 'MolWt'
oxide_mass_liq_hyd_df_redox.set_index('Sample_ID_Liq', inplace=True)


cation_num_liq_hyd_redox = {'SiO2_Liq': 1, 'MgO_Liq': 1, 'MnO_Liq': 1,
'FeO_Liq': 1, 'Fe2O3_Liq': 2, 'CaO_Liq': 1, 'Al2O3_Liq': 2, 'Na2O_Liq': 2,
'K2O_Liq': 2, 'TiO2_Liq': 1, 'P2O5_Liq': 2, 'H2O_Liq':2,
'Cr2O3_Liq':2}

cation_num_liq_hyd_df_redox = pd.DataFrame.from_dict(
    cation_num_liq_hyd_redox, orient='index').T
cation_num_liq_hyd_df_redox['Sample_ID_Liq'] = 'CatNum'
cation_num_liq_hyd_df_redox.set_index('Sample_ID_Liq', inplace=True)



def calculate_hydrous_mol_proportions_liquid_redox(liq_comps):
    '''Import Liq compositions using liq_comps=My_Liquids, returns anhydrous mole proportions

   Parameters
    -------


    liq_comps: pandas.DataFrame
        liquid compositions with column headings SiO2_Liq, TiO2_Liq etc.

    Returns
    -------
    pandas DataFrame
        anhydrous mole proportions for the liquid with column headings of the ..Liq_mol_prop

    '''
    # This makes the input match the columns in the oxide mass dataframe
    liq_wt = liq_comps.reindex(oxide_mass_liq_hyd_df_redox.columns, axis=1).fillna(0)
    # Combine the molecular weight and weight percent dataframes
    liq_wt_combo = pd.concat([oxide_mass_liq_hyd_df_redox, liq_wt],)
    # Drop the calculation column
    mol_prop_hyd = liq_wt_combo.div(
        liq_wt_combo.loc['MolWt', :], axis='columns').drop(['MolWt'])
    mol_prop_hyd.columns = [
        str(col) + '_mol_prop' for col in mol_prop_hyd.columns]
    return mol_prop_hyd

def calculate_hydrous_mol_fractions_liquid_redox(liq_comps):
    '''Import Liq compositions using liq_comps=My_Liquids, returns anhydrous mole fractions

   Parameters
    -------

    liq_comps: pandas.DataFrame
        liquid compositions with column headings SiO2_Liq, TiO2_Liq etc.



    Returns
    -------
    pandas DataFrame
        anhydrous mole fractions for the liquid with column headings of the form SiO2_Liq_mol_frac

    '''
    mol_prop = calculate_hydrous_mol_proportions_liquid_redox(liq_comps)
    mol_prop['sum'] = mol_prop.sum(axis='columns')
    mol_frac_hyd = mol_prop.div(mol_prop['sum'], axis='rows')
    mol_frac_hyd.drop(['sum'], axis='columns', inplace=True)
    mol_frac_hyd.columns = [str(col).replace('prop', 'frac')
                            for col in mol_frac_hyd.columns]
    return mol_frac_hyd

## Fractions with redox

oxide_mass_liq_hyd_redox = {'SiO2_Liq': 60.0843, 'MgO_Liq': 40.3044,
'MnO_Liq': 70.9375, 'FeO_Liq': 71.8464, 'Fe2O3_Liq': 159.69, 'CaO_Liq': 56.0774,
'Al2O3_Liq': 101.961, 'Na2O_Liq': 61.9789, 'K2O_Liq': 94.196,
'TiO2_Liq': 79.8788, 'P2O5_Liq': 141.944, 'H2O_Liq': 18.02 }

oxide_mass_liq_hyd_redox_df = pd.DataFrame.from_dict(
    oxide_mass_liq_hyd_redox, orient='index').T
oxide_mass_liq_hyd_redox_df['Sample_ID_Liq'] = 'MolWt'
oxide_mass_liq_hyd_redox_df.set_index('Sample_ID_Liq', inplace=True)


def calculate_hydrous_cat_proportions_liquid_redox(liq_comps):
    '''Import Liq compositions using liq_comps=My_Liquids, returns hydrous cation proportions (e.g., mole proportions * no of cations)

   Parameters
    -------

    liq_comps: pandas.DataFrame
                Panda DataFrame of liquid compositions with column headings SiO2_Liq, TiO2_Liq etc.


    Returns
    -------
    pandas DataFrame
        hydrous cation proportions for the liquid with column headings of the form S_Liq_cat_prop

    '''

    mol_prop_no_cat_num = calculate_hydrous_mol_proportions_liquid_redox(liq_comps)
    mol_prop_no_cat_num.columns = [str(col).replace(
        '_mol_prop', '') for col in mol_prop_no_cat_num.columns]
    ox_num_reindex = cation_num_liq_hyd_df_redox.reindex(
        oxide_mass_liq_hyd_redox_df.columns, axis=1).fillna(0)
    df_calc_comb = pd.concat([ox_num_reindex, mol_prop_no_cat_num])
    cation_prop_hyd = df_calc_comb.multiply(
        df_calc_comb.loc['CatNum', :], axis='columns').drop(['CatNum'])
    cation_prop_hyd.columns = [
        str(col) + '_cat_prop' for col in cation_prop_hyd.columns]



    return cation_prop_hyd

def calculate_hydrous_cat_fractions_liquid_redox(liq_comps):
    '''Import Liq compositions using liq_comps=My_Liquids, returns hydrous cation fractions

   Parameters
    -------

    liq_comps: pandas.DataFrame
        liquid compositions with column headings SiO2_Liq, TiO2_Liq etc.


    Returns
    -------
    pandas DataFrame
        hydrous cation fractions for the liquid with column headings of the form _Liq_cat_frac,
        as well as the initial dataframe of liquid compositions.


    '''
    cat_prop = calculate_hydrous_cat_proportions_liquid_redox(liq_comps=liq_comps)
    mol_prop = calculate_hydrous_mol_fractions_liquid_redox(liq_comps=liq_comps)

    cat_prop['sum'] = cat_prop.sum(axis='columns')
    cat_frac_hyd = cat_prop.div(cat_prop['sum'], axis='rows')
    cat_frac_hyd.drop(['sum'], axis='columns', inplace=True)
    cat_frac_hyd.columns = [str(col).replace('prop', 'frac')
                              for col in cat_frac_hyd.columns]
    cat_frac_hyd = pd.concat([liq_comps, mol_prop, cat_frac_hyd], axis=1)



    cation_frac_hyd2=cat_frac_hyd.rename(columns={
                        'SiO2_Liq_cat_frac': 'Si_Liq_cat_frac',
                        'TiO2_Liq_cat_frac': 'Ti_Liq_cat_frac',
                        'Al2O3_Liq_cat_frac': 'Al_Liq_cat_frac',
                        'FeO_Liq_cat_frac': 'Fe2_Liq_cat_frac',
                        'Fe2O3_Liq_cat_frac': 'Fe3_Liq_cat_frac',
                        'MnO_Liq_cat_frac': 'Mn_Liq_cat_frac',
                        'MgO_Liq_cat_frac': 'Mg_Liq_cat_frac',
                        'CaO_Liq_cat_frac': 'Ca_Liq_cat_frac',
                        'Na2O_Liq_cat_frac': 'Na_Liq_cat_frac',
                        'K2O_Liq_cat_frac': 'K_Liq_cat_frac',
                        'Cr2O3_Liq_cat_frac': 'Cr_Liq_cat_frac',
                        'P2O5_Liq_cat_frac': 'P_Liq_cat_frac',

                        })

    return cation_frac_hyd2



def convert_fe_partition_to_fo2(*, liq_comps, T_K, P_kbar,  model="Kress1991", renorm=False):
    '''
    Calculates delta fo2 relative to QFM and NNO buffer for liq compositions with FeO and Fe2O3

   Parameters
    -------

    liq_comps: pandas.DataFrame
        Liquid compositions with column headings SiO2_Liq, MgO_Liq, FeO_Liq and Fe2O3_Liq etc.

    T_K:  int, flt, pandas.Series
        Temperature in Kelvin (buffer positions are very T-sensitive)

    P_kbar: int, flt, pandas.Series
        Pressure in Kbar (Buffer positions are slightly sensitive to pressure)






    model: str
        "Kress1991" - Uses Kress and Carmichael 1991 to calculate XFe2Fe3 from fo2
        "Put2016_eq6b" - Uses Putirka (2016) expression to calculate XFe2Fe3 from fo2

    renorm: bool
        Following excel code of K. Iacovino.
        If True, renormalizes other oxide concentrations
        to account for change in total following partitioning of Fe into FeO and Fe2O3.

    Returns
    -------

    liquid compositions with calculated Fe3Fet_Liq, FeO_Liq, Fe2O3_Liq, and XFe3Fe2.

    '''

    liq_comps_c=liq_comps.copy()
    hyd_mol_frac_test=calculate_hydrous_mol_fractions_liquid(liq_comps=liq_comps_c)

    liq_comps_c['FeO_Liq']=liq_comps_c['FeOt_Liq']*(1-liq_comps_c['Fe3Fet_Liq'])
    liq_comps_c['Fe2O3_Liq']=liq_comps_c['FeOt_Liq']*(liq_comps_c['Fe3Fet_Liq'])*1.11111

    mol_frac_hyd_redox=calculate_hydrous_mol_fractions_liquid_redox(liq_comps=liq_comps_c)

    # Calculating buffer positions from Frost 1991

    To= 1673.15

    logfo2_NNO=(-24930/T_K) + 9.36 + 0.046 * ((P_kbar*1000)-1)/T_K
    fo2_NNO=10**logfo2_NNO



    logfo2_QFM=(-25096.3/T_K) + 8.735 + 0.11 * ((P_kbar*1000)-1)/T_K
    fo2_QFM=10**logfo2_QFM

    # This is Ln (XFe2O3/XFeO) from the Kress and Carmichael 1991 paper
    Z=np.log(mol_frac_hyd_redox['Fe2O3_Liq_mol_frac']/
         (mol_frac_hyd_redox['FeO_Liq_mol_frac']))

    # We've simplified the equatoin down to Z= a ln fo2 + rightside

    rightside=( (11492/T_K)-6.675+((-2.243*mol_frac_hyd_redox['Al2O3_Liq_mol_frac'])+(-1.828*hyd_mol_frac_test['FeOt_Liq_mol_frac'])
    +(3.201*mol_frac_hyd_redox['CaO_Liq_mol_frac'])+(5.854*mol_frac_hyd_redox['Na2O_Liq_mol_frac'])+(6.215*mol_frac_hyd_redox['K2O_Liq_mol_frac']))
    -3.36*(1-(To/T_K) - np.log(T_K/To)) -0.000000701*((P_kbar*100000000)/T_K)
    + -0.000000000154*(((T_K-1673)*(P_kbar*100000000))/T_K) + 0.0000000000000000385*((P_kbar*100000000)**2/T_K)
    )

    ln_fo2_calc=(Z-rightside)/0.196
    fo2_calc=np.exp(ln_fo2_calc)
    # and back to log base 10
    log_fo2_calc=np.log10(fo2_calc)
    DeltaQFM=log_fo2_calc-logfo2_QFM
    DeltaNNO=log_fo2_calc-logfo2_NNO


    liq_comps_c.insert(0, 'DeltaQFM', DeltaQFM)
    liq_comps_c.insert(1, 'DeltaNNO', DeltaNNO)
    liq_comps_c.insert(2, 'fo2_calc', fo2_calc)
    return liq_comps_c


    ## Converting between S, SO2, SO3 etc.
mol_mass_S=32.065
mol_mass_O=15.999
mol_mass_SO3=mol_mass_S+mol_mass_O*3
mol_mass_SO4=mol_mass_S+mol_mass_O*4
mol_mass_SO2=mol_mass_S+mol_mass_O*2

def convert_S_types(SO3_wt=None, SO3_ppm=None, S_wt=None, S_ppm=None, SO2_wt=None, SO2_ppm=None,
 SO4_wt=None, SO4_ppm=None):
    """ converts SO3 in wt% into S in ppm
    """
    params = {
        "a": SO3_wt,
        "ab": SO3_ppm,
        "b": S_wt,
        "c": S_ppm,
        "d": SO2_wt,
        "de": SO2_ppm,
        "e": SO4_wt,
        "f": SO4_ppm
    }

    not_none_params = {k:v for k, v in params.items() if v is not None}
    if len(not_none_params)>1:
        raise TypeError('Please only enter one input type, the function returns a dataframe of all the other outputs')

    else:
        if S_ppm is not None:
            moles_s=(S_ppm/10**4)/mol_mass_S
        if S_wt is not None:
            moles_s=S_wt/mol_mass_S
        if SO3_wt is not None:
            moles_s=SO3_wt/mol_mass_SO3
        if SO3_ppm is not None:
            moles_s=(SO3_ppm/10**4)/mol_mass_SO3

        if SO4_wt is not None:
            moles_s=SO4_wt/mol_mass_SO4
        if SO4_ppm is not None:
            moles_s=(SO4_ppm/10**4)/mol_mass_SO4

        if SO2_wt is not None:
            moles_s=SO2_wt/mol_mass_SO2
        if SO2_ppm is not None:
            moles_s=(SO2_ppm/10**4)/mol_mass_SO2

        S_wt2= moles_s*mol_mass_S
        S_ppm2=S_wt2*10**4
        SO2=moles_s*mol_mass_SO2
        SO3=moles_s*mol_mass_SO3
        SO4=moles_s*mol_mass_SO4


        df=pd.DataFrame(data={
                                'S_wt': S_wt2,
                                'S_ppm': S_ppm2,
                                'SO2_wt': SO2,
                                'SO2_ppm': SO2*10**4,
                                'SO3_wt': SO3,
                                'SO3_ppm': SO3*10**4,
                                'SO4_wt': SO4,
                                'SO4_ppm': SO4*10**4})



        return df

