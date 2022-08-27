import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.optimize as optimize
from scipy.special import erf

## Import data

df_ideal_liq = pd.DataFrame(columns=['SiO2_Liq', 'TiO2_Liq', 'Al2O3_Liq',
'FeOt_Liq', 'MnO_Liq', 'MgO_Liq', 'CaO_Liq', 'Na2O_Liq', 'K2O_Liq',
'Cr2O3_Liq', 'P2O5_Liq', 'H2O_Liq', 'Fe3Fet_Liq', 'Ni_Liq_ppm', 'Cu_Liq_ppm'])

def import_data(filename, sheet_name, sample_label=None, suffix=None):
    if 'csv' in filename:
        my_input = pd.read_csv(filename)

    elif 'xls' in filename:
        if sheet_name is not None:
            my_input = pd.read_excel(filename, sheet_name=sheet_name)
            #my_input[my_input < 0] = 0
        else:
            my_input = pd.read_excel(filename)
            #my_input[my_input < 0] = 0

    if suffix is not None:
        if any(my_input.columns.str.contains(suffix)):
            w.warn('We notice you have specified a suffix, but some of your columns already have this suffix. '
        'e.g., If you already have _Liq in the file, you shouldnt specify suffix="_Liq" during the import')

    my_input_c = my_input.copy()
    if suffix is not None:
        my_input_c=my_input_c.add_suffix(suffix)

    if any(my_input.columns.str.contains("FeO_")) and (all(my_input.columns.str.contains("FeOt_")==False)):
        raise ValueError("No FeOt found. You've got a column heading with FeO. To avoid errors based on common EPMA outputs"
    " thermobar only recognises columns with FeOt for all phases except liquid"
    " where you can also enter a Fe3Fet_Liq heading used for equilibrium tests")

    if any(my_input.columns.str.contains("Fe2O3_")) and (all(my_input.columns.str.contains("FeOt_")==False)):
        raise ValueError("No FeOt column found. You've got a column heading with Fe2O3. To avoid errors based on common EPMA outputs"
        " thermobar only recognises columns with FeOt for all phases except liquid"
        " where you can also enter a Fe3Fet_Liq heading used for equilibrium tests")

    if any(my_input.columns.str.contains("FeOT_")) and (all(my_input.columns.str.contains("FeOt_")==False)):
        raise ValueError("No FeOt column found. You've got a column heading with FeOT. Change to a lower case t")



    myLiquids1 = my_input_c.reindex(df_ideal_liq.columns, axis=1).fillna(0)
    myLiquids1 = myLiquids1.apply(pd.to_numeric, errors='coerce').fillna(0)
    myLiquids1[myLiquids1 < 0] = 0
    print('We have replaced all missing liquid oxides and strings with zeros. ')

    cols2=myLiquids1.columns
    my_input_c=my_input.copy()
    for col in cols2:
        if col in my_input_c.columns:
            my_input_c=my_input_c.drop(columns=col)

    out=pd.concat([myLiquids1, my_input_c], axis=1)
    return out






## Anhydrous Molar masses, used for ONeill (2021). Doesnt include P2O5, Cr2O3 or H2O in the calculation.
oxide_mass_liq_anhyd = {'SiO2_Liq': 60.0843, 'MgO_Liq': 40.3044,
'MnO_Liq': 70.9375, 'FeOt_Liq': 71.8464, 'CaO_Liq': 56.0774,
'Al2O3_Liq': 101.961, 'Na2O_Liq': 61.9789, 'K2O_Liq': 94.196,
'TiO2_Liq': 79.8788}


cation_num_liq_anhyd = {'SiO2_Liq': 1, 'MgO_Liq': 1, 'MnO_Liq': 1,
'FeOt_Liq': 1, 'CaO_Liq': 1, 'Al2O3_Liq': 2, 'Na2O_Liq': 2,
'K2O_Liq': 2, 'TiO2_Liq': 1}


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
'TiO2_Liq': 79.8788, 'P2O5_Liq': 141.944, 'H2O_Liq': 18}


cation_num_liq_hyd = {'SiO2_Liq': 1, 'MgO_Liq': 1, 'MnO_Liq': 1,
'FeOt_Liq': 1, 'CaO_Liq': 1, 'Al2O3_Liq': 2, 'Na2O_Liq': 2,
'K2O_Liq': 2, 'TiO2_Liq': 1, 'P2O5_Liq': 2, 'H2O_Liq':2}


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





## Anhydrous calculations from Thermobar (Wieser et al. in prep)


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
## Hydrous calculations for Fortin et al. (2015)

## hydrous calculations from Thermobar (Wieser et al. in prep)


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
    mol_frac_hyd = mol_prop.div(mol_prop['sum'], axis='rows')
    mol_frac_hyd.drop(['sum'], axis='columns', inplace=True)
    mol_frac_hyd.columns = [str(col).replace('prop', 'frac')
                              for col in mol_frac_hyd.columns]
    return mol_frac_hyd


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
    el_wt=convert_oxide_percent_to_element_weight_percent(df, Fe3Fet_Liq=None)

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



## Stuff for Oneill Expression


def calculate_Oneill2021_SCSS(*, df, T_K, P_kbar, Fe3Fet_Liq=None, Fe_FeNiCu_Sulf=None,
Cu_FeNiCu_Sulf=None, Ni_FeNiCu_Sulf=None, Fe_Sulf=None, Ni_Sulf=None, Cu_Sulf=None,
Ni_Liq=None, Cu_Liq=None):

    '''
    Calculates SCSS using the model of O'Neill (2021), with options for users to
    calculate sulfide composition from liquid composition using the approach of Smythe et al. (2017),
    the empirical parameterization of O'Neill (2021),  or input a sulfide composition.

    Parameters
    -------
    df: pandas.DataFrame
        Dataframe of liquid compositions. Needs to have the headings "SiO2_Liq", "TiO2_Liq" etc, with
        compositions in oxide wt%. all FeO should be entered as FeOt_Liq, which can then be partitioned
        using the Fe3Fet_Liq input. Heading order doesn't matter.

    T_K: int, float, pandas.Series
        Temperature in Kelvin.

    P_kbar: int, float, pandas.Series
        Pressure in kbar

    Fe3Fet_Liq: int, float, pandas.Series, or "Calc_ONeill"
        Proportion of Fe3+ in the liquid, as various parts of the calculation use only Fe2.
        If "Calc_ONeill" Calculates as a function of MgO content, using the MORB relationship.

    Sulfide Composition: Options to calculate from the liquid composition, enter the comp in el wt%,
    or enter the FeFeNiCu, Cu

      if you want to calculate sulfide composition:

            Fe_FeNiCu_Sulf = "Calc_Smythe"
            Calculates sulfide composition analagous to the Solver function in the Smythe et al. (2017) spreadsheet.
            Here, we use Scipy optimize to find the ideal Ni and Cu
            contents using Kds from Kiseeva et al. (2015) and the Ni and Cu content of the melt.
            Requires user to also enter Ni_Liq and Cu_Liq.

            Or

            Fe_FeNiCu_Sulf = "Calc_ONeill"
            Calculates sulfide composition using the empirical expression of O'Neill (2021), which depends on
            FeOt_Liq, Ni_Liq, Cu_Liq, and Fe3Fet_Liq. We allow users to enter their own Fe3Fet_Liq,
            as we believe the empirical model of Neill where Fe3Fet_Liq is a function of MgO content is not
            broadly applicable.

        if you want to input a Fe_FeNiCu_Sulf ratio:
            Fe_FeNiCu_Sulf = int, float, pandas series
            Calculates SCSS using this ratio.
            If you want the non-ideal SCSS to be returned, you also need to enter
            values for Cu_FeNiCu_Sulf and Ni_FeNiCu_Sulf

        if you want to input a measured sulfide composition in el wt%
            Fe_Sulf, Ni_Sulf, Cu_Sulf = int, float, pandas series




    Returns
    -------
    pandas.DataFrame:
        Contains column for SCSS ideal, 1 sigma, input T and P, and the various intermediate steps of the calculation.

    '''
    df_c=df.copy()
    # if 'P2O5_Liq' in df_c.columns:
    #     df_c['P2O5_Liq']=0
    liqs=calculate_anhydrous_cat_fractions_liquid(liq_comps=df_c)

    if isinstance(Fe3Fet_Liq, str) and Fe3Fet_Liq == "Calc_ONeill":
        Fe2O3_Calc=np.exp(1.46-0.177*df_c['MgO_Liq'])
        Fe3Fet_Liq=Fe2O3_Calc*0.8998/(df_c['FeOt_Liq'])

    liqs['Fe2_Liq_cat_frac']=liqs['Fet_Liq_cat_frac']*(1-Fe3Fet_Liq)

    # Then check options, basically need to have entered FeFeNiCu_Calc, or Fe_FeNiCu ratio, or all sulfide
    if Fe_Sulf is not None and Ni_Sulf is not None and Cu_Sulf is not None and Fe_FeNiCu_Sulf is not None:
        raise ValueError('You have entered both a Fe_FeNiCu_Sulf ratio, and the conc of Fe, Ni and Cu in your sulf. Please enter one set of inputs or another')

    elif isinstance(Fe_FeNiCu_Sulf, str) and (Fe_FeNiCu_Sulf=="Calc_Smythe" or Fe_FeNiCu_Sulf == "Calc_ONeill"):
        if Ni_Liq is None:
            raise ValueError('If you select Fe_FeNiCu_Sulf = Calc, you need to enter the concentration of Cu and Ni in the liquid in ppm"')
        if Cu_Liq is None:
            raise ValueError('If you select Fe_FeNiCu_Sulf = Calc, you need to enter the concentration of Cu and Ni in the liquid in ppm"')
        print('Calculating Sulf composition using Scipy minimisation of Kiseeva et al. (2015) Kd Values for Ni and Cu')

    elif Fe_FeNiCu_Sulf is not None:
        print('Using inputted Fe_FeNiCu_Sulf ratio for calculations.')

    elif Fe_Sulf is not None and Ni_Sulf is not None and Cu_Sulf is not None:
        print('Using inputted Sulf compositions to calculate Fe_FeNiCu_Sulf ratios.')

    else:
        raise ValueError('Input for sulfide composition not recognised.')


    if isinstance(Fe_FeNiCu_Sulf, str):
        if Fe_FeNiCu_Sulf=="Calc_Smythe":


            # This does the Scipy minimisation of Cu and Ni contents using Kiseeva et al. (2015)
            calc_sulf=calculate_Symthe_Sulf_Composition(FeOt_Liq=df_c['FeOt_Liq'], Fe3Fet_Liq=df_c['Fe3Fet_Liq'],
                                            T_K=T_K, Ni_Liq=Ni_Liq, Cu_Liq=Cu_Liq)

            # This feeds those result back into a simpler function to get the Fe, S and O content of the sulfide
            Sulf_All=calculate_Kiseeva_sulf_comp_kd(Ni_Sulf=calc_sulf['Ni_Sulf'],Cu_Sulf=calc_sulf['Cu_Sulf'],
                                FeOt_Liq=df_c['FeOt_Liq'], Fe3Fet_Liq=Fe3Fet_Liq,
                                            T_K=T_K, Ni_Liq=df_c['Ni_Liq (ppm)'], Cu_Liq=df_c['Cu_Liq (ppm)'])

            Fe_FeNiCu_Sulf=calculate_sulf_FeFeNiCu(Sulf_All['Ni_Sulf'], Sulf_All['Cu_Sulf'], Sulf_All['Fe_Sulf'])
            Cu_FeNiCu_Sulf=calculate_sulf_CuFeNiCu(Sulf_All['Ni_Sulf'], Sulf_All['Cu_Sulf'], Sulf_All['Fe_Sulf'])
            Ni_FeNiCu_Sulf=calculate_sulf_NiFeNiCu(Sulf_All['Ni_Sulf'], Sulf_All['Cu_Sulf'], Sulf_All['Fe_Sulf'])

            Smythe_calcs['Ni_Sulf_Calc']=Sulf_All['Ni_Sulf']
            Smythe_calcs['Cu_Sulf_Calc']=Sulf_All['Cu_Sulf']
            Smythe_calcs['Fe_Sulf_Calc']=Sulf_All['Fe_Sulf']
            Smythe_calcs['O_Sulf_Calc']=Sulf_All['O_Sulf']
            Smythe_calcs['S_Sulf_Calc']=Sulf_All['S_Sulf']

        if Fe_FeNiCu_Sulf=="Calc_ONeill":
            Fe_FeNiCu_Sulf=calculate_ONeill_sulf(FeOt_Liq=df_c['FeOt_Liq'],
            Ni_Liq=Ni_Liq, Cu_Liq=Cu_Liq, Fe3Fet_Liq=Fe3Fet_Liq)

    if Fe_Sulf is not None and Ni_Sulf is not None and Cu_Sulf is not None:
            Fe_FeNiCu_Sulf=calculate_sulf_FeFeNiCu(Ni_Sulf, Cu_Sulf, Fe_Sulf)
            Cu_FeNiCu_Sulf=calculate_sulf_CuFeNiCu(Ni_Sulf, Cu_Sulf, Fe_Sulf)
            Ni_FeNiCu_Sulf=calculate_sulf_NiFeNiCu(Ni_Sulf, Cu_Sulf, Fe_Sulf)




    liqs['LnCS2_calc']=(8.77-23590/T_K+(1673/T_K)*(6.7*(liqs['Na_Liq_cat_frac']+liqs['K_Liq_cat_frac'])
    +4.9*liqs['Mg_Liq_cat_frac']+8.1*liqs['Ca_Liq_cat_frac']+8.9*(liqs['Fet_Liq_cat_frac']+liqs['Mn_Liq_cat_frac'])
    +5*liqs['Ti_Liq_cat_frac']+1.8*liqs['Al_Liq_cat_frac']
    -22.2*liqs['Ti_Liq_cat_frac']*(liqs['Fet_Liq_cat_frac']+liqs['Mn_Liq_cat_frac'])
        +7.2*(liqs['Fet_Liq_cat_frac']*liqs['Si_Liq_cat_frac']))-2.06*erf(-7.2*(liqs['Fet_Liq_cat_frac']+liqs['Mn_Liq_cat_frac'])))

    liqs['DeltaG']=((137778-91.666*T_K+8.474*T_K*np.log(T_K))/(8.31441*T_K)+(-291*(P_kbar/10)+351*erf((P_kbar/10)))/T_K)

    liqs['Ln_a_FeS']=(np.log(Fe_FeNiCu_Sulf*(1-liqs['Fe2_Liq_cat_frac'])))

    liqs['Ln_a_FeO']=( np.log(liqs['Fe2_Liq_cat_frac'])+(((1-liqs['Fe2_Liq_cat_frac'])**2)*
            (28870-14710*liqs['Mg_Liq_cat_frac']+1960*liqs['Ca_Liq_cat_frac']+43300*liqs['Na_Liq_cat_frac']+95380*liqs['K_Liq_cat_frac']-76880*liqs['Ti_Liq_cat_frac'])
            +(1-liqs['Fe2_Liq_cat_frac'])*(-62190*liqs['Si_Liq_cat_frac']+31520*liqs['Si_Liq_cat_frac']*liqs['Si_Liq_cat_frac']))/(8.31441*T_K))

    liqs['LnS']=(liqs['LnCS2_calc']+liqs['DeltaG']+liqs['Ln_a_FeS']-liqs['Ln_a_FeO'])

    liqs['SCSS2_ppm']=np.exp(liqs['LnS'])

    cols_to_move = ['SCSS2_ppm', 'LnS', "Ln_a_FeO",
                    'Ln_a_FeS', 'DeltaG', 'LnCS2_calc']

    Liqs = liqs[cols_to_move +
                                    [col for col in liqs.columns if col not in cols_to_move]]

    return Liqs




## Calculating sulfide composition using the Smythe Minimisation

def Loop_Smythe_Sulf_Calc_Residual(single_argx0, FeO_Liq, T_K,  Ni_Liq, Cu_Liq):
    #, x1, *, FeO_Liq, T_K,  Ni_Liq, Cu_Liq):
    Ni_Sulf=single_argx0[0]
    Cu_Sulf=single_argx0[1]

    OCalc_CellAG12=0.24*FeO_Liq
    FeS_calc_AL19=((100-Ni_Sulf-Cu_Sulf-OCalc_CellAG12-Cu_Sulf*20.1442646/79.8557354
    -Ni_Sulf*35.32650016/64.67349984))*36.47119049/100
    FeCalc_CellAF12=FeS_calc_AL19*63.52880951/36.47119049

    AG6=Ni_Sulf/58.6934/(Ni_Sulf/58.6934+Cu_Sulf/63.546+FeCalc_CellAF12/55.845)
    AG8=AG6*0.97
    AG9=AG6*0.92

        # If FeO<13
    if FeO_Liq<13:
        O_Sulf=(FeO_Liq*0.24*((1-AG8)**2))
    elif FeO_Liq>=13:
        O_Sulf=(FeO_Liq*0.24*((1-AG9)**2))

    AL17=(Ni_Sulf*35.32650016/64.67349984)
    AL18=(Cu_Sulf*20.1442646/79.8557354)
    AM20=O_Sulf*77.72983506/22.27016494

    S_Sulf=((100-AL17-AL18-(O_Sulf*77.72983506/22.27016494)
             -O_Sulf-Cu_Sulf-Ni_Sulf)*36.47119049/100 + AL17 + AL18)

    Fe_Sulf=((100-Ni_Sulf-Cu_Sulf-O_Sulf-AL17-AL18-AM20)
    *(63.52880951/100)+O_Sulf*77.72983506/22.27016494)

    NiSS=Ni_Sulf*(35.32650016/64.67349984)
    CuS05=Cu_Sulf*20.1442646/79.8557354
    OCalc=O_Sulf*77.72983506/22.27016494
    S_Sulf=(100-NiSS-CuS05-OCalc-O_Sulf-Cu_Sulf-Ni_Sulf)*36.47119049/100+NiSS+CuS05

    # Cu Coefficients
    A_Cu=0.6995
    B_Cu=4200
    eps_Cu_Cu=-1.0310
    eps_FeO_Cu=2.2030
    eps_Ni_Cu=0.00
    # Ni Coefficients

    A_Ni=1.869
    B_Ni=3300
    eps_Cu_Ni=0.00
    eps_Ni_Ni=-0.904
    eps_FeO_Ni=0




    Fe_Corr=FeO_Liq/(Fe_Sulf/55.845/(Ni_Sulf/58.6934+Cu_Sulf/63.546+Fe_Sulf/55.845)) #AF7
    Ni_FeNiCu=Ni_Sulf/58.6934/(Ni_Sulf/58.6934+Cu_Sulf/63.546+Fe_Sulf/55.845)
    Cu_FeNiCu=Cu_Sulf/63.546/(Ni_Sulf/58.6934+Cu_Sulf/63.546+Fe_Sulf/55.845)
    Fe_FeNiCu=Fe_Sulf/55.845/(Ni_Sulf/58.6934+Cu_Sulf/63.546+Fe_Sulf/55.845)

    DNi=(10**(B_Ni/T_K+A_Ni-np.log10(Fe_Corr)+(1673*eps_FeO_Ni*np.log10(1-0.049*O_Sulf)/T_K)
            +(1673*eps_Ni_Ni*np.log10(1-Ni_FeNiCu)/T_K)+(1673*eps_Cu_Ni*np.log10(1-Cu_FeNiCu)/T_K)))

    DCu=(10**(B_Cu/T_K+A_Cu-0.5*np.log10(Fe_Corr)+(1673*eps_FeO_Cu*np.log10(1-0.049*O_Sulf)/T_K)
             +(1673*eps_Cu_Cu*np.log10(1-Cu_FeNiCu)/T_K)+(1673*eps_Ni_Cu*np.log10(1-Ni_FeNiCu)/T_K)))

    Ni_Melt_calc=Ni_Sulf*10000/DNi
    Cu_Melt_calc=Cu_Sulf*10000/DCu

    Residual=(Ni_Melt_calc-Ni_Liq)**2 + (Cu_Melt_calc-Cu_Liq)**2
    return Residual
## ONeill parameterization for sulfide composition
def calculate_ONeill_sulf(FeOt_Liq, Ni_Liq, Cu_Liq, MgO_Liq=None, Fe3Fet_Liq=None):
    '''
    Calculating the Fe_FeNiCu ratio in the sulfide using the empirical
    parameterizatin of ONeill.

    '''

    Fe_FeNiCu=1/(1+(Ni_Liq/(FeOt_Liq*(1-Fe3Fet_Liq)))*0.013+(Cu_Liq/(FeOt_Liq*(1-Fe3Fet_Liq)))*0.025)
    return Fe_FeNiCu
## Calculating sulfide stuff once minimised

def calculate_Kiseeva_sulf_comp_kd(Ni_Sulf, Cu_Sulf, FeOt_Liq, Fe3Fet_Liq, T_K,  Ni_Liq, Cu_Liq):

    FeO_Liq=FeOt_Liq*(1-Fe3Fet_Liq)
    OCalc_CellAG12=0.24*FeO_Liq
    FeS_calc_AL19=((100-Ni_Sulf-Cu_Sulf-OCalc_CellAG12-Cu_Sulf*20.1442646/79.8557354
    -Ni_Sulf*35.32650016/64.67349984))*36.47119049/100
    FeCalc_CellAF12=FeS_calc_AL19*63.52880951/36.47119049

    AG6=Ni_Sulf/58.6934/(Ni_Sulf/58.6934+Cu_Sulf/63.546+FeCalc_CellAF12/55.845)
    AG8=AG6*0.97
    AG9=AG6*0.92

        # If FeO<13
    O_Sulf=(FeO_Liq*0.24*((1-AG8)**2))
    O_Sulf.loc[FeO_Liq>=13]=(FeO_Liq*0.24*((1-AG9)**2))

    AL17=(Ni_Sulf*35.32650016/64.67349984)
    AL18=(Cu_Sulf*20.1442646/79.8557354)
    AM20=O_Sulf*77.72983506/22.27016494

    S_Sulf=((100-AL17-AL18-(O_Sulf*77.72983506/22.27016494)
             -O_Sulf-Cu_Sulf-Ni_Sulf)*36.47119049/100 + AL17 + AL18)

    Fe_Sulf=((100-Ni_Sulf-Cu_Sulf-O_Sulf-AL17-AL18-AM20)
    *(63.52880951/100)+O_Sulf*77.72983506/22.27016494)

    NiSS=Ni_Sulf*(35.32650016/64.67349984)
    CuS05=Cu_Sulf*20.1442646/79.8557354
    OCalc=O_Sulf*77.72983506/22.27016494
    S_Sulf=(100-NiSS-CuS05-OCalc-O_Sulf-Cu_Sulf-Ni_Sulf)*36.47119049/100+NiSS+CuS05

    # Cu Coefficients
    A_Cu=0.6995
    B_Cu=4200
    eps_Cu_Cu=-1.0310
    eps_FeO_Cu=2.2030
    eps_Ni_Cu=0.00
    # Ni Coefficients

    A_Ni=1.869
    B_Ni=3300
    eps_Cu_Ni=0.00
    eps_Ni_Ni=-0.904
    eps_FeO_Ni=0




    Fe_Corr=FeO_Liq/(Fe_Sulf/55.845/(Ni_Sulf/58.6934+Cu_Sulf/63.546+Fe_Sulf/55.845)) #AF7
    Ni_FeNiCu=Ni_Sulf/58.6934/(Ni_Sulf/58.6934+Cu_Sulf/63.546+Fe_Sulf/55.845)
    Cu_FeNiCu=Cu_Sulf/63.546/(Ni_Sulf/58.6934+Cu_Sulf/63.546+Fe_Sulf/55.845)
    Fe_FeNiCu=Fe_Sulf/55.845/(Ni_Sulf/58.6934+Cu_Sulf/63.546+Fe_Sulf/55.845)

    DNi=(10**(B_Ni/T_K+A_Ni-np.log10(Fe_Corr)+(1673*eps_FeO_Ni*np.log10(1-0.049*O_Sulf)/T_K)
            +(1673*eps_Ni_Ni*np.log10(1-Ni_FeNiCu)/T_K)+(1673*eps_Cu_Ni*np.log10(1-Cu_FeNiCu)/T_K)))

    DCu=(10**(B_Cu/T_K+A_Cu-0.5*np.log10(Fe_Corr)+(1673*eps_FeO_Cu*np.log10(1-0.049*O_Sulf)/T_K)
             +(1673*eps_Cu_Cu*np.log10(1-Cu_FeNiCu)/T_K)+(1673*eps_Ni_Cu*np.log10(1-Ni_FeNiCu)/T_K)))

    Ni_Melt_calc=Ni_Sulf*10000/DNi
    Cu_Melt_calc=Cu_Sulf*10000/DCu

    df_out=pd.DataFrame(data={
    'S_Sulf': S_Sulf, 'O_Sulf': O_Sulf, 'Fe_Sulf': Fe_Sulf, 'Ni_Sulf': Ni_Sulf, 'Cu_Sulf': Cu_Sulf})
    return df_out





## Calculating sulfide composition once you have iterated towards a respond.


def calculate_Symthe_Sulf_Composition(FeOt_Liq, Fe3Fet_Liq, T_K, Ni_Liq, Cu_Liq):
    Ni_Sulf=np.empty(len(FeOt_Liq), dtype=float)
    Cu_Sulf=np.empty(len(FeOt_Liq), dtype=float)
    bnds=((0, 30), (0, 30))

    for i in range(0, len(FeOt_Liq)):
        Calc_Sulf=scipy.optimize.minimize(Loop_Smythe_Sulf_Calc_Residual, x0=(0, 0), bounds=bnds,
                        args=(FeOt_Liq.iloc[i]*(1-Fe3Fet_Liq.iloc[i]),
                              T_K.iloc[i], Ni_Liq.iloc[i], Cu_Liq.iloc[i])).get('x')


        Ni_Sulf[i]=Calc_Sulf[0]
        Cu_Sulf[i]=Calc_Sulf[1]
    df_out=pd.DataFrame(data={'Ni_Sulf': Ni_Sulf, 'Cu_Sulf': Cu_Sulf})


    return df_out


def calculate_sulf_FeFeNiCu(Ni_Sulf, Cu_Sulf, Fe_Sulf):
    Ni_moles=Ni_Sulf/58.6934
    Cu_moles=Cu_Sulf/63.546
    Fe_moles=Fe_Sulf/55.8457
    FeFeNiCu=Fe_moles/(Fe_moles+Cu_moles+Ni_moles)
    return FeFeNiCu

def calculate_sulf_CuFeNiCu(Ni_Sulf, Cu_Sulf, Fe_Sulf):
    Ni_moles=Ni_Sulf/58.6934
    Cu_moles=Cu_Sulf/63.546
    Fe_moles=Fe_Sulf/55.8457
    CuFeNiCu=Cu_moles/(Fe_moles+Cu_moles+Ni_moles)
    return CuFeNiCu

def calculate_sulf_NiFeNiCu(Ni_Sulf, Cu_Sulf, Fe_Sulf):
    Ni_moles=Ni_Sulf/58.6934
    Cu_moles=Cu_Sulf/63.546
    Fe_moles=Fe_Sulf/55.8457
    NiFeNiCu=Ni_moles/(Fe_moles+Cu_moles+Ni_moles)
    return NiFeNiCu
## Smythe Parameterization

def calculate_Smythe2017_SCSS(*, df, T_K, P_kbar, Fe3Fet_Liq, Fe_FeNiCu_Sulf=None,
Cu_FeNiCu_Sulf=None, Ni_FeNiCu_Sulf=None, Fe_Sulf=None, Ni_Sulf=None, Cu_Sulf=None, Ni_Liq=None, Cu_Liq=None):
    '''
    Calculates SCSS using the model of Smythe et al. (2017), with options for users to
    calculate sulfide composition from liquid composition, or input sulfide composition.

    Parameters
    -------
    df: pandas.DataFrame
        Dataframe of liquid compositions. Needs to have the headings "SiO2_Liq", "TiO2_Liq" etc, with
        compositions in oxide wt%. all FeO should be entered as FeOt_Liq, which can then be partitioned
        using the Fe3Fet_Liq input. Heading order doesn't matter.

    T_K: int, float, pandas.Series
        Temperature in Kelvin.

    P_kbar: int, float, pandas.Series
        Pressure in kbar

    Fe3Fet_Liq: int, float, pandas.Series
        Proportion of Fe3+ in the liquid.
        Various parts of the calculation use only Fe2.

    Sulfide Composition: Options to calculate from the liquid composition, enter the comp in el wt%,
    or enter the FeFeNiCu, Cu

      if you want to calculate sulfide composition:

            Fe_FeNiCu_Sulf = "Calc_Smythe"
            Calculates sulfide composition analagous to the Solver function in the Smythe et al. (2017) spreadsheet.
            Here, we use Scipy optimize to find the ideal Ni and Cu
            contents using Kds from Kiseeva et al. (2015) and the Ni and Cu content of the melt.
            Requires user to also enter Ni_Liq and Cu_Liq.

            Or

            Fe_FeNiCu_Sulf = "Calc_ONeill"
            Calculates sulfide composition using the empirical expression of O'Neill (2021), which depends on
            FeOt_Liq, Ni_Liq, Cu_Liq, and Fe3Fet_Liq. We allow users to enter their own Fe3Fet_Liq,
            as we believe the empirical model of Neill where Fe3Fet_Liq is a function of MgO content is not
            broadly applicable.

        if you want to input a Fe_FeNiCu_Sulf ratio:
            Fe_FeNiCu_Sulf = int, float, pandas series
            Calculates SCSS using this ratio.
            If you want the non-ideal SCSS to be returned, you also need to enter
            values for Cu_FeNiCu_Sulf and Ni_FeNiCu_Sulf

        if you want to input a measured sulfide composition in el wt%
            Fe_Sulf, Ni_Sulf, Cu_Sulf = int, float, pandas series




    Returns
    -------
    pandas.DataFrame:
        Contains column for SCSS ideal, 1 sigma, input T and P, and the various intermediate steps of the calculation.

    '''

    # First, calculate silicate hydrous mole fractions, as true regardless of choice of sulfide composition
    Smythe_calcs=calculate_Smythe_silicate_mole_fractions(df, Fe3Fet_Liq)

    # Then check options, basically need to have entered FeFeNiCu_Calc, or Fe_FeNiCu ratio, or all sulfide
    if Fe_Sulf is not None and Ni_Sulf is not None and Cu_Sulf is not None and Fe_FeNiCu_Sulf is not None:
        raise ValueError('You have entered both a Fe_FeNiCu_Sulf ratio, and the conc of Fe, Ni and Cu in your sulf. Please enter one set of inputs or another')

    elif isinstance(Fe_FeNiCu_Sulf, str) and (Fe_FeNiCu_Sulf=="Calc_Smythe" or Fe_FeNiCu_Sulf == "Calc_ONeill"):
        if Ni_Liq is None:
            raise ValueError('If you select Fe_FeNiCu_Sulf = Calc, you need to enter the concentration of Cu and Ni in the liquid in ppm')
        if Cu_Liq is None:
            raise ValueError('If you select Fe_FeNiCu_Sulf = Calc, you need to enter the concentration of Cu and Ni in the liquid in ppm')
        print('Calculating Sulf composition using Scipy minimisation of Kiseeva et al. (2015) Kd Values for Ni and Cu')

    elif Fe_FeNiCu_Sulf is not None:
        print('Using inputted Fe_FeNiCu_Sulf ratio for calculations.')

    elif Fe_Sulf is not None and Ni_Sulf is not None and Cu_Sulf is not None:
        print('Using inputted Sulf compositions to calculate Fe_FeNiCu_Sulf ratios.')

    else:
        raise ValueError('Input for sulfide composition not recognised.')


    if isinstance(Fe_FeNiCu_Sulf, str):
        if Fe_FeNiCu_Sulf=="Calc_Smythe":


            # This does the Scipy minimisation of Cu and Ni contents using Kiseeva et al. (2015)
            calc_sulf=calculate_Symthe_Sulf_Composition(FeOt_Liq=df['FeOt_Liq'], Fe3Fet_Liq=df['Fe3Fet_Liq'],
                                            T_K=T_K, Ni_Liq=Ni_Liq, Cu_Liq=Cu_Liq)

            # This feeds those result back into a simpler function to get the Fe, S and O content of the sulfide
            Sulf_All=calculate_Kiseeva_sulf_comp_kd(Ni_Sulf=calc_sulf['Ni_Sulf'],Cu_Sulf=calc_sulf['Cu_Sulf'],
                                FeOt_Liq=df['FeOt_Liq'], Fe3Fet_Liq=Fe3Fet_Liq,
                                            T_K=T_K, Ni_Liq=df['Ni_Liq (ppm)'], Cu_Liq=df['Cu_Liq (ppm)'])

            Fe_FeNiCu_Sulf=calculate_sulf_FeFeNiCu(Sulf_All['Ni_Sulf'], Sulf_All['Cu_Sulf'], Sulf_All['Fe_Sulf'])
            Cu_FeNiCu_Sulf=calculate_sulf_CuFeNiCu(Sulf_All['Ni_Sulf'], Sulf_All['Cu_Sulf'], Sulf_All['Fe_Sulf'])
            Ni_FeNiCu_Sulf=calculate_sulf_NiFeNiCu(Sulf_All['Ni_Sulf'], Sulf_All['Cu_Sulf'], Sulf_All['Fe_Sulf'])

            Smythe_calcs['Ni_Sulf_Calc']=Sulf_All['Ni_Sulf']
            Smythe_calcs['Cu_Sulf_Calc']=Sulf_All['Cu_Sulf']
            Smythe_calcs['Fe_Sulf_Calc']=Sulf_All['Fe_Sulf']
            Smythe_calcs['O_Sulf_Calc']=Sulf_All['O_Sulf']
            Smythe_calcs['S_Sulf_Calc']=Sulf_All['S_Sulf']

        if isinstance(Fe_FeNiCu_Sulf, str) and Fe_FeNiCu_Sulf=="Calc_ONeill":
            Fe_FeNiCu_Sulf=calculate_ONeill_sulf(FeOt_Liq=df['FeOt_Liq'],
            Ni_Liq=df['Ni_Liq (ppm)'], Cu_Liq=df['Cu_Liq (ppm)'], Fe3Fet_Liq=Fe3Fet_Liq)

    if Fe_Sulf is not None and Ni_Sulf is not None and Cu_Sulf is not None:
            Fe_FeNiCu_Sulf=calculate_sulf_FeFeNiCu(Ni_Sulf, Cu_Sulf, Fe_Sulf)
            Cu_FeNiCu_Sulf=calculate_sulf_CuFeNiCu(Ni_Sulf, Cu_Sulf, Fe_Sulf)
            Ni_FeNiCu_Sulf=calculate_sulf_NiFeNiCu(Ni_Sulf, Cu_Sulf, Fe_Sulf)



    # Calculating the different liquid components
    Smythe_calcs['Si_XA_ideal']=Smythe_calcs['Si_wt_atom']*(-27561.044)
    Smythe_calcs['Ti_XA_ideal']=Smythe_calcs['Ti_wt_atom']*(-11220.488)
    Smythe_calcs['Al_XA_ideal']=Smythe_calcs['Al_wt_atom']*(-18450.292)
    Smythe_calcs['Mg_XA_ideal']=Smythe_calcs['Mg_wt_atom']*(-13969.67)
    Smythe_calcs['Fe2_XA_ideal']=Smythe_calcs['Fe2_wt_atom']*(-34274.174)
    Smythe_calcs['Ca_XA_ideal']=Smythe_calcs['Ca_wt_atom']*(-7830.853)
    Smythe_calcs['Na_XA_ideal']=Smythe_calcs['Na_wt_atom']*(-13246.751)
    Smythe_calcs['K_XA_ideal']=Smythe_calcs['K_wt_atom']*(-29014.575)
    Smythe_calcs['H_XA_ideal']=Smythe_calcs['H_wt_atom']*(-17495.266)
    Smythe_calcs['Si*Fe_ideal']=(Smythe_calcs['Si_wt_atom']*Smythe_calcs['Fe2_wt_atom'])*116567.625


    Smythe_calcs['Si_XA_non_ideal']=Smythe_calcs['Si_wt_atom']*(-27996.431)
    Smythe_calcs['Ti_XA_non_ideal']=Smythe_calcs['Ti_wt_atom']*(-10714.991)
    Smythe_calcs['Al_XA_non_ideal']=Smythe_calcs['Al_wt_atom']*(-18999.945)
    Smythe_calcs['Mg_XA_non_ideal']=Smythe_calcs['Mg_wt_atom']*(-14512.488)
    Smythe_calcs['Fe2_XA_non_ideal']=Smythe_calcs['Fe2_wt_atom']*(-34895.294)
    Smythe_calcs['Ca_XA_non_ideal']=Smythe_calcs['Ca_wt_atom']*(-8831.616)
    Smythe_calcs['Na_XA_non_ideal']=Smythe_calcs['Na_wt_atom']*(-13712.715)
    Smythe_calcs['K_XA_non_ideal']=Smythe_calcs['K_wt_atom']*(-28583.983)
    Smythe_calcs['H_XA_non_ideal']=Smythe_calcs['H_wt_atom']*(-17766.114)
    Smythe_calcs['Si*Fe_non_ideal']=(Smythe_calcs['Si_wt_atom']*Smythe_calcs['Fe2_wt_atom'])*117815.515


    Smythe_calcs['log_SCSS_ideal']=(

(122175-80.28*T_K+8.474*T_K*np.log(T_K))/(8.314*T_K)+9.087+(Smythe_calcs['Si_XA_ideal']+Smythe_calcs['Ti_XA_ideal']
+Smythe_calcs['Al_XA_ideal']+Smythe_calcs['Mg_XA_ideal']+Smythe_calcs['Fe2_XA_ideal']+Smythe_calcs['Ca_XA_ideal']
+Smythe_calcs['Na_XA_ideal']+Smythe_calcs['K_XA_ideal']+Smythe_calcs['H_XA_ideal']+Smythe_calcs['Si*Fe_ideal'])/T_K
+np.log(Fe_FeNiCu_Sulf)-np.log(Smythe_calcs['Fe2_wt_atom'])-269.4*0.1*P_kbar/T_K


    )


    Smythe_calcs['SCSS_ideal_ppm_Smythe2017']=np.exp(Smythe_calcs['log_SCSS_ideal'])
    Smythe_calcs['SCSS_ideal_ppm_Smythe2017_1sigma']=Smythe_calcs['SCSS_ideal_ppm_Smythe2017']*0.273169775211857
    Smythe_calcs['T_Input_K']=T_K
    Smythe_calcs['P_Input_kbar']=P_kbar
    Smythe_calcs['Fe_FeNiCu_Sulf']=Fe_FeNiCu_Sulf
    Smythe_calcs['Fe3Fet_Liq_input']=Fe3Fet_Liq


    if Cu_FeNiCu_Sulf is not None and Ni_FeNiCu_Sulf is not None:
        Smythe_calcs['log_SCSS_non_ideal']=(
        (122175-80.28*T_K+8.474*T_K*np.log(T_K))/(8.314*T_K)+9.352+(Smythe_calcs['Si_XA_non_ideal']+Smythe_calcs['Ti_XA_non_ideal']
    +Smythe_calcs['Al_XA_non_ideal']+Smythe_calcs['Mg_XA_non_ideal']+Smythe_calcs['Fe2_XA_non_ideal']+Smythe_calcs['Ca_XA_non_ideal']
    +Smythe_calcs['Na_XA_non_ideal']+Smythe_calcs['K_XA_non_ideal']+Smythe_calcs['H_XA_non_ideal']+Smythe_calcs['Si*Fe_non_ideal'])
    /T_K+np.log(Fe_FeNiCu_Sulf)-np.log(Smythe_calcs['Fe2_wt_atom'])-264.85*0.1*P_kbar/T_K+546.362*((Cu_FeNiCu_Sulf**2 +Cu_FeNiCu_Sulf*Ni_FeNiCu_Sulf)/T_K)
        )
        Smythe_calcs['SCSS_non_ideal_ppm_Smythe2017']=np.exp(Smythe_calcs['log_SCSS_non_ideal'])
        Smythe_calcs['SCSS_non_ideal_ppm_Smythe2017_1sigma']=Smythe_calcs['SCSS_non_ideal_ppm_Smythe2017']*0.267299081373473

        cols_to_move = ['SCSS_ideal_ppm_Smythe2017', 'SCSS_ideal_ppm_Smythe2017_1sigma',
        'SCSS_non_ideal_ppm_Smythe2017', 'SCSS_non_ideal_ppm_Smythe2017_1sigma',
        'T_Input_K', "P_Input_kbar",'Fe_FeNiCu_Sulf']

    else:
        print('You havent entered a value for Ni_FeNiCu_Sulf and Cu_FeNiCu_Sulf so we cant calculate the non-ideal SCSS')
        cols_to_move = ['SCSS_ideal_ppm_Smythe2017', 'SCSS_ideal_ppm_Smythe2017_1sigma',
        'T_Input_K', "P_Input_kbar",'Fe_FeNiCu_Sulf']

    Smythe_calcs = Smythe_calcs[cols_to_move +
                                    [col for col in Smythe_calcs.columns if col not in cols_to_move]]


    return Smythe_calcs












