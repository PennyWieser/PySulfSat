import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.optimize as optimize
from scipy.special import erf
from io import StringIO

## Import data

df_ideal_liq = pd.DataFrame(columns=['SiO2_Liq', 'TiO2_Liq', 'Al2O3_Liq',
'FeOt_Liq', 'MnO_Liq', 'MgO_Liq', 'CaO_Liq', 'Na2O_Liq', 'K2O_Liq',
'P2O5_Liq', 'H2O_Liq', 'Fe3Fet_Liq', 'Ni_Liq_ppm', 'Cu_Liq_ppm'])

def import_dataframe(df, suffix=None):
    """ This function takes a user dataframe, and reforms the columns into the format required by PySulfSat,
    In many cases this involves renaming columns to get into the format SiO2_Liq, TiO2_Liq, etc

    Parameters
    --------------

    df: pandas Dataframe

    suffix: str
        Suffix which is appended onto all columns. E.g. if you have EPMA data, and you cant be bothered to add
    '_Liq' to each column manually, it will do so here (but for all columns, not just the oxides)



    """
    my_input=df

    if suffix is not None:
        if any(my_input.columns.str.contains(suffix)):
            w.warn('We notice you have specified a suffix, but some of your columns already have this suffix. '
        'e.g., If you already have _Liq in the file, you shouldnt specify suffix="_Liq" during the import')

    my_input_c = my_input.copy()

    if "Sample_ID_Liq" not in my_input_c:
        my_input_c['Sample_ID_Liq'] = my_input.index

    if suffix is not None:
        if any(my_input.columns.str.contains("FeOT")) and (all(my_input.columns.str.contains("FeOt")==False)):
            raise ValueError("No FeOt column found. You've got a column heading with FeOT. Change to a lower case t")

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
    #my_input_c=my_input.copy()
    for col in cols2:
        if col in my_input_c.columns:
            my_input_c=my_input_c.drop(columns=col)

    out=pd.concat([myLiquids1, my_input_c], axis=1)
    return out





def import_data(filename, sheet_name=None, Petrolog=False, MELTS=False, MELTS_txt=False, suffix=None):
    """ This function takes a user input, and reforms the columns into the format required by PySulfSat,
    In many cases this involves renaming columns to get into the format SiO2_Liq, TiO2_Liq, etc

    Parameters
    --------------
    filename: str
        File name (e.g. 'Test1.xlsx')


    Petrolog: bool
        True if output from Petrolog3 software, False (default) if not

    MELTS: bool
        True if loading a MELTS tbl file, that is, a file that ends in .tbl

    MELTS_txt: bool
        True if loading a MELTS file that has tbl in its name, but is actually  a txt file.



    Optional:
    sheet_name: str
        Name of sheet if filename refers to an Excel spreadsheet ('Sheet1')

    suffix: str
        Suffix which is appended onto all columns. E.g. if you have EPMA data, and you cant be bothered to add
    '_Liq' to each column manually, it will do so here (but for all columns, not just the oxides)



    Returns
    --------------
    pd.DataFrame
        DataFrame of all input data, renamed columns for liquids at start, other columns appended onto the back.



    """

    if 'csv' in filename:
        my_input = pd.read_csv(filename)

    if Petrolog is True:

        df=pd.read_excel(filename, sheet_name='Petrolog_Output_FRAC')
        df.columns= df.columns.str.replace(' ','',regex=True)
        df.columns= df.columns.str.replace('_melt','_Liq',regex=True)
        if sum(df.columns.str.contains('Ni_Liq'))>0:
            df['Ni_Liq_ppm']=df['Ni_Liq']
        else:
            df['Ni_Liq_ppm']=0
            print('We didnt find a Ni column in your Petrolog input')

        if sum(df.columns.str.contains('Cu_Liq'))>0:
            df['Cu_Liq_ppm']=df['Cu_Liq']
        else:
            df['Cu_Liq_ppm']=0
            print('We didnt find a Ni column in your Petrolog input')

        df['FeOt_Liq']=df['FeO_Liq']+df['Fe2O3_Liq']*0.89998
        df['Fe3Fet_Liq']=1-(df['FeO_Liq']/df['FeOt_Liq'])
        df2=df.drop(['FeO_Liq', 'Fe2O3_Liq'], axis=1)
        df2['T_K']=df2['Temperature']+273.15
        df2['P_kbar']=df2['Pressure(kbar)']

        if 'Melt_%_magma' in df2.columns:
            df2['Fraction_melt'] = df2['Melt_%_magma'] / 100
        my_input=df2


    elif MELTS is True:
        df=pd.read_table(filename, sep=',')
        df.columns=df.columns.str.replace('wt% ','',regex=True)
        df2=df.rename(columns={
                                'SiO2': 'SiO2_Liq',
                                'TiO2': 'TiO2_Liq',
                                'CaO': 'CaO_Liq',
                                'Al2O3': 'Al2O3_Liq',
                                'FeO': 'FeO_Liq',
                                'Fe2O3': 'Fe2O3_Liq',
                                'Cr2O3': 'Cr2O3_Liq',
                                'MnO': 'MnO_Liq',
                                'MgO': 'MgO_Liq',
                                'NiO': 'NiO_Liq',
                                'Na2O': 'Na2O_Liq',
                                'K2O': 'K2O_Liq',
                                'P2O5': 'P2O5_Liq',
                                'H2O': 'H2O_Liq'
                                })

        df2['FeOt_Liq']=df2['FeO_Liq']+df2['Fe2O3_Liq']*0.89998
        df2['Fe3Fet_Liq']=1-(df2['FeO_Liq']/df2['FeOt_Liq'])
        df2['T_K']=df2['T (C)']+273.15
        df2['P_kbar']=df2['P (kbars)']
        my_input=df2

    elif MELTS_txt is True:
        from io import StringIO
        with open(filename, 'r') as file:
            lines = file.readlines()

        # Find the index where the "Pressure" word is located
        pressure_index = 0
        for i, line in enumerate(lines):
            if "Pressure" in line:
                pressure_index = i
                break

        # Create a new list of lines starting from the "Pressure" line
        filtered_lines = lines[pressure_index:]

        # Create a Pandas DataFrame from the filtered lines
        df = pd.read_csv(StringIO(''.join(filtered_lines)), delim_whitespace=True)

        # Create a new DataFrame with the desired columns
        df2 = df[['SiO2', 'TiO2', 'CaO', 'Al2O3', 'FeO', 'Fe2O3', 'Cr2O3', 'MnO', 'MgO', 'Na2O', 'K2O', 'P2O5', 'H2O']].copy()

        # Rename the columns in the new DataFrame
        df2.columns = ['SiO2_Liq', 'TiO2_Liq', 'CaO_Liq', 'Al2O3_Liq', 'FeO_Liq', 'Fe2O3_Liq', 'Cr2O3_Liq', 'MnO_Liq', 'MgO_Liq',  'Na2O_Liq', 'K2O_Liq', 'P2O5_Liq', 'H2O_Liq']


        df2['FeOt_Liq']=df['FeO']+df['Fe2O3']*0.89998
        df2['Fe3Fet_Liq']=1-(df2['FeO_Liq']/df2['FeOt_Liq'])
        df2['T_K']=df['Temperature']+273.15
        df2['P_kbar']=df['Pressure']/1000
        my_input=df2






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

    if "Sample_ID_Liq" not in my_input_c:
        my_input_c['Sample_ID_Liq'] = my_input.index

    if suffix is not None:
        if any(my_input.columns.str.contains("FeOT")) and (all(my_input.columns.str.contains("FeOt")==False)):
            raise ValueError("No FeOt column found. You've got a column heading with FeOT. Change to a lower case t")

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
    #my_input_c=my_input.copy()
    for col in cols2:
        if col in my_input_c.columns:
            my_input_c=my_input_c.drop(columns=col)

    out=pd.concat([myLiquids1, my_input_c], axis=1)
    return out



df_ideal_liq_noise = pd.DataFrame(columns=['SiO2_Liq_Err', 'TiO2_Liq_Err', 'Al2O3_Liq_Err',
'FeOt_Liq_Err', 'MnO_Liq_Err', 'MgO_Liq_Err', 'CaO_Liq_Err', 'Na2O_Liq_Err', 'K2O_Liq_Err',
'P2O5_Liq_Err', 'H2O_Liq_Err', 'Fe3Fet_Liq_Err', 'Ni_Liq_ppm_Err', 'Cu_Liq_ppm_Err'])

def import_data_noise(filename, sheet_name=None, sample_label=None):
    """ This function takes a user input for errors, and reforms the columns into the format required by PySulfSat,
    Input columns should have _Err after them.

    Parameters
    --------------
    filename: str
        File name (e.g. 'Test1.xlsx')

    sheet_name: str
        Name of sheet if filename refers to an Excel spreadsheet ('Sheet1')



    Returns
    --------------
    pd.DataFrame
        DataFrame of all input data, renamed columns for liquids at start, other columns appended onto the back.



    """
    if 'csv' in filename:
        my_input = pd.read_csv(filename)




    if 'xls' in filename:
        if sheet_name is not None:
            my_input = pd.read_excel(filename, sheet_name=sheet_name)
            #my_input[my_input < 0] = 0
        else:
            my_input = pd.read_excel(filename)
            #my_input[my_input < 0] = 0



    my_input_c = my_input.copy()

    if "Sample_ID_Liq" not in my_input_c:
        my_input_c['Sample_ID_Liq'] = my_input.index



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



    myLiquids1 = my_input_c.reindex(df_ideal_liq_noise.columns, axis=1).fillna(0)
    myLiquids1 = myLiquids1.apply(pd.to_numeric, errors='coerce').fillna(0)
    myLiquids1[myLiquids1 < 0] = 0
    print('We have replaced all missing liquid oxides and strings with zeros. ')

    cols2=myLiquids1.columns
    #my_input_c=my_input.copy()
    for col in cols2:
        if col in my_input_c.columns:
            my_input_c=my_input_c.drop(columns=col)

    out=pd.concat([myLiquids1, my_input_c], axis=1)
    return out



