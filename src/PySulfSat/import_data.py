import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.optimize as optimize
from scipy.special import erf

## Import data

df_ideal_liq = pd.DataFrame(columns=['SiO2_Liq', 'TiO2_Liq', 'Al2O3_Liq',
'FeOt_Liq', 'MnO_Liq', 'MgO_Liq', 'CaO_Liq', 'Na2O_Liq', 'K2O_Liq',
'P2O5_Liq', 'H2O_Liq', 'Fe3Fet_Liq', 'Ni_Liq_ppm', 'Cu_Liq_ppm'])

def import_data(filename, sheet_name=None, Petrolog=False, MELTS=False, sample_label=None, suffix=None):
    if 'csv' in filename:
        my_input = pd.read_csv(filename)

    if Petrolog is True:

        df=pd.read_excel(filename, sheet_name='Petrolog_Output_FRAC')
        df.columns= df.columns.str.replace(' ','',regex=True)
        df.columns= df.columns.str.replace('_melt','_Liq',regex=True)
        if sum(df.columns.str.contains('Ni_Liq'))>0:
            df['Ni_Liq_ppm']=df['Ni_Liq']
        else:
            df['Ni_Liq_ppm']=df['Ni_Liq']
            print('We didnt find a Ni column in your Petrolog input')

        if sum(df.columns.str.contains('Cu_Liq'))>0:
            df['Cu_Liq_ppm']=df['Cu_Liq']
        else:
            df['Cu_Liq_ppm']=df['Cu_Liq']
            print('We didnt find a Ni column in your Petrolog input')

        df['FeOt_Liq']=df['FeO_Liq']+df['Fe2O3_Liq']*0.899

        98
        df['Fe3Fet_Liq']=1-(df['FeO_Liq']/df['FeOt_Liq'])
        df2=df.drop(['FeO_Liq', 'Fe2O3_Liq'], axis=1)
        df2['T_K']=df2['Temperature']+273.15
        df2['P_kbar']=df2['Pressure(kbar)']
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


