import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.optimize as optimize
from scipy.special import erf
import warnings as w


from PySulfSat.core_calcs import *

## Li and Zhang 2022
df_ideal_liq_lizhang = pd.DataFrame(columns=['SiO2_Liq', 'TiO2_Liq', 'Al2O3_Liq',
'FeOt_Liq', 'MnO_Liq', 'MgO_Liq', 'CaO_Liq', 'Na2O_Liq', 'K2O_Liq',
'P2O5_Liq'])

def norm_liqs_excl_H2O(Liqs):
    Liqs_c=Liqs.copy()
    Liqs_c=Liqs_c.fillna(0)
    Liqs2=Liqs_c.reindex(
        df_ideal_liq_lizhang.columns, axis=1).fillna(0)

    sum_rows=Liqs2.sum(axis=1)

    Liqs_norm=Liqs2.divide(sum_rows/100, axis='rows')
    return Liqs_norm



def calculate_LiZhang2022_SCSS(*, df, T_K, P_kbar, H2O_Liq=None, Fe_FeNiCu_Sulf=None, Cu_FeNiCu_Sulf=None, Ni_FeNiCu_Sulf=None, Fe3Fet_Liq=None, logfo2=None,
Ni_Liq=None, Cu_Liq=None, Ni_Sulf_init=5, Cu_Sulf_init=5, Fe_Sulf=None, Cu_Sulf=None, Ni_Sulf=None, T_K_transition=True,
highT=False, lowT=False):

    '''
    Calculates SCSS using the model of Liu and Zhang (2022), doi: https://doi.org/10.1016/j.gca.2022.03.0080
    with options for users to
    calculate sulfide composition from liquid composition using the approach of Smythe et al. (2017),
    the empirical parameterization of O'Neill (2021), or input a sulfide composition.

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

    logfo2: int, float, or pandas.Series
        logfo2, used to convert to Fe3 following the method of Li and Zhang (2022)

    T_K_transition : bool
        The model shows a flip at 1200C. If T_K_transition is False, this doesnt happen. This generates more coherent results.

    Fe3Fet_Liq: int, float, pandas.Series, or "Calc_ONeill"
        Proportion of Fe3+ in the liquid, as various parts of the calculation use only Fe2.
        If "Calc_ONeill" Calculates as a function of MgO content, using the MORB relationship.

    Sulfide Composition: Options to calculate from the liquid composition, enter the comp in el wt%,
    or enter the Fe_FeNiCu_Sulf, Cu

      if you want to calculate sulfide composition:

            Fe_FeNiCu_Sulf = "Calc_Smythe", also needs Ni_Sulf_init, and Cu_Sulf_init
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
        Contains column for SCSS ideal, and user inputted dataframe

    '''



    T_K2=T_K-0.15 # As they use 273, not 273.15 for their conversion

    Liqs=df.copy()

    if H2O_Liq is not None:
        Liqs['H2O_Liq']=H2O_Liq
    if Fe3Fet_Liq is not None:
        Liqs['Fe3Fet_Liq']=Fe3Fet_Liq
    else:
        Liqs['Fe3Fet_Liq']=0


    Liqs_norm=norm_liqs_excl_H2O(Liqs)

    Liqs=calculate_sulfide_comp_generic(
    Fe_Sulf=Fe_Sulf, Ni_Sulf=Ni_Sulf, Cu_Sulf=Cu_Sulf, Fe_FeNiCu_Sulf=Fe_FeNiCu_Sulf,
    Cu_FeNiCu_Sulf=Cu_FeNiCu_Sulf, Ni_FeNiCu_Sulf=Ni_FeNiCu_Sulf,
    Ni_Liq=Ni_Liq, Cu_Liq=Cu_Liq, df_c=Liqs, T_K=T_K)

    moles=calculate_anhydrous_mol_proportions_liquid(liq_comps=Liqs_norm)
    mol_frac=calculate_anhydrous_mol_fractions_liquid(liq_comps=Liqs_norm)

    # Cation calculations

    if logfo2 is not None:
        logfo2=logfo2


        logXFe2O3XFeO=(0.196*logfo2/0.4343+11492/(T_K2)-6.675-2.243*mol_frac['Al2O3_Liq_mol_frac']
                    -1.828*mol_frac['FeOt_Liq_mol_frac']+3.201*mol_frac['CaO_Liq_mol_frac']
                    +5.854*mol_frac['Na2O_Liq_mol_frac']+6.215*mol_frac['K2O_Liq_mol_frac']
                    -3.36*(1-1673/(T_K2)-np.log((T_K2)/1673))-0.0701*(P_kbar*1000)/(T_K2)
                    -0.0000154*(P_kbar*1000)/(T_K2)*(T_K2-1673)+0.000000385*(P_kbar*1000)**2/(T_K2))
        FeO_mol_frac=mol_frac['FeOt_Liq_mol_frac']/(1+2*np.exp(logXFe2O3XFeO))
        Fe2O3_mol_frac=mol_frac['FeOt_Liq_mol_frac']/(1+2*np.exp(logXFe2O3XFeO))*np.exp(logXFe2O3XFeO)

        mol_cat=pd.DataFrame(data={'Si_cat': moles['SiO2_Liq_mol_prop'],
                            'Ti_cat': moles['TiO2_Liq_mol_prop'],
                                'Al_cat': 2*moles['Al2O3_Liq_mol_prop'],
                            'Fe_cat':moles['FeOt_Liq_mol_prop']/(1+2*np.exp(logXFe2O3XFeO)),
                            'Mn_cat': moles['MnO_Liq_mol_prop'],
                            'Mg_cat': moles['MgO_Liq_mol_prop'],
                            'Ca_cat': moles['CaO_Liq_mol_prop'],
                            'Na_cat':2* moles['Na2O_Liq_mol_prop'],
                            'K_cat': 2* moles['K2O_Liq_mol_prop'],
                            'P_cat': 2* moles['P2O5_Liq_mol_prop'],
                            'H_cat': 0,
                            'Fe3_cat':moles['FeOt_Liq_mol_prop']/(1+2*np.exp(logXFe2O3XFeO))*2*np.exp(logXFe2O3XFeO) })

    elif Fe3Fet_Liq is not None:
        mol_cat=pd.DataFrame(data={'Si_cat': moles['SiO2_Liq_mol_prop'],
                            'Ti_cat': moles['TiO2_Liq_mol_prop'],
                            'Al_cat': 2*moles['Al2O3_Liq_mol_prop'],
                            'Fe_cat':moles['FeOt_Liq_mol_prop']*(1-Fe3Fet_Liq),
                            'Mn_cat': moles['MnO_Liq_mol_prop'],
                            'Mg_cat': moles['MgO_Liq_mol_prop'],
                            'Ca_cat': moles['CaO_Liq_mol_prop'],
                            'Na_cat':2* moles['Na2O_Liq_mol_prop'],
                            'K_cat': 2* moles['K2O_Liq_mol_prop'],
                            'P_cat': 2* moles['P2O5_Liq_mol_prop'],
                            'H_cat': 0,
                            'Fe3_cat':moles['FeOt_Liq_mol_prop']*Fe3Fet_Liq})
    else:
        w.warn('The Li and Zhang (2022) model is sensitive to redox. we prefer you enter either logfo2, or Fe3Fet_Liq, we have performed calculations using Fe3Fet_Liq=0')
        mol_cat=pd.DataFrame(data={'Si_cat': moles['SiO2_Liq_mol_prop'],
                            'Ti_cat': moles['TiO2_Liq_mol_prop'],
                            'Al_cat': 2*moles['Al2O3_Liq_mol_prop'],
                            'Fe_cat':moles['FeOt_Liq_mol_prop']*(1-0),
                            'Mn_cat': moles['MnO_Liq_mol_prop'],
                            'Mg_cat': moles['MgO_Liq_mol_prop'],
                            'Ca_cat': moles['CaO_Liq_mol_prop'],
                            'Na_cat':2* moles['Na2O_Liq_mol_prop'],
                            'K_cat': 2* moles['K2O_Liq_mol_prop'],
                            'P_cat': 2* moles['P2O5_Liq_mol_prop'],
                            'H_cat': 0,
                            'Fe3_cat':moles['FeOt_Liq_mol_prop']*0})

    sum_mol_cat=mol_cat.sum(axis=1)
    mol_cat_norm=mol_cat.divide(sum_mol_cat, axis='rows')

    Fe3Fet_Liq=mol_cat_norm['Fe3_cat']/(mol_cat_norm['Fe_cat']+mol_cat_norm['Fe3_cat'])

    # print('Fe3_cat')
    # print(mol_cat_norm['Fe3_cat'])
    # print('Fe_cat')
    # print(mol_cat_norm['Fe_cat'])


    # Sulfide composition bit










    NaKAl=mol_cat_norm['Na_cat']+mol_cat_norm['K_cat']-mol_cat_norm['Al_cat']
    DeltaGRT=(137778-91.666*(T_K2)+8.474*(T_K2)*np.log(T_K2))/8.314/(T_K2)

    SumXMAM=(1673/(T_K2)*(6.7*(mol_cat_norm['Na_cat']+mol_cat_norm['K_cat'])
    +1.8*(mol_cat_norm['Al_cat']+mol_cat_norm['Fe3_cat'])+4.9*mol_cat_norm['Mg_cat']
    +8.1*mol_cat_norm['Ca_cat']+5*mol_cat_norm['Ti_cat']+8.9*(mol_cat_norm['Fe_cat']
    +mol_cat_norm['Mn_cat'])-22.2*(mol_cat_norm['Fe_cat']+
    mol_cat_norm['Mn_cat'])*mol_cat_norm['Ti_cat']+7.2*(mol_cat_norm['Fe_cat']
    +mol_cat_norm['Mn_cat'])*mol_cat_norm['Si_cat'])-2.06*erf(-7.2*(mol_cat_norm['Fe_cat']
    +mol_cat_norm['Mn_cat'])))

    lnCs=-23590/(T_K2)+8.77+SumXMAM
    lnXFeO=-np.log(mol_cat_norm['Fe_cat'])
    LnrFeO=-((1-mol_cat_norm['Fe_cat'])**2*(28870-14710*mol_cat_norm['Mg_cat']
    +1960*mol_cat_norm['Ca_cat']+43300*mol_cat_norm['Na_cat']+95380*mol_cat_norm['K_cat']
    -76880*mol_cat_norm['Ti_cat'])+(1-mol_cat_norm['Fe_cat'])*(-62190*mol_cat_norm['Si_cat']
    +31520*mol_cat_norm['Si_cat']**2))/8.314/(T_K2)
    # Need to sort out for low T<1200 C these two flip


    lnaFeS_lowT=-(31464-(T_K2)*21.506)/8.314/(T_K2)+np.log(Liqs['Fe_FeNiCu_Sulf_calc'])
    lnaFeS_HighT=np.log(1-mol_cat_norm['Fe_cat'])+np.log(Liqs['Fe_FeNiCu_Sulf_calc'])


    C1PC2erf_lowT=(-0.0291*(P_kbar)*1000+351*erf((P_kbar)*1000/10000))/(T_K2)+0.04*(P_kbar)*1000/8.314/(T_K2)
    C1PC2erf_highT=(-0.0291*(P_kbar)*1000+351*erf((P_kbar)*1000/10000))/(T_K2)


    if T_K_transition is True:
        C1PC2erf=C1PC2erf_highT
        C1PC2erf[T_K2<(1200+273)] =C1PC2erf_lowT
        lnaFeS=lnaFeS_HighT
        lnaFeS[T_K2<(1200+273)] =lnaFeS_lowT
    else:
        if highT is False and lowT is False:
            raise Exception('You have turned off the default T_K_transition, you now must decide if you want to use highT or lowT behavoir. Set highT or lowT=True')
        if highT is True:
            C1PC2erf=C1PC2erf_highT
            lnaFeS=lnaFeS_HighT
        if lowT is True:
            C1PC2erf=C1PC2erf_lowT
            lnaFeS=lnaFeS_lowT



    lnS=DeltaGRT+lnCs+lnXFeO+LnrFeO+lnaFeS+C1PC2erf

    S2_calc=np.exp(lnS)
    S2_calc[(Liqs['FeOt_Liq']<5)&(Liqs['H2O_Liq']>0)]=0

    SumMoles_H2O=(moles['FeOt_Liq_mol_prop']+moles['MnO_Liq_mol_prop']+moles['MgO_Liq_mol_prop']
                +moles['CaO_Liq_mol_prop']+moles['Na2O_Liq_mol_prop']+moles['K2O_Liq_mol_prop'])
    XH2Ot=(Liqs['H2O_Liq']/18)/(Liqs['H2O_Liq']/18+(moles['SiO2_Liq_mol_prop']*2+moles['TiO2_Liq_mol_prop']*2+
    moles['Al2O3_Liq_mol_prop']*3+SumMoles_H2O+moles['P2O5_Liq_mol_prop']*5)*(100-Liqs['H2O_Liq'])/100)


    lnXH2Ot=np.log(XH2Ot)
    KOH=np.exp(2.6*mol_cat_norm['Si_cat']-4339*mol_cat_norm['Si_cat']/(T_K2))
    XOH=(0.5-(0.25-(KOH-4)/KOH*(XH2Ot-XH2Ot**2))**0.5)/(KOH-4)*2*KOH
    lnXOH=np.log(XOH)
    XH2Om=XH2Ot-0.5*XOH
    lnXH2Om=np.log(XH2Om)
    lnXOH_XH2O=np.log(XOH+XH2Om)

    lnCHScalc=-19748*(1/T_K2)+7.81+SumXMAM+lnXOH_XH2O

    HScal=(np.exp(DeltaGRT+lnCHScalc+lnXFeO+LnrFeO+lnaFeS+C1PC2erf))/(100+Liqs['H2O_Liq'])*100

    # If NaKAl>-0.015
    NaKAlterm=(1673/(T_K2))*(19.634*NaKAl+0.2542)
    #Else
    LowNaKAl= NaKAl<-0.015
    NaKAlterm[LowNaKAl] =(1673/(T_K))*(26.365*NaKAl+0.9587)

    # If NaKAlterm>0
    lnCHS_NKA_term=lnCHScalc+NaKAlterm
    # Else
    lnCHS_NKA_term[NaKAlterm<=0]=lnCHScalc

    HScal2=(np.exp(DeltaGRT+lnCHS_NKA_term+lnXFeO+LnrFeO+lnaFeS+C1PC2erf))/(100+Liqs['H2O_Liq'])*100

    # If H2O>0

    # If HScal2>HScal
    Stot_calc_H2O=HScal2+S2_calc
    #If HScal2<=HScal
    calc2gcal=HScal2>HScal

    Stot_calc_H2O[~calc2gcal]=HScal+S2_calc



    # If H2O=0
    noH2O=~(Liqs['H2O_Liq']>0)
    Stot_calc_H2O[noH2O]=S2_calc

    #
    # Liqs.insert(1, 'Fe_FeNiCu_Sulf', Fe_FeNiCu_Sulf_calc)
    # Liqs.insert(2, 'HS cal', HScal)
    # Liqs.insert(3, 'S2_calc', S2_calc)
    # Liqs.insert(4, 'lnCHS_NKA_term', lnCHS_NKA_term)

    params_S=pd.DataFrame(data={'lnCHS_NKA_term': lnCHS_NKA_term,
    'NaKAl':NaKAl,
    'DeltaGRT': DeltaGRT,
    'SumXMAM': SumXMAM ,
    'lnCs': lnCs ,
    'lnXFeO': lnXFeO ,
    'LnrFeO': LnrFeO ,
    'lnaFeS': lnaFeS ,
    'C1PC2erf': C1PC2erf ,
    'lnS': lnS ,
    'S2_calc': S2_calc ,
    'lnXH2Ot': lnXH2Ot ,
    'KOH':KOH  ,
    'XOH': XOH ,
    'lnXOH': lnXOH ,
    'XH2Om': XH2Om ,
    'lnXH2Om':  lnXH2Om,
    'lnXOH_XH2O': lnXOH_XH2O ,
    'lnCHScalc':  lnCHScalc,
    'HScal': HScal ,
    'NaKAlterm':NaKAlterm  ,
    'HScal2':HScal2  })



    if any(Liqs.columns.str.contains('Fe3Fet_Liq')):
        if Fe3Fet_Liq is not None:
            print('replacing Fe3Fet_Liq in the original dataframe with that input into the function')
            Liqs['Fe3Fet_Liq']=Fe3Fet_Liq

    else:
        Liqs.insert(2, 'Fe3Fet_Liq', Fe3Fet_Liq)

    df_out=pd.concat([Liqs, params_S, mol_cat_norm], axis=1)
    df_out.insert(0, 'SCSS_Tot', Stot_calc_H2O)
    return df_out

## Blanchard et al. 2021 - SCSS calculations

def calculate_B2021_SCSS(*, df, T_K, P_kbar, Fe_FeNiCu_Sulf=None,  Cu_FeNiCu_Sulf=None, Ni_FeNiCu_Sulf=None, H2O_Liq=None,
Fe_Sulf=None, Ni_Sulf=None, Cu_Sulf=None,
Ni_Liq=None, Cu_Liq=None):
    '''
    Calculates SCSS using the model of Blanchard et al. 2021
    doi: https://doi.org/10.2138/am-2021-7649
    designed for calculating sulfide saturation in peridotitic melt at upper         mantle conditions


    Parameters
    -----------
    df: pandas.DataFrame
        Dataframe of liquid compositions from the import_data function.
        Needs to have the headings "SiO2_Liq", "TiO2_Liq" etc, with
        compositions in oxide wt%. all FeO should be entered as FeOt_Liq, which can then be partitioned
        using the Fe3Fet_Liq input. Heading order doesn't matter.

    T_K: int, float, pandas.Series
        Temperature in Kelvin.

    P_kbar: int, float, pandas.Series
        Pressure in kbar


    H2O_Liq: int, float, pandas series (optional)
        Overwrites H2O in the user entered dataframe

    Sulfide Composition: Options to calculate from the liquid composition, enter the comp in el wt%,
    or enter the FeFeNiCu, Cu

      if you want to calculate sulfide composition:

            Fe_FeNiCu_Sulf = "Calc_Smythe", also needs Ni_Sulf_init, and Cu_Sulf_init
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
    -----------
    pandas.DataFrame
        Calculated SCSS using eq11, eq12, mol fractions, and input compositions.



    '''
    # Coefficients for model 1 (Equation 11)
    a_m1=27
    B_m1=-4621
    C_m1=-193
    A_Si_m1=-25
    A_Ti_m1=-13
    A_Al_m1=-18
    A_Mg_m1=-16
    A_Fe_m1=-32
    A_Ca_m1=-14
    A_Na_m1=-17
    A_K_m1=-27
    A_H_m1=-19
    A_SiFe_m1=76

    # Coefficients for model 2 (Equation 12)



    a_m2=7.95
    B_m2=18159
    C_m2=-190
    A_Si_m2=-32677
    A_Ti_m2=-15014
    A_Al_m2=-23071
    A_Mg_m2=-18258
    A_Fe_m2=-41706
    A_Ca_m2=-14668
    A_Na_m2=-19529
    A_K_m2=-34641
    A_H_m2=-22677
    A_SiFe_m2=120662
    Liqs=df.copy()

    Liqs=calculate_sulfide_comp_generic(
    Fe_Sulf=Fe_Sulf, Ni_Sulf=Ni_Sulf, Cu_Sulf=Cu_Sulf, Fe_FeNiCu_Sulf=Fe_FeNiCu_Sulf,
    Cu_FeNiCu_Sulf=Cu_FeNiCu_Sulf, Ni_FeNiCu_Sulf=Ni_FeNiCu_Sulf,
    Ni_Liq=Ni_Liq, Cu_Liq=Cu_Liq, df_c=Liqs, T_K=T_K)


    calcs=calculate_hydrous_cat_fractions_liquid(liq_comps=df)

    P_GPa =P_kbar/10

    # Calculating  cation sum term
    catsum_term_m1=(+calcs['Si_Liq_cat_frac']*A_Si_m1 + calcs['Ti_Liq_cat_frac']*A_Ti_m1
    +calcs['Al_Liq_cat_frac']*A_Al_m1 + calcs['Mg_Liq_cat_frac']*A_Mg_m1
    +calcs['Fet_Liq_cat_frac']*A_Fe_m1 + calcs['Ca_Liq_cat_frac']*A_Ca_m1
    +calcs['Na_Liq_cat_frac']*A_Na_m1 + calcs['K_Liq_cat_frac']*A_K_m1
    +calcs['H2O_Liq_cat_frac']*A_H_m1
    +A_SiFe_m1*(calcs['Fet_Liq_cat_frac']*calcs['Si_Liq_cat_frac'] )
                    + calcs['K_Liq_cat_frac']*A_K_m1)

    catsum_term_m2=(+calcs['Si_Liq_cat_frac']*A_Si_m2 + calcs['Ti_Liq_cat_frac']*A_Ti_m2
    +calcs['Al_Liq_cat_frac']*A_Al_m2 + calcs['Mg_Liq_cat_frac']*A_Mg_m2
    +calcs['Fet_Liq_cat_frac']*A_Fe_m2 + calcs['Ca_Liq_cat_frac']*A_Ca_m2
    +calcs['Na_Liq_cat_frac']*A_Na_m2 + calcs['K_Liq_cat_frac']*A_K_m2
    +calcs['H2O_Liq_cat_frac']*A_H_m2
    +A_SiFe_m2*(calcs['Fet_Liq_cat_frac']*calcs['Si_Liq_cat_frac'] )
                    + calcs['K_Liq_cat_frac']*A_K_m2)




    lnSCSS_equation11=(a_m1+B_m1/T_K + (C_m1 * P_GPa)/T_K +
                    catsum_term_m1 +np.log(Liqs['Fe_FeNiCu_Sulf_calc'])
                    -np.log(calcs['Fet_Liq_cat_frac'])
                    )
    SCSS_eq11=np.exp(lnSCSS_equation11)
    lnSCSS_equation12=(a_m2+B_m2/T_K + (C_m2 * P_GPa)/T_K +
                    catsum_term_m2/T_K +np.log(Liqs['Fe_FeNiCu_Sulf_calc'])
                    -np.log(calcs['Fet_Liq_cat_frac'])
                    )
    SCSS_eq12=np.exp(lnSCSS_equation12)

    Liqs.insert(0, 'SCSS_eq11', SCSS_eq11)
    Liqs.insert(1, 'SCSS_eq12', SCSS_eq12)

    return Liqs

## Fortin et al. (2015) SCSS Calculation

def calculate_F2015_SCSS(df, T_K, P_kbar, H2O_Liq=None):
    '''
    Calculates SCSS using the model of Fortin et al. (2015).
    doi: http://dx.doi.org/10.1016/j.gca.2015.03.0220
    Unlike the models of Symthe et al. (2017) or O'Neill (2021) this model doesn't
    require the sulfide composition, or the proportion of Fe3 in the liqiud.

    Parameters
    -----------
    df: pandas.DataFrame
        Dataframe of liquid compositions from the import_data function.
        Needs to have the headings "SiO2_Liq", "TiO2_Liq" etc, with
        compositions in oxide wt%. all FeO should be entered as FeOt_Liq, which can then be partitioned
        using the Fe3Fet_Liq input. Heading order doesn't matter.

    T_K: int, float, pandas.Series
        Temperature in Kelvin.

    P_kbar: int, float, pandas.Series
        Pressure in kbar

    Returns
    -----------
    pandas.DataFrame
        Calculated SCSS, mol fractions, and input compositions.

    '''

    df_c=df.copy()
    if H2O_Liq is not None:
        df_c['H2O_Liq']=H2O_Liq
    mol_fracs=calculate_hydrous_mol_fractions_liquid(df_c)
    Ln_SCSS_Calc=(34.7837+(1/T_K*-5772.322)+
((P_kbar/10)/T_K*-346.5377)+
(mol_fracs['SiO2_Liq_mol_frac']*-25.4986)+
(mol_fracs['TiO2_Liq_mol_frac']*-18.3435)+
(mol_fracs['Al2O3_Liq_mol_frac']*-27.3807)+
(mol_fracs['FeOt_Liq_mol_frac']*-17.2752)+
(mol_fracs['MnO_Liq_mol_frac']*0)+
(mol_fracs['MgO_Liq_mol_frac']*-22.3975)+
(mol_fracs['CaO_Liq_mol_frac']*-20.3778)+
(mol_fracs['Na2O_Liq_mol_frac']*-18.9539)+
(mol_fracs['K2O_Liq_mol_frac']*-32.1944)+
(mol_fracs['P2O5_Liq_mol_frac']*0)+
(mol_fracs['H2O_Liq_mol_frac']*-20.3934))
    SCSS_Calc=np.exp(Ln_SCSS_Calc)
    out=pd.concat([mol_fracs, df_c], axis=1)
    out.insert(0, 'SCSS_ppm_Fortin2015', SCSS_Calc)
    return out
## Liu et al. 2021 SCSS calculations
def calculate_Liu2021_SCSS(df, T_K, P_kbar, Fe_FeNiCu_Sulf=None, H2O_Liq=None, Cu_FeNiCu_Sulf=None, Ni_FeNiCu_Sulf=None,
Ni_Liq=None, Cu_Liq=None, Fe_Sulf=None, Cu_Sulf=None, Ni_Sulf=None, Ni_Sulf_init=5, Cu_Sulf_init=5, Fe3Fet_Liq=None):
    '''
    Calculates the SCSS using the model of Liu et al. (2021), doesnt depend on silicate melt composition apart from H2O_Liq.
    Doi - https://doi.org/10.1016/j.chemgeo.2020.119913

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

            Fe_FeNiCu_Sulf = "Calc_Smythe", also needs Ni_Sulf_init, and Cu_Sulf_init
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



    '''
    df_c=df.copy()

    Fe_FeNiCu_Sulf_calc=calculate_sulfide_comp_generic(
    Fe_Sulf=Fe_Sulf, Ni_Sulf=Ni_Sulf, Cu_Sulf=Cu_Sulf, Fe_FeNiCu_Sulf=Fe_FeNiCu_Sulf,
    Cu_FeNiCu_Sulf=Cu_FeNiCu_Sulf, Ni_FeNiCu_Sulf=Ni_FeNiCu_Sulf,
    Ni_Liq=Ni_Liq, Cu_Liq=Cu_Liq, df_c=df_c)

    scss=Fe_FeNiCu_Sulf_calc*np.exp(13.88-9744/T_K-328*(P_kbar/10)/T_K)+104*H2O_Liq

    df_c.insert(0, 'SCSS_calc', scss)
    return df_c

## ONeill 2021 SCSS calculations

def calculate_ONeill_sulf(FeOt_Liq, Ni_Liq, Cu_Liq, MgO_Liq=None, Fe3Fet_Liq=None):
    '''
    Calculating the Fe_FeNiCu ratio in the sulfide using the empirical
    parameterizatin of ONeill et al. (2022)
    doi: https://doi.org/10.1002/9781119473206.ch10

    '''

    Fe_FeNiCu=1/(1+(Ni_Liq/(FeOt_Liq*(1-Fe3Fet_Liq)))*0.013+(Cu_Liq/(FeOt_Liq*(1-Fe3Fet_Liq)))*0.025)
    return Fe_FeNiCu
# New 2022 Spreadsheet with Mavrogenes that is circulating

def calculate_OM2022_SCSS(*, df, T_K, P_kbar, Fe3Fet_Liq=None, Fe_FeNiCu_Sulf=None,
Cu_FeNiCu_Sulf=None, Ni_FeNiCu_Sulf=None, Fe_Sulf=None, Ni_Sulf=None, Cu_Sulf=None,
Ni_Liq=None, Cu_Liq=None, Ni_Sulf_init=5, Cu_Sulf_init=5):

    '''
    Calculates SCSS using the model of O'Neill (2021) as implemented in the supporting spreadsheet of ONeill and Mavrogenes, 2022
    doi: https://doi.org/10.1016/j.gca.2022.06.020
    The only difference between OM2022 and O2021 is a 7.2*Fe*Si term in 2021, 7.2*(Mn+Fe)*Si in 2022.
    This function has options for users to
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

            Fe_FeNiCu_Sulf = "Calc_Smythe", also needs Ni_Sulf_init, and Cu_Sulf_init
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
    df_c['P2O5_Liq']=0 # Doesnt have P2O5 in the input
    liqs=calculate_anhydrous_cat_fractions_liquid(liq_comps=df_c)

    if isinstance(Fe3Fet_Liq, str) and Fe3Fet_Liq == "Calc_ONeill":
        Fe2O3_Calc=np.exp(1.46-0.177*df_c['MgO_Liq'])
        Fe3Fet_Liq=Fe2O3_Calc*0.8998/(df_c['FeOt_Liq'])
        df_c['Fe3Fet_Liq']=Fe3Fet_Liq
        liqs['Fe3Fet_Liq_calc']=Fe3Fet_Liq
    if Fe3Fet_Liq is not None and not isinstance(Fe3Fet_Liq, str):
          df_c['Fe3Fet_Liq']=Fe3Fet_Liq

    if Fe3Fet_Liq is None:
        Fe3Fet_Liq=0

    df_c2=calculate_sulfide_comp_generic(
    Fe_Sulf=Fe_Sulf, Ni_Sulf=Ni_Sulf, Cu_Sulf=Cu_Sulf, Fe_FeNiCu_Sulf=Fe_FeNiCu_Sulf,
    Cu_FeNiCu_Sulf=Cu_FeNiCu_Sulf, Ni_FeNiCu_Sulf=Ni_FeNiCu_Sulf,
    Ni_Liq=Ni_Liq, Cu_Liq=Cu_Liq, df_c=df_c, T_K=T_K)

    liqs=calculate_anhydrous_cat_fractions_liquid(liq_comps=df_c)

    liqs['Fe2_Liq_cat_frac']=liqs['Fet_Liq_cat_frac']*(1-Fe3Fet_Liq)





    liqs['LnCS2_calc']=(8.77-23590/T_K+(1673/T_K)*(6.7*(liqs['Na_Liq_cat_frac']+liqs['K_Liq_cat_frac'])
    +4.9*liqs['Mg_Liq_cat_frac']+8.1*liqs['Ca_Liq_cat_frac']+8.9*(liqs['Fet_Liq_cat_frac']+liqs['Mn_Liq_cat_frac'])
    +5*liqs['Ti_Liq_cat_frac']+1.8*liqs['Al_Liq_cat_frac']
    -22.2*liqs['Ti_Liq_cat_frac']*(liqs['Fet_Liq_cat_frac']+liqs['Mn_Liq_cat_frac'])
        +7.2*((liqs['Fet_Liq_cat_frac']+liqs['Mn_Liq_cat_frac'])*liqs['Si_Liq_cat_frac']))-2.06*erf(-7.2*(liqs['Fet_Liq_cat_frac']+liqs['Mn_Liq_cat_frac'])))

    liqs['DeltaG']=((137778-91.666*T_K+8.474*T_K*np.log(T_K))/(8.31441*T_K)+(-291*(P_kbar/10)+351*erf((P_kbar/10)))/T_K)

    # print('Cation fractions')
    # print()


    liqs['Ln_a_FeS']=(np.log(df_c2['Fe_FeNiCu_Sulf_calc']*(1-liqs['Fe2_Liq_cat_frac'])))

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

# Old 2021 spreadsheet that is circulating

def calculate_O2021_SCSS(*, df, T_K, P_kbar, Fe3Fet_Liq=None, Fe_FeNiCu_Sulf=None, Cu_FeNiCu_Sulf=None, Ni_FeNiCu_Sulf=None,
 Fe_Sulf=None, Ni_Sulf=None, Cu_Sulf=None,
Ni_Liq=None, Cu_Liq=None, Ni_Sulf_init=5, Cu_Sulf_init=5):

    '''
    Calculates SCSS using the model of O'Neill (2021).
    doi:  https://doi.org/10.1002/9781119473206.ch10
    Has  options for users to
    calculate sulfide composition from liquid composition using the approach of Smythe et al. (2017), the empirical parameterization of O'Neill (2021),  or input a sulfide composition.

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

            Fe_FeNiCu_Sulf = "Calc_Smythe", also needs Ni_Sulf_init, and Cu_Sulf_init
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
    df_c['P2O5_Liq']=0 # Doesnt have P2O5 in the input




    if isinstance(Fe3Fet_Liq, str) and Fe3Fet_Liq == "Calc_ONeill":
        Fe2O3_Calc=np.exp(1.46-0.177*df_c['MgO_Liq'])
        Fe3Fet_Liq=Fe2O3_Calc*0.8998/(df_c['FeOt_Liq'])
        df_c['Fe3Fet_Liq']=Fe3Fet_Liq
        df_c['Fe3Fet_Liq_calc']=Fe3Fet_Liq
    if Fe3Fet_Liq is not None and not isinstance(Fe3Fet_Liq, str):
          df_c['Fe3Fet_Liq']=Fe3Fet_Liq

    if Fe3Fet_Liq is None:
        Fe3Fet_Liq=0

    df_c2=calculate_sulfide_comp_generic(
    Fe_Sulf=Fe_Sulf, Ni_Sulf=Ni_Sulf, Cu_Sulf=Cu_Sulf, Fe_FeNiCu_Sulf=Fe_FeNiCu_Sulf,
    Cu_FeNiCu_Sulf=Cu_FeNiCu_Sulf, Ni_FeNiCu_Sulf=Ni_FeNiCu_Sulf,
    Ni_Liq=Ni_Liq, Cu_Liq=Cu_Liq, df_c=df_c, T_K=T_K)

    liqs=calculate_anhydrous_cat_fractions_liquid(liq_comps=df_c)

    liqs['Fe2_Liq_cat_frac']=liqs['Fet_Liq_cat_frac']*(1-Fe3Fet_Liq)




    df_c2['LnCS2_calc']=(8.77-23590/T_K+(1673/T_K)*(6.7*(liqs['Na_Liq_cat_frac']+liqs['K_Liq_cat_frac'])
    +4.9*liqs['Mg_Liq_cat_frac']+8.1*liqs['Ca_Liq_cat_frac']+8.9*(liqs['Fet_Liq_cat_frac']+liqs['Mn_Liq_cat_frac'])
    +5*liqs['Ti_Liq_cat_frac']+1.8*liqs['Al_Liq_cat_frac']
    -22.2*liqs['Ti_Liq_cat_frac']*(liqs['Fet_Liq_cat_frac']+liqs['Mn_Liq_cat_frac'])
        +7.2*(liqs['Fet_Liq_cat_frac']*liqs['Si_Liq_cat_frac']))-2.06*erf(-7.2*(liqs['Fet_Liq_cat_frac']+liqs['Mn_Liq_cat_frac'])))

    df_c2['DeltaG']=((137778-91.666*T_K+8.474*T_K*np.log(T_K))/(8.31441*T_K)+(-291*(P_kbar/10)+351*erf((P_kbar/10)))/T_K)

    # print('Cation fractions')
    # print()


    df_c2['Ln_a_FeS']=(np.log(df_c2['Fe_FeNiCu_Sulf_calc']*(1-liqs['Fe2_Liq_cat_frac'])))

    df_c2['Ln_a_FeO']=( np.log(liqs['Fe2_Liq_cat_frac'])+(((1-liqs['Fe2_Liq_cat_frac'])**2)*
            (28870-14710*liqs['Mg_Liq_cat_frac']+1960*liqs['Ca_Liq_cat_frac']+43300*liqs['Na_Liq_cat_frac']+95380*liqs['K_Liq_cat_frac']-76880*liqs['Ti_Liq_cat_frac'])
            +(1-liqs['Fe2_Liq_cat_frac'])*(-62190*liqs['Si_Liq_cat_frac']+31520*liqs['Si_Liq_cat_frac']*liqs['Si_Liq_cat_frac']))/(8.31441*T_K))

    df_c2['LnS']=(df_c2['LnCS2_calc']+df_c2['DeltaG']+df_c2['Ln_a_FeS']-df_c2['Ln_a_FeO'])

    df_c2['SCSS2_ppm']=np.exp(df_c2['LnS'])

    cols_to_move = ['SCSS2_ppm', 'LnS', "Ln_a_FeO",
                    'Ln_a_FeS', 'DeltaG', 'LnCS2_calc']

    Liqs = df_c2[cols_to_move +
                                    [col for col in df_c2.columns if col not in cols_to_move]]

    return Liqs



## Smythe SCSS calculations - More steps, as iteratively calculate sulfide composition.


def Loop_Smythe_sulf_calc_residual(single_argx0, FeO_Liq, T_K,  Ni_Liq, Cu_Liq):
    '''
    Calculates the residual between the Ni and Cu content of the Liquid and that of the sulfide for a
    Ni and Cu Kd calculated from Kiseeva et al. (2015).
    This function is used by the calculate_Smythe_sulf_minimisation function to calculate the sulfide compositoin,
    analagous to the Solver function in the supporting spreadsheet of Smythe et al.(2017).

    Parameters
    -----------
    single_argx0: float
        arguement for Ni and Cu in the sulfide to be solved by the scipy minimisation function function.

    FeO_Liq: int, float
        FeO (2+) content of the liquid in wt%

    Ni_Liq: int, float
        Ni content of the Liquid in ppm

    Cu_Liq: int, float
        Cu content of the Liquid in ppm

    T_K: int, float
        Temperature in Kelvin.


    Returns
    -----------
    float
        Calculated residual between measured Ni and Cu in the liquid, and that predicted
        from the sulfide composition and the Kd.
    '''
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

    else:

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


def calculate_sulf_kds(Ni_Sulf, Cu_Sulf, FeOt_Liq,  T_K, Fe3Fet_Liq=None):
    '''
    Calculates the Fe, O, S content of the sulfide using Kiseeva et al. (2015). Users can enter Ni, Cu and FeOt
    contents they have measured, or use this function once the
    Ni and Cu content of the sulfide have been calculated using scipy minimisation.
    Returns Kds, and S-O-Fe of sulfide.
    Also returns Se and Te calculated using Brenan, 2015

    Parameters
    -----------
    Ni_Sulf: pandas.Series
        Ni content of the sulfide in wt%, from the calculate_Smythe_sulf_minimisation function

    Cu_Sulf: pandas.Series
        Cu content of the sulfide in wt%, from the calculate_Smythe_sulf_minimisation function

    FeOt_Liq: pandas.Series
        FeOt content of the liquid in wt%, from the read_excel function.


    T_K: int, float, pandas.Series
        Temperature in Kelvin.

    P_kbar: int, float, pandas.Series
        Pressure in kbar

    Returns
    -----------
    pandas.DataFrame
        Calculated Cu, Ni, O, S and Fe content of sulfide.

    '''
    if Fe3Fet_Liq is not None:
        FeO_Liq=FeOt_Liq*(1-Fe3Fet_Liq)
    else:
        FeO_Liq=FeOt_Liq
    OCalc_CellAG12=0.24*FeO_Liq
    FeS_calc_AL19=((100-Ni_Sulf-Cu_Sulf-OCalc_CellAG12-Cu_Sulf*20.1442646/79.8557354
    -Ni_Sulf*35.32650016/64.67349984))*36.47119049/100
    FeCalc_CellAF12=FeS_calc_AL19*63.52880951/36.47119049

    AG6=Ni_Sulf/58.6934/(Ni_Sulf/58.6934+Cu_Sulf/63.546+FeCalc_CellAF12/55.845)
    AG8=AG6*0.97
    AG9=AG6*0.92

        # If FeO<13
    if isinstance(FeO_Liq, pd.Series):
        O_Sulf=(FeO_Liq*0.24*((1-AG8)**2))
        O_Sulf.loc[FeO_Liq>=13]=(FeO_Liq*0.24*((1-AG9)**2))
    if isinstance(FeO_Liq, float) or isinstance(FeO_Liq, int):
        if FeO_Liq>=13:
            O_Sulf=(FeO_Liq*0.24*((1-AG9)**2))
        else:
            O_Sulf=(FeO_Liq*0.24*((1-AG8)**2))


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
    eps_FeO_Cu=2.2030
    eps_Ni_Cu=0.00
    eps_Cu_Cu=-1.0310
    B_Cu=4200
    A_Cu=0.6995

    # Ni Coefficients
    eps_FeO_Ni=0
    eps_Ni_Ni=-0.904
    eps_Cu_Ni=0.00
    A_Ni=1.869
    B_Ni=3300

    # Pb Coefficients
    eps_FeO_Pb=0.454
    eps_Ni_Pb=0
    eps_Cu_Pb=0.00
    B_Pb=1260
    A_Pb=1.834


    # Ag Coefficients
    eps_FeO_Ag=1.901
    eps_Ni_Ag=0
    eps_Cu_Ag=0.00
    B_Ag=4300
    A_Ag=0.724

    # Zn Coefficients
    eps_FeO_Zn=-0.673
    eps_Ni_Zn=0.691
    eps_Cu_Zn=-1.023
    B_Zn=-990
    A_Zn=1.915

    # Cd Coefficients
    eps_FeO_Cd=0
    eps_Ni_Cd=0.48
    eps_Cu_Cd=-1.06
    B_Cd=1420
    A_Cd=1.919

    # Tl Coefficients
    eps_FeO_Tl=0.99
    eps_Ni_Tl=0.95
    eps_Cu_Tl=0.98
    B_Tl=0
    A_Tl=1.678

    # Mn Coefficients
    eps_FeO_Mn=-1.9
    eps_Ni_Mn=0.76
    eps_Cu_Mn=-0.41
    B_Mn=-1520
    A_Mn=1.629

    # In Coefficients
    eps_FeO_In=-1.61
    eps_Ni_In=-0.54
    eps_Cu_In=-2.36
    B_In=0
    A_In=2.502

    # Ti Coefficients
    eps_FeO_Ti=-12.91
    eps_Ni_Ti=-2.22
    eps_Cu_Ti=-2.04
    B_Ti=-2740
    A_Ti=1.027

    # Ga Coefficients
    eps_FeO_Ga=-5.09
    eps_Ni_Ga=0
    eps_Cu_Ga=-1.89
    B_Ga=-5470
    A_Ga=3.347


    # Sb Coefficients
    eps_FeO_Sb=-1.52
    eps_Ni_Sb=-2.32
    eps_Cu_Sb=-2.95
    B_Sb=-2670
    A_Sb=4.303

    # Co Coefficients
    eps_FeO_Co=0.60
    eps_Ni_Co=-0.28
    eps_Cu_Co=0.00
    B_Co=1280
    A_Co=1.964

    # V Cofficients
    eps_FeO_V=-5.34
    eps_Ni_V=0
    eps_Cu_V=0.00
    B_V=-2840
    A_V=2.524

    # Ge coefficients
    eps_FeO_Ge=-4.94
    eps_Ni_Ge=-1.65
    eps_Cu_Ge=-4.07
    B_Ge=-5000
    A_Ge=4.635

    # Cr coefficients
    eps_FeO_Cr=-0.76
    eps_Ni_Cr=0.44
    eps_Cu_Cr=0.00
    B_Cr=-1810
    A_Cr=2.356


    Fe_Corr=FeO_Liq/(Fe_Sulf/55.845/(Ni_Sulf/58.6934+Cu_Sulf/63.546+Fe_Sulf/55.845)) #AF7
    Ni_FeNiCu=Ni_Sulf/58.6934/(Ni_Sulf/58.6934+Cu_Sulf/63.546+Fe_Sulf/55.845)
    Cu_FeNiCu=Cu_Sulf/63.546/(Ni_Sulf/58.6934+Cu_Sulf/63.546+Fe_Sulf/55.845)
    Fe_FeNiCu=Fe_Sulf/55.845/(Ni_Sulf/58.6934+Cu_Sulf/63.546+Fe_Sulf/55.845)

    DNi=(10**(B_Ni/T_K+A_Ni-np.log10(Fe_Corr)+(1673*eps_FeO_Ni*np.log10(1-0.049*O_Sulf)/T_K)
            +(1673*eps_Ni_Ni*np.log10(1-Ni_FeNiCu)/T_K)+(1673*eps_Cu_Ni*np.log10(1-Cu_FeNiCu)/T_K)))

    DCu=(10**(B_Cu/T_K+A_Cu-0.5*np.log10(Fe_Corr)+(1673*eps_FeO_Cu*np.log10(1-0.049*O_Sulf)/T_K)
             +(1673*eps_Cu_Cu*np.log10(1-Cu_FeNiCu)/T_K)+(1673*eps_Ni_Cu*np.log10(1-Ni_FeNiCu)/T_K)))

    DPb=(10**(B_Pb/T_K+A_Pb-np.log10(Fe_Corr)
    +(1673*eps_FeO_Pb*np.log10(1-0.049*O_Sulf)/T_K)
    +(1673*eps_Ni_Pb*np.log10(1-Ni_FeNiCu)/T_K)
    +(1673*eps_Cu_Pb*np.log10(1-Cu_FeNiCu)/T_K)))

    DAg=(10**(B_Ag/T_K+A_Ag-0.5*np.log10(Fe_Corr)
    +(1673*eps_FeO_Ag*np.log10(1-0.049*O_Sulf)/T_K)
    +(1673*eps_Ni_Ag*np.log10(1-Ni_FeNiCu)/T_K)
    +(1673*eps_Cu_Ag*np.log10(1-Cu_FeNiCu)/T_K)))


    DZn=(10**(B_Zn/T_K+A_Zn-np.log10(Fe_Corr)
    +(1673*eps_FeO_Zn*np.log10(1-0.049*O_Sulf)/T_K)
     +(1673*eps_Ni_Zn*np.log10(1-Ni_FeNiCu)/T_K)
     +(1673*eps_Cu_Zn*np.log10(1-Cu_FeNiCu)/T_K)))

    DCd=(10**(B_Cd/T_K+A_Cd-np.log10(Fe_Corr)
    +(1673*eps_FeO_Cd*np.log10(1-0.049*O_Sulf)/T_K)
    +(1673*eps_Ni_Cd*np.log10(1-Ni_FeNiCu)/T_K)
    +(1673*eps_Cu_Cd*np.log10(1-Cu_FeNiCu)/T_K)))

    DTl=(10**(B_Tl/T_K+A_Tl-0.5*np.log10(Fe_Corr)
    +(1673*eps_FeO_Tl*np.log10(1-0.049*O_Sulf)/T_K)
     +(1673*eps_Ni_Tl*np.log10(1-Ni_FeNiCu)/T_K)
     +(1673*eps_Cu_Tl*np.log10(1-Cu_FeNiCu)/T_K)))

    DMn=(10**(B_Mn/T_K+A_Mn-np.log10(Fe_Corr)
    +(1673*eps_FeO_Mn*np.log10(1-0.049*O_Sulf)/T_K)
    +(1673*eps_Ni_Mn*np.log10(1-Ni_FeNiCu)/T_K)
    +(1673*eps_Cu_Mn*np.log10(1-Cu_FeNiCu)/T_K)))

    DIn=(10**(B_In/T_K+A_In-1.5*np.log10(Fe_Corr)
    +(1673*eps_FeO_In*np.log10(1-0.049*O_Sulf)/T_K)
     +(1673*eps_Ni_In*np.log10(1-Ni_FeNiCu)/T_K)
     +(1673*eps_Cu_In*np.log10(1-Cu_FeNiCu)/T_K)))

    DTi=(10**(B_Ti/T_K+A_Ti-2*np.log10(Fe_Corr)
    +(1673*eps_FeO_Ti*np.log10(1-0.049*O_Sulf)/T_K)
     +(1673*eps_Ni_Ti*np.log10(1-Ni_FeNiCu)/T_K)
     +(1673*eps_Cu_Ti*np.log10(1-Cu_FeNiCu)/T_K)))

    DGa=(10**(B_Ga/T_K+A_Ga-1.5*np.log10(Fe_Corr)
    +(1673*eps_FeO_Ga*np.log10(1-0.049*O_Sulf)/T_K)
     +(1673*eps_Ni_Ga*np.log10(1-Ni_FeNiCu)/T_K)
     +(1673*eps_Cu_Ga*np.log10(1-Cu_FeNiCu)/T_K)))

    DSb=(10**(B_Sb/T_K+A_Sb-1.5*np.log10(Fe_Corr)
    +(1673*eps_FeO_Sb*np.log10(1-0.049*O_Sulf)/T_K)
     +(1673*eps_Ni_Sb*np.log10(1-Ni_FeNiCu)/T_K)
     +(1673*eps_Cu_Sb*np.log10(1-Cu_FeNiCu)/T_K)))

    DCo=(10**(B_Co/T_K+A_Co-np.log10(Fe_Corr)
    +(1673*eps_FeO_Co*np.log10(1-0.049*O_Sulf)/T_K)
    +(1673*eps_Ni_Co*np.log10(1-Ni_FeNiCu)/T_K)
    +(1673*eps_Cu_Co*np.log10(1-Cu_FeNiCu)/T_K)))

    DV=(10**(B_V/T_K+A_V-1.5*np.log10(Fe_Corr)
    +(1673*eps_FeO_V*np.log10(1-0.049*O_Sulf)/T_K)
     +(1673*eps_Ni_V*np.log10(1-Ni_FeNiCu)/T_K)
     +(1673*eps_Cu_V*np.log10(1-Cu_FeNiCu)/T_K)))

    DGe=(10**(B_Ge/T_K+A_Ge-2*np.log10(Fe_Corr)
    +(1673*eps_FeO_Ge*np.log10(1-0.049*O_Sulf)/T_K)
     +(1673*eps_Ni_Ge*np.log10(1-Ni_FeNiCu)/T_K)
     +(1673*eps_Cu_Ge*np.log10(1-Cu_FeNiCu)/T_K)))

    DCr=(10**(B_Cr/T_K+A_Cr-np.log10(Fe_Corr)
    +(1673*eps_FeO_Cr*np.log10(1-0.049*O_Sulf)/T_K)
    +(1673*eps_Ni_Cr*np.log10(1-Ni_FeNiCu)/T_K)
    +(1673*eps_Cu_Cr*np.log10(1-Cu_FeNiCu)/T_K)))

    # Te and Se calculations using Brenan (2015) assuming only FeO
    DSe=(10**(3.47-3.07*10**(-2)*FeO_Liq - 9.13*10**(-4)*FeO_Liq**2
    ))
    DTe=(10**(4.3 -2.58*10**(-3)*FeO_Liq**2)
    )

    Ni_Melt_calc=Ni_Sulf*10000/DNi
    Cu_Melt_calc=Cu_Sulf*10000/DCu

    if isinstance(DNi, pd.Series):
        df_out=pd.DataFrame(data={
    'S_Sulf': S_Sulf, 'O_Sulf': O_Sulf, 'Fe_Sulf': Fe_Sulf, 'Ni_Sulf': Ni_Sulf, 'Cu_Sulf': Cu_Sulf, 'DNi': DNi,
    'DCu': DCu, 'DAg': DAg,  'DPb': DPb, 'DZn': DZn, 'DCd': DCd, 'DTl': DTl, 'DMn': DMn, 'DIn': DIn, 'DTi': DTi,
    'DGa': DGa, 'DSb': DSb, 'DCo': DCo, 'DV': DV, 'DGe': DGe, 'DCr': DCr, 'DSe_B2015': DSe,  'DTe_B2015': DTe})
    else:
        df_out=pd.DataFrame(data={
    'S_Sulf': S_Sulf, 'O_Sulf': O_Sulf, 'Fe_Sulf': Fe_Sulf, 'Ni_Sulf': Ni_Sulf, 'Cu_Sulf': Cu_Sulf, 'DNi': DNi,
    'DCu': DCu,  'DAg': DAg,  'DPb': DPb, 'DZn': DZn, 'DCd': DCd, 'DTl': DTl, 'DMn': DMn, 'DIn': DIn, 'DTi': DTi,
    'DGa': DGa, 'DSb': DSb, 'DCo': DCo, 'DV': DV, 'DGe': DGe , 'DCr': DCr, 'DSe_B2015': DSe,  'DTe_B2015': DTe}, index=[0])

    return df_out





## Calculating sulfide composition once you have iterated towards a respond.


def calculate_Symthe_sulf_minimisation(FeOt_Liq, Fe3Fet_Liq, T_K, Ni_Liq, Cu_Liq, Cu_Sulf_init=5, Ni_Sulf_init=5):

    '''
    Uses the Scipy optimize minimise function to find the Cu and Ni content of the sulfide
    that minimises the residual between liquid compositions, sulfide compositions and Kd,
    this minimis the Solver function in the supporting spreadsheet of Smythe et al.(2017).

    Parameters
    -----------
    FeOt_Liq: int, float, pd.Series
        FeOt content of the liquid in wt%

    Fe3Fet_Liq,: int, float, pd.Series
        Proportion of Fe3 in the liquid.

    Cu_Liq: int, float, pd.Series
        Cu content of the Liquid in ppm

    Ni_Liq: int, float, pd.Series
        Ni content of the Liquid in ppm

    T_K: int, float
        Temperature in Kelvin.


    Returns
    -----------
    float
        Calculated residual between measured Ni and Cu in the liquid, and that predicted
        from the sulfide composition and the Kd.
    '''

    Ni_Sulf=np.empty(len(FeOt_Liq), dtype=float)
    Cu_Sulf=np.empty(len(FeOt_Liq), dtype=float)


    bnds=((0, 30), (0, 30))

    for i in range(0, len(FeOt_Liq)):
        Calc_Sulf=scipy.optimize.minimize(Loop_Smythe_sulf_calc_residual, x0=(Ni_Sulf_init, Cu_Sulf_init), bounds=bnds,
                        args=(FeOt_Liq.iloc[i]*(1-Fe3Fet_Liq.iloc[i]),
                              T_K.iloc[i], Ni_Liq.iloc[i], Cu_Liq.iloc[i])).get('x')



        Ni_Sulf[i]=Calc_Sulf[0]
        Cu_Sulf[i]=Calc_Sulf[1]
    if len(Ni_Sulf)>1:
        df_out=pd.DataFrame(data={'Ni_Sulf': Ni_Sulf, 'Cu_Sulf': Cu_Sulf})
    if len(Ni_Sulf)==1:
        df_out=pd.DataFrame(data={'Ni_Sulf': Ni_Sulf, 'Cu_Sulf': Cu_Sulf}, index=[0])

    return df_out


def calculate_sulf_FeFeNiCu(Ni_Sulf, Cu_Sulf, Fe_Sulf):
    '''
    Calculates the atomic Fe/(Fe+Ni+Cu) ratio of the sulfide.

    Parameters
    -----------
    Ni_Sulf: int, float, pd.Series
        Ni content of sulfide in wt%

    Cu_Sulf: int, float, pd.Series
        Cu content of sulfide in wt%

    Fe_Sulf: int, float, pd.Series
        Cu content of sulfide in wt%

    Returns
    -----------
    pd.Series, float
       Fe/(Fe+Ni+Cu) ratio of the sulfide


    '''

    Ni_moles=Ni_Sulf/58.6934
    Cu_moles=Cu_Sulf/63.546
    Fe_moles=Fe_Sulf/55.8457
    FeFeNiCu=Fe_moles/(Fe_moles+Cu_moles+Ni_moles)
    return FeFeNiCu

def calculate_sulf_CuFeNiCu(Ni_Sulf, Cu_Sulf, Fe_Sulf):
    '''
    Calculates the atomic Cu/(Fe+Ni+Cu) ratio of the sulfide.

    Parameters
    -----------
    Ni_Sulf: int, float, pd.Series
        Ni content of sulfide in wt%

    Cu_Sulf: int, float, pd.Series
        Cu content of sulfide in wt%

    Fe_Sulf: int, float, pd.Series
        Cu content of sulfide in wt%

    Returns
    -----------
    pd.Series, float
       Cu/(Fe+Ni+Cu) ratio of the sulfide


    '''
    Ni_moles=Ni_Sulf/58.6934
    Cu_moles=Cu_Sulf/63.546
    Fe_moles=Fe_Sulf/55.8457
    CuFeNiCu=Cu_moles/(Fe_moles+Cu_moles+Ni_moles)
    return CuFeNiCu

def calculate_sulf_NiFeNiCu(Ni_Sulf, Cu_Sulf, Fe_Sulf):
    '''
    Calculates the atomic Ni/(Fe+Ni+Cu) ratio of the sulfide.

    Parameters
    -----------
    Ni_Sulf: int, float, pd.Series
        Ni content of sulfide in wt%

    Cu_Sulf: int, float, pd.Series
        Cu content of sulfide in wt%

    Fe_Sulf: int, float, pd.Series
        Cu content of sulfide in wt%

    Returns
    -----------
    pd.Series, float
       Ni/(Fe+Ni+Cu) ratio of the sulfide


    '''
    Ni_moles=Ni_Sulf/58.6934
    Cu_moles=Cu_Sulf/63.546
    Fe_moles=Fe_Sulf/55.8457
    NiFeNiCu=Ni_moles/(Fe_moles+Cu_moles+Ni_moles)
    return NiFeNiCu
## Smythe Parameterization

def calculate_S2017_SCSS(*, df, T_K, P_kbar, Fe3Fet_Liq=None, Fe_FeNiCu_Sulf=None, Cu_FeNiCu_Sulf=None, Ni_FeNiCu_Sulf=None,
Fe_Sulf=None, Ni_Sulf=None, Cu_Sulf=None, Ni_Liq=None, Cu_Liq=None,
Ni_Sulf_init=5, Cu_Sulf_init=5):
    '''
    Calculates SCSS using the model of Smythe et al. (2017).
    doi: https://doi.org/10.2138/am-2017-5800CCBY

    Has with options for users to calculate sulfide composition from liquid composition, or input sulfide  composition.

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

            Fe_FeNiCu_Sulf = "Calc_Smythe", also needs Ni_Sulf_init, and Cu_Sulf_init
            Calculates sulfide composition analagous to the Solver function in the Smythe et al. (2017) spreadsheet.
            Here, we use Scipy optimize to find the ideal Ni and Cu
            contents using Kds from Kiseeva et al. (2015) and the Ni and Cu content of the melt.
            Requires user to also enter Ni_Liq and Cu_Liq.

            Also can enter Ni_Sulf_init, and Cu_Sulf_init to help convergence (Default 5 for both)

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
    if Fe3Fet_Liq is not None:
        df_c['Fe3Fet_Liq']=Fe3Fet_Liq
    # First, calculate silicate hydrous mole fractions, as true regardless of choice of sulfide composition
    Smythe_calcs=calculate_Smythe_silicate_mole_fractions(df_c, Fe3Fet_Liq)


    df_c=calculate_sulfide_comp_generic(
    Fe_Sulf=Fe_Sulf, Ni_Sulf=Ni_Sulf, Cu_Sulf=Cu_Sulf, Fe_FeNiCu_Sulf=Fe_FeNiCu_Sulf,
    Cu_FeNiCu_Sulf=Cu_FeNiCu_Sulf, Ni_FeNiCu_Sulf=Ni_FeNiCu_Sulf,
    Ni_Liq=Ni_Liq, Cu_Liq=Cu_Liq, df_c=df_c, T_K=T_K)


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
+np.log(df_c['Fe_FeNiCu_Sulf_calc'])-np.log(Smythe_calcs['Fe2_wt_atom'])-269.4*0.1*P_kbar/T_K


    )


    Smythe_calcs['SCSS_ideal_ppm_Smythe2017']=np.exp(Smythe_calcs['log_SCSS_ideal'])
    Smythe_calcs['SCSS_ideal_ppm_Smythe2017_1sigma']=Smythe_calcs['SCSS_ideal_ppm_Smythe2017']*0.273169775211857
    Smythe_calcs['T_Input_K']=T_K
    Smythe_calcs['P_Input_kbar']=P_kbar
    Smythe_calcs['Fe_FeNiCu_Sulf']=df_c['Fe_FeNiCu_Sulf_calc']
    Smythe_calcs['Fe3Fet_Liq_input']=Fe3Fet_Liq

    # if (isinstance(Fe_FeNiCu_Sulf, float) or isinstance(Fe_FeNiCu_Sulf, int)) or (isinstance(Fe_FeNiCu_Sulf, pd.Series)) and (Cu_FeNiCu_Sulf is None and Ni_FeNiCu_Sulf is None):
    #     non_ideal=False
    #     print('no non ideal SCSS as no Cu/CuFeNiCu')
    # elif not (isinstance(Fe_FeNiCu_Sulf, pd.Series)):
    #     if Fe_FeNiCu_Sulf == 'Calc_ONeill':
    #         non_ideal=False
    #     if Fe_FeNiCu_Sulf == 'Calc_Smythe':
    #         non_ideal=True
    #
    #
    # elif (isinstance(Fe_FeNiCu_Sulf, pd.Series)) and (Cu_FeNiCu_Sulf is not None and Ni_FeNiCu_Sulf is not None):
    #     non_ideal=True
    #
    #
    # if
    if 'Ni_FeNiCu_Sulf_calc' in df_c.columns and 'Cu_FeNiCu_Sulf_calc' in df_c.columns:
        Smythe_calcs['log_SCSS_non_ideal']=(
        (122175-80.28*T_K+8.474*T_K*np.log(T_K))/(8.314*T_K)+9.352+(Smythe_calcs['Si_XA_non_ideal']+Smythe_calcs['Ti_XA_non_ideal']
    +Smythe_calcs['Al_XA_non_ideal']+Smythe_calcs['Mg_XA_non_ideal']+Smythe_calcs['Fe2_XA_non_ideal']+Smythe_calcs['Ca_XA_non_ideal']
    +Smythe_calcs['Na_XA_non_ideal']+Smythe_calcs['K_XA_non_ideal']+Smythe_calcs['H_XA_non_ideal']+Smythe_calcs['Si*Fe_non_ideal'])
    /T_K+np.log(df_c['Fe_FeNiCu_Sulf_calc'])-np.log(Smythe_calcs['Fe2_wt_atom'])-264.85*0.1*P_kbar/T_K+546.362*((df_c['Cu_FeNiCu_Sulf_calc']**2 +df_c['Cu_FeNiCu_Sulf_calc']*df_c['Ni_FeNiCu_Sulf_calc'])/T_K)
        )
        Smythe_calcs['SCSS_non_ideal_ppm_Smythe2017']=np.exp(Smythe_calcs['log_SCSS_non_ideal'])
        Smythe_calcs['SCSS_non_ideal_ppm_Smythe2017_1sigma']=Smythe_calcs['SCSS_non_ideal_ppm_Smythe2017']*0.267299081373473

        cols_to_move = ['SCSS_ideal_ppm_Smythe2017', 'SCSS_ideal_ppm_Smythe2017_1sigma',
        'SCSS_non_ideal_ppm_Smythe2017', 'SCSS_non_ideal_ppm_Smythe2017_1sigma',
        'T_Input_K', "P_Input_kbar",'Fe_FeNiCu_Sulf', 'Fe3Fet_Liq_input']

    else:
        print('no non ideal SCSS as no Cu/CuFeNiCu')


        cols_to_move = ['SCSS_ideal_ppm_Smythe2017', 'SCSS_ideal_ppm_Smythe2017_1sigma',
        'T_Input_K', "P_Input_kbar",'Fe_FeNiCu_Sulf', 'Fe3Fet_Liq_input']

    Smythe_calcs = Smythe_calcs[cols_to_move +
                                    [col for col in Smythe_calcs.columns if col not in cols_to_move]]

    Smythe_calcs['Fe_FeNiCu_Sulf_calc']=df_c['Fe_FeNiCu_Sulf_calc']
    Concat=pd.concat([Smythe_calcs, df_c], axis=1)


    return Concat


## Generic function for calculating sulfide composition

def calculate_sulfide_comp_generic(*, T_K,
Fe_Sulf=None, Ni_Sulf=None, Cu_Sulf=None, Fe_FeNiCu_Sulf=None, Cu_FeNiCu_Sulf=None, Ni_FeNiCu_Sulf=None,
Ni_Liq=None, Cu_Liq=None, df_c=None, Ni_Sulf_init=5, Cu_Sulf_init=5 ):
    """ Generic function for calculating sulfide composition,
    either from input Fe-Ni-Cu contents in wt%, or using models of
    ONeill and Smythe

    Sulfide Composition: Options to calculate from the liquid composition, enter the comp in el wt%,
    or enter the FeFeNiCu, Cu

      if you want to calculate sulfide composition:

            Fe_FeNiCu_Sulf = "Calc_Smythe", also needs Ni_Sulf_init, and Cu_Sulf_init
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


    """
    # First, if the user entered an integer or float for Ni and Cu, turn into a panda series
    if isinstance(Ni_Liq, int) is True or isinstance(Ni_Liq, float) is True:
        Ni_Liq=pd.Series(Ni_Liq, index=range(len(df_c)))
    if isinstance(Cu_Liq, int) is True or isinstance(Cu_Liq, float) is True:
        Cu_Liq=pd.Series(Cu_Liq, index=range(len(df_c)))

    # Checking no contradictory intputs
    if Fe_Sulf is not None and Ni_Sulf is not None and Cu_Sulf is not None and Fe_FeNiCu_Sulf is not None:
        raise ValueError('You have entered both a Fe_FeNiCu_Sulf ratio, and the conc of Fe, Ni and Cu in your sulf. Please enter one set of inputs or another')

    # IF they have entered a string,
    elif isinstance(Fe_FeNiCu_Sulf, str) and (Fe_FeNiCu_Sulf=="Calc_Smythe" or Fe_FeNiCu_Sulf == "Calc_ONeill"):
        if Ni_Liq is None:
            raise ValueError('If you select Fe_FeNiCu_Sulf = model, you need to enter the concentration of Cu and Ni in the liquid in ppm"')
        if Cu_Liq is None:
            raise ValueError('If you select Fe_FeNiCu_Sulf = model, you need to enter the concentration of Cu and Ni in the liquid in ppm"')


    elif Fe_FeNiCu_Sulf is not None:
        print('Using inputted Fe_FeNiCu_Sulf ratio for calculations.')

    elif Fe_Sulf is not None and Ni_Sulf is not None and Cu_Sulf is not None:
        print('Using inputted Sulf compositions to calculate Fe_FeNiCu_Sulf ratios.')
        if isinstance(Fe_Sulf, int) is True:
            Fe_Sulf=float(Fe_Sulf)
        if isinstance(Ni_Sulf, int) is True:
            Ni_Sulf=float(Ni_Sulf)
        if isinstance(Cu_Sulf, int) is True:
            Cu_Sulf=float(Cu_Sulf)
    else:
        raise ValueError('Input for sulfide composition not recognised.')

    # Start with these none
    Cu_FeNiCu_Sulf_calc=None
    Ni_FeNiCu_Sulf_calc=None

    # If its a model, do this.
    if isinstance(Fe_FeNiCu_Sulf, str):
        if Fe_FeNiCu_Sulf=="Calc_Smythe":


            # This does the Scipy minimisation of Cu and Ni contents using Kiseeva et al. (2015)
            calc_sulf=calculate_Symthe_sulf_minimisation(FeOt_Liq=df_c['FeOt_Liq'], Fe3Fet_Liq=df_c['Fe3Fet_Liq'],
                                            T_K=T_K, Ni_Liq=Ni_Liq, Cu_Liq=Cu_Liq, Ni_Sulf_init=Ni_Sulf_init, Cu_Sulf_init=Cu_Sulf_init)

            # This feeds those result back into a simpler function to get the Fe, S and O content of the sulfide
            Sulf_All=calculate_sulf_kds(Ni_Sulf=calc_sulf['Ni_Sulf'],Cu_Sulf=calc_sulf['Cu_Sulf'],
                                FeOt_Liq=df_c['FeOt_Liq'], Fe3Fet_Liq=df_c['Fe3Fet_Liq'],
                                            T_K=T_K)

            Fe_FeNiCu_Sulf_calc=calculate_sulf_FeFeNiCu(Sulf_All['Ni_Sulf'], Sulf_All['Cu_Sulf'], Sulf_All['Fe_Sulf'])
            Cu_FeNiCu_Sulf_calc=calculate_sulf_CuFeNiCu(Sulf_All['Ni_Sulf'], Sulf_All['Cu_Sulf'], Sulf_All['Fe_Sulf'])
            Ni_FeNiCu_Sulf_calc=calculate_sulf_NiFeNiCu(Sulf_All['Ni_Sulf'], Sulf_All['Cu_Sulf'], Sulf_All['Fe_Sulf'])

            df_c['Ni_Sulf_Calc']=Sulf_All['Ni_Sulf']
            df_c['Cu_Sulf_Calc']=Sulf_All['Cu_Sulf']
            df_c['Fe_Sulf_Calc']=Sulf_All['Fe_Sulf']
            df_c['O_Sulf_Calc']=Sulf_All['O_Sulf']
            df_c['S_Sulf_Calc']=Sulf_All['S_Sulf']

        if Fe_FeNiCu_Sulf=="Calc_ONeill":
            Fe_FeNiCu_Sulf_calc=calculate_ONeill_sulf(FeOt_Liq=df_c['FeOt_Liq'],
            Ni_Liq=Ni_Liq, Cu_Liq=Cu_Liq, Fe3Fet_Liq=df_c['Fe3Fet_Liq'])
            Fe_FeNiCu_Sulf=Fe_FeNiCu_Sulf_calc

    # If they have given the sulfide chemistry
    elif Fe_Sulf is not None and Ni_Sulf is not None and Cu_Sulf is not None:
            Fe_FeNiCu_Sulf_calc=calculate_sulf_FeFeNiCu(Ni_Sulf, Cu_Sulf, Fe_Sulf)
            Cu_FeNiCu_Sulf_calc=calculate_sulf_CuFeNiCu(Ni_Sulf, Cu_Sulf, Fe_Sulf)
            Ni_FeNiCu_Sulf_calc=calculate_sulf_NiFeNiCu(Ni_Sulf, Cu_Sulf, Fe_Sulf)

    # If they have given the FeFeNiCu ratio, but nothing else
    else:
        Fe_FeNiCu_Sulf_calc=Fe_FeNiCu_Sulf
        if Cu_FeNiCu_Sulf is not None:
            Cu_FeNiCu_Sulf_calc=Cu_FeNiCu_Sulf
        if Ni_FeNiCu_Sulf is not None:
            Ni_FeNiCu_Sulf_calc=Ni_FeNiCu_Sulf

    # if isinstance(Fe_FeNiCu_Sulf_calc, float) is True:
    #     Fe_FeNiCu_Sulf_calc=Fe_FeNiCu_Sulf_calc
    # elif isinstance(Fe_FeNiCu_Sulf_calc, int) is True:
    #     Fe_FeNiCu_Sulf_calc=float(Fe_FeNiCu_Sulf_calc)
    # else:
    #     Fe_FeNiCu_Sulf_calc=Fe_FeNiCu_Sulf_calc.astype(float)

    df_c['Fe_FeNiCu_Sulf_calc']=Fe_FeNiCu_Sulf_calc
    if Cu_FeNiCu_Sulf_calc is not None:
        df_c['Cu_FeNiCu_Sulf_calc']=Cu_FeNiCu_Sulf_calc
    if Ni_FeNiCu_Sulf_calc is not None:
        df_c['Ni_FeNiCu_Sulf_calc']=Ni_FeNiCu_Sulf_calc

    return df_c











