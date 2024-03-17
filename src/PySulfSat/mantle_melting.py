import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings as w

## Mantle melting, adapted from Lee et al. (2012) and Wieser et al. (2020)

def Lee_Wieser_sulfide_melting(*,  M_Max=0.01, Modes, KDs,
                        N=3001, F=None,
                        S_Sulf=363636, elem_Per=30,
                        S_Mantle=200,
                        S_Melt_SCSS_2_ppm=980.7872088,
                        Prop_S6=0
                        ):

    """ Function for calculating trajectory of S during mantle melting,
    adapted from the Excel workbook of Lee et al. (2012) and the Python code
    of Wieser et al. (2020)

    Parameters
    -----------
    N: int
        Number of iterative steps for melting.

    Either

    M_Max: float
        Amount of peridotite left at final step. E.g. M_Max=0.7
        calculates 0-30% melting (F=1-M)

    F: np.array
        F you want

    Modes: Pandas.Dataframe
        Dataframe of ol-opx-cpx-sp modes, with column headings 'ol',
        'opx', 'cpx', 'sp'
        (can be 1 row or multiple for each step of melting, in which case N is adjusted to be length of this dataframe)

    KDs: pandas.DataFrame
        Dataframe of KD values. Can be 1 row or multiple for each step of melting

    S_Sulf: float, int
        S content of mantle sulfides in ppm

    S_Mantle: float, int
        S content of mantle in ppm

    elem_Per: float, int
        Content of element of interest in mantle prior to melting (in ppm).

    S_Melt_SCSS_2_ppm: float, int, pd.Series
        S2- content of mantle in ppm

    Prop_S6: float, int, pd.Series
        Proportion of S6, increases S_Melt_SCSS if S6+ is present


    Returns
    -------------
    pd.DataFrame: dataframe with melt extents, instantaneous and aggregate melts, and residual mantle source composition.


    """
    elx=KDs['element'][0]

    # Mass of peridotite




    if F is not None:
        M=1-F
    else:
        if len(Modes)==1:
            M=np.linspace(1, M_Max, N)
        else:
            if len(Modes)!=N:
                w.warn('You have inputted a dataframe of silicate modes that doesnt match the number of steps you asked for. We are changing the number of steps to match the length of your dataframe')
                M=M=np.linspace(1, M_Max, len(Modes))
            else:
                M=M=np.linspace(1, M_Max, len(Modes))




    if isinstance(S_Melt_SCSS_2_ppm, float) or isinstance(S_Melt_SCSS_2_ppm, int):
        S_Melt_SCSS_2_ppm=np.repeat(S_Melt_SCSS_2_ppm, len(M))

    S_Melt = S_Melt_SCSS_2_ppm/(1-Prop_S6) # Was wrong in previous version.

    # Setting up variables to be filled in loop
    #These ones have only 1 dimension for different M steps, as they do not vary with different mantle S contents
    DeltaM= np.empty(len(M), dtype=float)
    X_FStep=np.empty(len(M), dtype=float)
    X_FStart=np.empty( len(M), dtype=float)

    len_S_Mantle=1

# These ones have 2 dimensions, one for different M steps, one for different mantle S contents
    S_residue= np.empty((len_S_Mantle, len(M)), dtype=float)
    XSulf= np.empty((len_S_Mantle, len(M)), dtype=float)
    elem_residue= np.empty((len_S_Mantle, len(M)), dtype=float)
    elem_Melt_Inst= np.empty((len_S_Mantle, len(M)), dtype=float)
    elem_Melt_Agg= np.empty((len_S_Mantle, len(M)), dtype=float)
    S_Melt_Agg= np.empty((len_S_Mantle, len(M)), dtype=float)
    S_Melt_Inst= np.empty((len_S_Mantle, len(M)), dtype=float)
    DeltaXSulf=np.empty((len_S_Mantle, len(M)), dtype=float)
    DeltaXSil=np.empty((len_S_Mantle, len(M)), dtype=float)
    P=np.empty((len_S_Mantle, len(M)), dtype=float)
    Kd_el=np.empty((len_S_Mantle, len(M)), dtype=float)
    DeltaXMeltTot=np.empty((len_S_Mantle, len(M)), dtype=float)


    # Loop for element
    for i in range(0, len(M)):

        # If length of KDs>1, need to extract right one for each loop
        #print(KDs)
        if len(KDs)>1:
            KDs_trim=KDs.iloc[i].to_frame().T.reset_index(drop=True)
        else:
            KDs_trim=KDs
        KdelSulf=KDs_trim['sulf']

        # If length of Modes is >1, need to extract the right one for each loop
        if len(Modes)>1:
            KdelSil=Modes.iloc[i].mul(KDs_trim).sum(axis=1).iloc[0]
        else:
            KdelSil=Modes.mul(KDs_trim).sum(axis=1)







        if i==0: # Setting up the system before any melting takes place
            S_residue[:,i]=S_Mantle # Initial S contents of mantle before any melting set at those defined in S_Mantle
            XSulf[:,i]=S_residue[:,i]/S_Sulf # Proportion (by mass) of sulfides in the source
            elem_residue[:,i]=elem_Per # Initial el content of peridotite mantle
            elem_Melt_Inst[:,i]=0 # No melting, so no el in the instantaneous melt
            elem_Melt_Agg[:,i]=0 # No melting, so no el in the aggregated melt
            DeltaM[i]=0 # No change in the amount of peridotite residue
            DeltaXSulf[:,i]=0 # No change in the proportion of sulfide
            DeltaXSil[:,i]=0  # No change in the proportion of silicate
            P[:,i]=0          # No melting, so no preferential melting of different minerals
            Kd_el[:,i]=KdelSil* (1- XSulf[:,i]) + KdelSulf* (XSulf[:,i]) # Bulk Kd for the mantle before any melting occurs
            X_FStep[i]=0 # No melt, so no melt proportion made in this step
            X_FStart[i]=0 # as above
            DeltaXMeltTot[:,i]=0 # as above
            S_Melt_Inst[:, i]=S_Melt[i]
        else:
            DeltaM[i]=M[i-1]-M[i]  # Equation 1
            X_FStep[i]= DeltaM[i]/M[i-1] # Equation 2
            S_residue[:,i] = (S_residue[:,i-1] -S_Melt[i]* X_FStep[i])/(1-X_FStep[i]) # Equation 3
            S_residue[S_residue<0]=0 # To avoid negative numbers.
            XSulf[:,i] = S_residue[:,i]/S_Sulf # Equation 4

            DeltaXSulf[:,i]=(S_residue[:,i-1]-S_residue[:,i])/S_Sulf # Equation 5
            DeltaXSil[:,i]=DeltaM[i]- DeltaXSulf[:,i] # Equation 6

            Kd_el[:,i]=KdelSil* (1- XSulf[:,i]) + KdelSulf* (XSulf[:,i]) # Equation 8
            P[:,i] = KdelSil* ( DeltaXSil[:,i]/DeltaM[i] ) + (KdelSulf * (DeltaXSulf[:,i]/DeltaM[i])) # Equation 9

            elem_Melt_Inst[:,i]=elem_residue[:,i-1]/ ( Kd_el[:,i-1] + (1-P[:,i])*X_FStep[i] ) # Equation 7
            elem_residue[:,i]=( elem_residue[:,i-1] - elem_Melt_Inst[:,i] * X_FStep[i] )/( 1- X_FStep[i] ) # Equation 10
            X_FStart[i]=DeltaM[i]/(1-M[i]) #Equation 13
            elem_Melt_Agg[:,i]=(elem_Melt_Agg[:,i-1]*(1-X_FStart[i])) + (elem_Melt_Inst[:,i]* X_FStart[i]) # Equation 12

            # Sulfur
            if S_residue[:, i]>0:
                S_Melt_Inst[:, i]=S_Melt[i]
            else:
                S_Melt_Inst[:, i]=0

            S_Melt_Agg[:, i]=(S_Melt_Agg[:,i-1]*(1-X_FStart[i])) + (S_Melt_Inst[:, i]* X_FStart[i]) # Equation 12





    df_out_S0=pd.DataFrame(data={'F': 1-M,
                        'M': M,
                    '{}_KD'.format(elx): Kd_el[0,:],
                    '{}_Melt_Agg'.format(elx): elem_Melt_Agg[0,:],
                '{}_Melt_Inst'.format(elx): elem_Melt_Inst[0,:],
                '{}_Residue'.format(elx): elem_residue[0,:],
                'S_Residue': S_residue[0,:],
                'S_Melt_Inst': S_Melt_Inst[0,:],
                'S_Melt_Agg': S_Melt_Agg[0,:],
                  'S_Melt_input': S_Melt,
                    'XSulf': XSulf[0,:]
                              })


                #   ,
                #

    return df_out_S0


def calculate_inst_melts_Thermocalc(Thermocalc_df):
    """ Calculates instantaneous melts from Thermocalc aggregated melt outputs
    Parameters
    ----------------
    Thermocalc_df: pd.DataFrame
        Dataframe of thermocalc outputs (aggregated melts)

    Returns
    ----------------
    pd.DataFrame: Instantaneous liquids.

    """
    SiO2_EJ=Thermocalc_df['SiO2']
    #TiO2_EJ=Thermocalc_df['TiO2'] - No Ti, ask Ellie
    Al2O3_EJ=Thermocalc_df['Al2O3']
    CaO_EJ=Thermocalc_df['CaO']
    MgO_EJ=Thermocalc_df['MgO']
    FeO_EJ=Thermocalc_df['FeO']
    Na2O_EJ=Thermocalc_df['Na2O']
    Fe2O3_EJ=Thermocalc_df['Fe2O3']
    Cr2O3_EJ=Thermocalc_df['Cr2O3']
    F_EJ=Thermocalc_df['liq']

    SiO2_Inst_EJ = [0]*len(SiO2_EJ)
    #TiO2_Inst_EJ = [0]*len(SiO2_EJ)
    Al2O3_Inst_EJ = [0]*len(SiO2_EJ)
    CaO_Inst_EJ = [0]*len(SiO2_EJ)
    MgO_Inst_EJ = [0]*len(SiO2_EJ)
    FeO_Inst_EJ = [0]*len(SiO2_EJ)
    Na2O_Inst_EJ = [0]*len(SiO2_EJ)
    Fe2O3_Inst_EJ = [0]*len(SiO2_EJ)
    Cr2O3_Inst_EJ = [0]*len(SiO2_EJ)
    X = [0]*len(SiO2_EJ)

    for i in range(len(SiO2_EJ)):
        if i == 0:
            SiO2_Inst_EJ[i] = SiO2_EJ[i]
            #TiO2_Inst_EJ[i] = TiO2_EJ[i]
            Al2O3_Inst_EJ[i] = Al2O3_EJ[i]
            CaO_Inst_EJ[i] = CaO_EJ[i]
            MgO_Inst_EJ[i] = MgO_EJ[i]
            FeO_Inst_EJ[i] = FeO_EJ[i]
            Na2O_Inst_EJ[i] = Na2O_EJ[i]
            Fe2O3_Inst_EJ[i] = Fe2O3_EJ[i]
            Cr2O3_Inst_EJ[i] = Cr2O3_EJ[i]
            X[i] = 1
        elif i > 0:
            X[i] = (F_EJ[i] - F_EJ[i-1]) / F_EJ[i]
            SiO2_Inst_EJ[i] = (SiO2_EJ[i] - (1 - X[i])*SiO2_EJ[i-1]) / X[i]
            #TiO2_Inst_EJ[i] = (TiO2_EJ[i] - (1 - X[i])*TiO2_EJ[i-1]) / X[i]
            Al2O3_Inst_EJ[i] = (Al2O3_EJ[i] - (1 - X[i])*Al2O3_EJ[i-1]) /X[i]
            CaO_Inst_EJ[i] = (CaO_EJ[i] - (1 - X[i])*CaO_EJ[i-1]) / X[i]
            MgO_Inst_EJ[i] = (MgO_EJ[i] - (1 - X[i])*MgO_EJ[i-1]) / X[i]
            FeO_Inst_EJ[i] = (FeO_EJ[i] - (1 - X[i])*FeO_EJ[i-1]) / X[i]
            Na2O_Inst_EJ[i] = (Na2O_EJ[i] - (1 - X[i])*Na2O_EJ[i-1]) / X[i]
            Fe2O3_Inst_EJ[i] = (Fe2O3_EJ[i] - (1 - X[i])*Fe2O3_EJ[i-1]) / X[i]
            Cr2O3_Inst_EJ[i] = (Cr2O3_EJ[i] - (1 - X[i])*Cr2O3_EJ[i-1]) / X[i]
    Inst_Liq=pd.DataFrame(data={'SiO2_Liq':SiO2_Inst_EJ,
                                'Al2O3_Liq':Al2O3_Inst_EJ,
                                'CaO_Liq':CaO_Inst_EJ,
                                'MgO_Liq':MgO_Inst_EJ,
                                'FeO_Liq':FeO_Inst_EJ,
                                'Na2O_Liq':Na2O_Inst_EJ,
                                'Fe2O3_Liq':Fe2O3_Inst_EJ})
    Inst_Liq['FeOt_Liq']=Inst_Liq['FeO_Liq'] + 0.8998*Inst_Liq['Fe2O3_Liq']

    return Inst_Liq