import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings as w

## Mantle melting, adapted from Lee et al. (2012) and Wieser et al. (2020)

def Lee_Wieser_sulfide_melting(*,  M_Max=0.01, Modes, KDs,
                        N=3001,
                        S_Sulf=363636, elem_Per=30,
                        S_Mantle=[300, 200, 100, 50, 0],
                        S_Melt_SCSS_2=980.7872088,
                        Prop_S6=0,
                        KdCuSulf=800,
                        KdBaSulf=0,
                        Trace1_per=6.85,
                        ):
    """ Function for calculating trajectory of S during mantle melting,
    adapted from the Excel workbook of Lee et al. (2012) and the Python code
    of Wieser et al. (2020)

    Parameters
    -----------
    N: int
        Number of iterative steps for melting.
    M_Max: float
        Amount of peridotite left at final step. E.g. M_Max=0.7
        calculates 0-30% melting (F=1-M)
    Modes: Pandas.Dataframe
        Dataframe of ol-opx-cpx-sp modes, with column headings 'ol',
        'opx', 'cpx', 'sp'


    """
    elx=KDs['element'][0]
    KdCuSulf=KDs['sulf']
    # Mass of peridotite
    if len(Modes)==1:
        M=np.linspace(1, M_Max, N)
    else:
        if len(Modes)!=N:
            w.warn('You have inputted a dataframe of silicate modes that doesnt match the number of steps you asked for. We are changing the number of steps to match the length of your dataframe')
            M=M=np.linspace(1, M_Max, len(Modes))
        else:
            M=M=np.linspace(1, M_Max, len(Modes))




    if isinstance(S_Melt_SCSS_2, float) or isinstance(S_Melt_SCSS_2, int):
        S_Melt_SCSS_2=np.repeat(S_Melt_SCSS_2, len(M))

    S_Melt = S_Melt_SCSS_2/(1+Prop_S6)

    # Setting up variables to be filled in loop
    #These ones have only 1 dimension for different M steps, as they do not vary with different mantle S contents
    DeltaM= np.empty(len(M), dtype=float)
    X_FStep=np.empty(len(M), dtype=float)
    X_FStart=np.empty( len(M), dtype=float)

# These ones have 2 dimensions, one for different M steps, one for different mantle S contents
    S_residue= np.empty( (len(S_Mantle), len(M)), dtype=float)
    XSulf= np.empty( (len(S_Mantle), len(M)), dtype=float)
    elem_residue= np.empty( (len(S_Mantle), len(M)), dtype=float)
    elem_Melt_Inst= np.empty( (len(S_Mantle), len(M)), dtype=float)
    elem_Melt_Agg= np.empty( (len(S_Mantle), len(M)), dtype=float)
    DeltaXSulf=np.empty( (len(S_Mantle), len(M)), dtype=float)
    DeltaXSil=np.empty( (len(S_Mantle), len(M)), dtype=float)
    P=np.empty( (len(S_Mantle), len(M)), dtype=float)
    Kd_Cu=np.empty( (len(S_Mantle), len(M)), dtype=float)
    DeltaXMeltTot=np.empty( (len(S_Mantle), len(M)), dtype=float)

    # Loop for Cu
    for i in range(0, len(M)):
        if len(Modes)>1:
            KdCuSil=Modes.iloc[i].mul(KDs).sum(axis=1)
            KdBaSil=Modes.iloc[i].mul(KDs).sum(axis=1)
        else:
            KdCuSil=Modes.mul(KDs).sum(axis=1)
            KdBaSil=Modes.mul(KDs).sum(axis=1)


        if i==0: # Setting up the system before any melting takes place
            S_residue[:,i]=S_Mantle # Initial S contents of mantle before any melting set at those defined in S_Mantle
            XSulf[:,i]=S_residue[:,i]/S_Sulf # Proportion (by mass) of sulfides in the source
            elem_residue[:,i]=elem_Per # Initial Cu content of peridotite mantle
            elem_Melt_Inst[:,i]=0 # No melting, so no Cu in the instantaneous melt
            elem_Melt_Agg[:,i]=0 # No melting, so no Cu in the aggregated melt
            DeltaM[i]=0 # No change in the amount of peridotite residue
            DeltaXSulf[:,i]=0 # No change in the proportion of sulfide
            DeltaXSil[:,i]=0  # No change in the proportion of silicate
            P[:,i]=0          # No melting, so no preferential melting of different minerals
            Kd_Cu[:,i]=KdCuSil* (1- XSulf[:,i]) + KdCuSulf* (XSulf[:,i]) # Bulk Kd for the mantle before any melting occurs
            X_FStep[i]=0 # No melt, so no melt proportion made in this step
            X_FStart[i]=0 # as above
            DeltaXMeltTot[:,i]=0 # as above
        else:
            DeltaM[i]=M[i-1]-M[i]  # Equation 1
            X_FStep[i]= DeltaM[i]/M[i-1] # Equation 2
            S_residue[:,i] = (S_residue[:,i-1] -S_Melt[i]* X_FStep[i])/(1-X_FStep[i]) # Equation 3
            S_residue[S_residue<0]=0 # To avoid negative numbers.
            XSulf[:,i] = S_residue[:,i]/S_Sulf # Equation 4

            DeltaXSulf[:,i]=(S_residue[:,i-1]-S_residue[:,i])/S_Sulf # Equation 5
            DeltaXSil[:,i]=DeltaM[i]- DeltaXSulf[:,i] # Equation 6

            Kd_Cu[:,i]=KdCuSil* (1- XSulf[:,i]) + KdCuSulf* (XSulf[:,i]) # Equation 8
            P[:,i] = KdCuSil* ( DeltaXSil[:,i]/DeltaM[i] ) + (KdCuSulf * (DeltaXSulf[:,i]/DeltaM[i])) # Equation 9

            elem_Melt_Inst[:,i]=elem_residue[:,i-1]/ ( Kd_Cu[:,i-1] + (1-P[:,i])*X_FStep[i] ) # Equation 7
            elem_residue[:,i]=( elem_residue[:,i-1] - elem_Melt_Inst[:,i] * X_FStep[i] )/( 1- X_FStep[i] ) # Equation 10
            X_FStart[i]=DeltaM[i]/(1-M[i]) #Equation 13
            elem_Melt_Agg[:,i]=(elem_Melt_Agg[:,i-1]*(1-X_FStart[i])) + (elem_Melt_Inst[:,i]* X_FStart[i]) # Equation 12




    df_out_S0=pd.DataFrame(data={'F': 1-M,
                        'M': M,
                    '{}_KD'.format(elx): Kd_Cu[0,:],
                    '{}_Melt_Agg'.format(elx): elem_Melt_Agg[0,:],
                '{}_Melt_Inst'.format(elx): elem_Melt_Inst[0,:],
                '{}_Residue'.format(elx): elem_residue[0,:],
                'S_Residue': S_residue[0,:],
                  'S_Melt': S_Melt,
                    'XSulf': XSulf[0,:]
                              })


                #   ,
                #

    return df_out_S0
