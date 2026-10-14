import inspect
import warnings as w
import numpy as np
import pandas as pd

from PySulfSat.scss_calcs2 import (
calculate_S2017_SCSS,
calculate_O2021_SCSS,
calculate_OM2022_SCSS,
calculate_LiZhang2022_SCSS,
calculate_Liu2021_SCSS,
calculate_B2021_SCSS,
calculate_F2015_SCSS)

from PySulfSat.scas_calc import (calculate_CD2019_SCAS, calculate_ZT2019_SCAS, calculate_MK2015_SCAS)


## Run every SCSS or SCAS model on the same melts, with one column per model.
# Each entry: output column name -> (function, column it returns, extra keyword arguments)
_SCSS_MODELS = {
    'SCSS2_ppm_Smythe2017': (calculate_S2017_SCSS, 'SCSS2_ppm_ideal_Smythe2017', {}),
    'SCSS2_ppm_ONeill2021': (calculate_O2021_SCSS, 'SCSS2_ppm', {}),
    'SCSS2_ppm_ONeillMavrogenes2022': (calculate_OM2022_SCSS, 'SCSS2_ppm', {}),
    'SCSS_Tot_LiZhang2022': (calculate_LiZhang2022_SCSS, 'SCSS_Tot', {}),
    'SCSS2_ppm_Liu2021': (calculate_Liu2021_SCSS, 'SCSS2_ppm', {}),
    'SCSS2_ppm_Blanchard2021_eq11': (calculate_B2021_SCSS, 'SCSS2_ppm_eq11', {}),
    'SCSS2_ppm_Blanchard2021_eq12': (calculate_B2021_SCSS, 'SCSS2_ppm_eq12', {}),
    'SCSS2_ppm_Fortin2015': (calculate_F2015_SCSS, 'SCSS2_ppm', {}),
}

_SCAS_MODELS = {
    'SCAS6_ppm_ChowdhuryDasgupta2019': (calculate_CD2019_SCAS, 'SCAS6_ppm', {}),
    'SCAS6_ppm_ZajaczTsay2019': (calculate_ZT2019_SCAS, 'SCAS6_ppm', {}),
    'SCAS6_ppm_MasottaKeppler2015': (calculate_MK2015_SCAS, 'SCAS6_ppm', {}),
}


def _run_models(models, df, kwargs, include, quiet):
    """Calls each model with only the arguments it accepts. A model that fails gives NaN and a warning."""
    if include is not None:
        unknown = [m for m in include if m not in models]
        if unknown:
            raise ValueError(f'Unknown model(s) {unknown}. Choose from {list(models)}')
        models = {k: v for k, v in models.items() if k in include}
    out = pd.DataFrame(index=range(len(df)))
    failed = {}
    for name, (fn, col, extra) in models.items():
        params = inspect.signature(fn).parameters
        kw = {k: v for k, v in {**kwargs, **extra}.items() if k in params and v is not None}
        try:
            with w.catch_warnings():
                if quiet:
                    w.simplefilter('ignore')
                res = fn(df=df.copy(), **kw)
            out[name] = np.asarray(res[col], dtype=float)
        except Exception as e:  # e.g. a model needing sulfide or H2O inputs that weren't given
            out[name] = np.nan
            failed[name] = f'{type(e).__name__}: {e}'
    for name, msg in failed.items():
        w.warn(f'{name} could not be calculated ({msg}); filled with NaN.')
    return out


def calculate_all_SCSS(*, df, T_K, P_kbar, H2O_Liq=None, Fe3Fet_Liq=None,
                       Fe_FeNiCu_Sulf=None, Cu_FeNiCu_Sulf=None, Ni_FeNiCu_Sulf=None,
                       Fe_Sulf=None, Ni_Sulf=None, Cu_Sulf=None, Ni_Liq=None, Cu_Liq=None,
                       T_K_transition=True, highT=False, lowT=False, models=None, quiet=True):
    '''
    Calculates the SCSS with every SCSS model in PySulfSat, for the same melts and conditions.

    Parameters
    -------
    df: pandas.DataFrame
        Liquid compositions with column headings SiO2_Liq, FeOt_Liq etc. (e.g. from import_data).
    T_K, P_kbar: float or pandas.Series
        Temperature in Kelvin and pressure in kbar.
    H2O_Liq, Fe3Fet_Liq: float or pandas.Series, optional
        Overwrite the H2O_Liq and Fe3Fet_Liq columns of df. Each model uses them only if it needs them.
    Fe_FeNiCu_Sulf, Cu_FeNiCu_Sulf, Ni_FeNiCu_Sulf, Fe_Sulf, Ni_Sulf, Cu_Sulf, Ni_Liq, Cu_Liq:
        Sulfide composition, as for the individual SCSS functions: a fixed Fe_FeNiCu_Sulf, Fe_FeNiCu_Sulf='Calc_ONeill'
        or 'Calc_Smythe' with Ni_Liq and Cu_Liq, or a measured sulfide (Fe_Sulf, Ni_Sulf, Cu_Sulf).
        Fortin et al. (2015) does not use the sulfide composition.
    T_K_transition, highT, lowT: bool
        Passed to Li and Zhang (2022), as in calculate_LiZhang2022_SCSS.
    models: list, optional
        Only calculate these (column names below). Default: all.
    quiet: bool
        If True (default), hides the warnings each model prints; models that fail are still reported.

    Returns
    -------
    pandas.DataFrame
        One column per model, followed by df:
        SCSS2_ppm_Smythe2017, SCSS2_ppm_ONeill2021, SCSS2_ppm_ONeillMavrogenes2022, SCSS_Tot_LiZhang2022,
        SCSS2_ppm_Liu2021, SCSS2_ppm_Blanchard2021_eq11, SCSS2_ppm_Blanchard2021_eq12, SCSS2_ppm_Fortin2015.
        Note Li and Zhang (2022) returns total S at sulfide saturation (SCSS_Tot); the others return S2-.
        A model that cannot be calculated with the inputs given is filled with NaN, with a warning.
    '''
    df_c = df.reset_index(drop=True)
    kwargs = dict(T_K=T_K, P_kbar=P_kbar, H2O_Liq=H2O_Liq, Fe3Fet_Liq=Fe3Fet_Liq,
                  Fe_FeNiCu_Sulf=Fe_FeNiCu_Sulf, Cu_FeNiCu_Sulf=Cu_FeNiCu_Sulf, Ni_FeNiCu_Sulf=Ni_FeNiCu_Sulf,
                  Fe_Sulf=Fe_Sulf, Ni_Sulf=Ni_Sulf, Cu_Sulf=Cu_Sulf, Ni_Liq=Ni_Liq, Cu_Liq=Cu_Liq,
                  T_K_transition=T_K_transition, highT=highT, lowT=lowT)
    out = _run_models(_SCSS_MODELS, df_c, kwargs, models, quiet)
    return pd.concat([out, df_c], axis=1)


def calculate_all_SCAS(*, df, T_K, H2O_Liq=None, Fe3Fet_Liq=None, P_kbar=None, models=None, quiet=True):
    '''
    Calculates the SCAS with every SCAS model in PySulfSat, for the same melts and conditions.

    Parameters
    -------
    df: pandas.DataFrame
        Liquid compositions with column headings SiO2_Liq, FeOt_Liq etc.
    T_K: float or pandas.Series
        Temperature in Kelvin.
    H2O_Liq, Fe3Fet_Liq, P_kbar: optional
        Passed to each model that uses them.
    models: list, optional
        Only calculate these (column names below). Default: all.
    quiet: bool
        If True (default), hides the warnings each model prints; models that fail are still reported.

    Returns
    -------
    pandas.DataFrame
        SCAS6_ppm_ChowdhuryDasgupta2019, SCAS6_ppm_ZajaczTsay2019, SCAS6_ppm_MasottaKeppler2015, followed by df.
    '''
    df_c = df.reset_index(drop=True)
    kwargs = dict(T_K=T_K, H2O_Liq=H2O_Liq, Fe3Fet_Liq=Fe3Fet_Liq, P_kbar=P_kbar)
    out = _run_models(_SCAS_MODELS, df_c, kwargs, models, quiet)
    return pd.concat([out, df_c], axis=1)
