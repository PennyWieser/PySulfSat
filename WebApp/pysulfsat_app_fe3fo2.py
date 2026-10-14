""" PySulfSat web app: sulfide and sulfate saturation, S6+/ST and total S for silicate melts, from spreadsheets,
Petrolog3 outputs or MELTS .tbl files. Follows the PySulfSat example notebooks in docs/Examples.
"""

import io
import os
import tempfile
import warnings
import contextlib
from pathlib import Path

import numpy as np
import pandas as pd
import altair as alt
import streamlit as st
import PySulfSat as ss

st.set_page_config(page_title="PySulfSat", page_icon="🌋", initial_sidebar_state="expanded")

HERE = Path(__file__).resolve().parent
EXAMPLES = HERE.parent / "docs" / "Examples"

# PySulfSat's PyPI release (setup.py package_data only lists *.pkl) leaves out Cali_datasets_CSV, so
# return_cali_dataset fails when installed from PyPI. The deployed repository has the CSVs in src/PySulfSat,
# so point PySulfSat there when its own copy is missing.
import PySulfSat.cali_plots as _cali_plots
if not (Path(_cali_plots.PySulfSat_dir) / "Cali_datasets_CSV").is_dir():
    _repo_pkg = HERE.parent / "src" / "PySulfSat"
    if (_repo_pkg / "Cali_datasets_CSV").is_dir():
        _cali_plots.PySulfSat_dir = _repo_pkg

# PySulfSat 1.0.13 (PyPI) fails in calculate_S_Total_SCSS_SCAS when S6St_Liq is given without a model
# ("'GW2026' in model" with model=None); fixed on GitHub after that release. Patch that one line if present.
import inspect as _inspect
import textwrap as _textwrap
import PySulfSat.s6_corrections as _s6c
try:
    _src = _textwrap.dedent(_inspect.getsource(_s6c.calculate_S_Total_SCSS_SCAS))
    _bad = "if model=='OM2022' or model =='BW2023' or 'GW2026' in model:"
    if _bad in _src:
        _ns = {}
        exec(_src.replace(_bad, "if model is not None and (model=='OM2022' or model =='BW2023' or 'GW2026' in model):"),
             _s6c.__dict__, _ns)
        _s6c.calculate_S_Total_SCSS_SCAS = ss.calculate_S_Total_SCSS_SCAS = _ns["calculate_S_Total_SCSS_SCAS"]
except (OSError, TypeError):
    pass
MAX_ROWS = 50_000
SAMPLE_COLS = ["Sample_ID_Liq", "Sample_ID", "Sample", "SampleID", "Sample Name", "Sample ID"]
OXIDES = ["SiO2", "TiO2", "Al2O3", "FeOt", "MnO", "MgO", "CaO", "Na2O", "K2O", "P2O5", "H2O"]

try:
    import Thermobar  # noqa: F401  # Gorojovsky and Wood (2026) uses its Kress and Carmichael (1991) redox model
    HAS_THERMOBAR = True
except ImportError:
    HAS_THERMOBAR = False

CITATION = "Wieser and Gleeson (2023), Volcanica, 6(1), 107-127. doi: 10.30909/vol.06.01.107127."

# %% ----------------------------------------------------------------
# models

SCSS_MODELS = {
    "Smythe et al. (2017)": dict(fn=ss.calculate_S2017_SCSS, out="SCSS2_ppm_ideal_Smythe2017",
                                 err="SCSS2_ppm_ideal_Smythe2017_1sigma", h2o=True, sulf=True, cali="S2017_SCSS"),
    "O'Neill (2021)": dict(fn=ss.calculate_O2021_SCSS, out="SCSS2_ppm", h2o=False, sulf=True, cali="O2021_SCSS"),
    "O'Neill and Mavrogenes (2022)": dict(fn=ss.calculate_OM2022_SCSS, out="SCSS2_ppm", h2o=False, sulf=True),
    "Li and Zhang (2022)": dict(fn=ss.calculate_LiZhang2022_SCSS, out="SCSS_Tot", h2o=True, sulf=True,
                                cali="LiZhang2022_SCSS"),
    "Liu et al. (2021)": dict(fn=ss.calculate_Liu2021_SCSS, out="SCSS2_ppm", h2o=True, sulf=True, cali="Liu2021_SCSS"),
    "Blanchard et al. (2021), eq. 11": dict(fn=ss.calculate_B2021_SCSS, out="SCSS2_ppm_eq11", h2o=True, sulf=True,
                                            cali="B2021_SCSS"),
    "Blanchard et al. (2021), eq. 12": dict(fn=ss.calculate_B2021_SCSS, out="SCSS2_ppm_eq12", h2o=True, sulf=True,
                                            cali="B2021_SCSS"),
    "Fortin et al. (2015)": dict(fn=ss.calculate_F2015_SCSS, out="SCSS2_ppm", h2o=True, sulf=False, cali="F2015_SCSS"),
}
SCAS_MODELS = {
    "Chowdhury and Dasgupta (2019)": dict(fn=ss.calculate_CD2019_SCAS, cali="CD2019_SCAS"),
    "Zajacz and Tsay (2019)": dict(fn=ss.calculate_ZT2019_SCAS, cali="ZT2022_SCAS"),
    "Masotta and Keppler (2015)": dict(fn=ss.calculate_MK2015_SCAS, cali="MK_2015"),
}
# S6+/ST models: the redox inputs each accepts, and the S6St column it returns
S6_MODELS = {
    "Nash et al. (2019)": dict(redox=["Fe3"], needs_df=False),
    "Jugo et al. (2010)": dict(redox=["dQFM"], needs_df=False),
    "O'Neill and Mavrogenes (2022)": dict(redox=["logfo2", "Fe3"], needs_df=True),
    "Boulliung and Wood (2022) with O'Neill and Mavrogenes (2022)": dict(redox=["logfo2", "Fe3"], needs_df=True),
    "Gorojovsky and Wood (2026), 1 atm": dict(redox=["logfo2", "Fe3", "dQFM"], needs_df=True, col="S6St_1atm"),
    "Gorojovsky and Wood (2026), P after Boulliung and Wood (2023)": dict(redox=["logfo2", "Fe3", "dQFM"],
                                                                          needs_df=True, col="S6St_BW23"),
    "Gorojovsky and Wood (2026), P after Muth and Cottrell (2025)": dict(redox=["logfo2", "Fe3", "dQFM"],
                                                                         needs_df=True, col="S6St_MC25"),
    "Gorojovsky and Wood (2026), P after Thomas and Wood (2026)": dict(redox=["logfo2", "Fe3", "dQFM"],
                                                                       needs_df=True, col="S6St_TW26"),
}
FROM_FO2 = "Calculate from log fO₂ column"
FROM_FE3 = "Calculate from Fe³⁺/ΣFe (Kress and Carmichael 1991)"
REDOX_LABEL = {"Fe3": "Fe³⁺/ΣFe", "logfo2": "log₁₀ fO₂", "dQFM": "ΔQFM (Frost 1991)"}

# Inputs that can be a single value or a column: label, default value, likely column names, step
PARAMS = {
    "T": ("Temperature", 1200.0, ["T_K", "T_C", "T (C)", "Temperature", "T_K_calc"], 10.0),
    "P": ("Pressure (kbar)", 1.0, ["P_kbar", "Pressure(kbar)", "P_kbar_calc"], 0.5),
    "H2O": ("H₂O in the liquid (wt%)", 0.0, ["H2O_Liq"], 0.1),
    "Fe3": ("Fe³⁺/ΣFe", 0.15, ["Fe3Fet_Liq"], 0.01),
    "logfo2": ("log₁₀ fO₂", -8.0, ["Lg(fO2)", "logfo2", "logfO2", "log_fO2"], 0.1),
    "dQFM": ("ΔQFM (Frost 1991)", 0.0, ["DeltaQFM", "deltaQFM", "dQFM", "DeltaQFM_Frost1991"], 0.1),
    "S6": ("S⁶⁺/ΣS", 0.1, ["S6St_Liq", "S6_ST", "S6/ST"], 0.01),
    "FeFeNiCu": ("Fe/(Fe+Ni+Cu) in the sulfide", 0.65, ["Fe_FeNiCu_Sulf"], 0.01),
    "Ni": ("Ni in the liquid (ppm)", 100.0, ["Ni_Liq_ppm", "Ni_Liq", "Ni_ppm"], 10.0),
    "Cu": ("Cu in the liquid (ppm)", 100.0, ["Cu_Liq_ppm", "Cu_Liq", "Cu_ppm"], 10.0),
    "Fe_Sulf": ("Fe in the sulfide (wt%)", 50.0, ["Fe_Sulf", "Fe_sulf_wt%"], 1.0),
    "Ni_Sulf": ("Ni in the sulfide (wt%)", 5.0, ["Ni_Sulf", "Ni_sulf_wt%"], 1.0),
    "Cu_Sulf": ("Cu in the sulfide (wt%)", 5.0, ["Cu_Sulf", "Cu_sulf_wt%"], 1.0),
}

CALCS = ("Sulfide saturation (SCSS)", "Sulfate saturation (SCAS)", "Total S at sulfide and sulfate saturation",
         "S⁶⁺/ΣS in the melt", "Compare with calibration datasets", "Convert S units")
EXAMPLE_FILES = {
    "Petrolog3 fractionation path (Hawaii)": ("Petrolog3 output", "Intro_Example_Petrolog_FC/PetrologCalculations.xlsx", None),
    "MELTS liquid .tbl": ("MELTS .tbl file", "Data_Input/melts-liquid.tbl", None),
    "Glass analyses (spreadsheet)": ("Spreadsheet", "Data_Input/Glass_input_example.xlsx", "Glass_input"),
}


# %% ----------------------------------------------------------------
# helpers

def call(fn, **kw):
    """Calls a PySulfSat function with only the arguments it accepts, capturing warnings and printed messages."""
    import inspect
    params = inspect.signature(fn).parameters
    kw = {k: v for k, v in kw.items() if k in params and v is not None}
    printed = io.StringIO()
    with warnings.catch_warnings(record=True) as caught, contextlib.redirect_stdout(printed):
        warnings.simplefilter("always")
        out = fn(**kw)
    msgs = sorted({str(w.message) for w in caught if issubclass(w.category, UserWarning)})
    return out, msgs, printed.getvalue()


def clean(df):
    df.columns = df.columns.astype(str).str.strip()
    df = df.loc[:, ~df.columns.str.startswith("Unnamed")]
    return df.dropna(how="all").reset_index(drop=True)


def load_data(fmt, data, name, sheet):
    """Liquid compositions through PySulfSat's own importers: Petrolog3 and MELTS files need the file on disk."""
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        if fmt == "Spreadsheet":
            raw = clean(pd.read_csv(io.BytesIO(data), encoding="utf-8-sig") if name.lower().endswith(".csv")
                        else pd.read_excel(io.BytesIO(data), sheet_name=sheet or 0))
            has_suffix = any(f"{ox}_Liq" in raw for ox in OXIDES)
            if not has_suffix and not any(ox in raw for ox in OXIDES):
                raise ValueError(f"No oxide columns found (e.g. SiO2_Liq or SiO2). Columns in your file: {list(raw.columns)}")
            if not has_suffix:
                raw = raw.rename(columns=lambda c: c if c.endswith("_Liq") else c + "_Liq")
            return ss.import_dataframe(raw)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, Path(name).name)
            Path(path).write_bytes(data)
            if fmt == "Petrolog3 output":
                return ss.import_data(path, sheet_name=sheet, Petrolog=True)
            return ss.import_data(path, MELTS=True)


def numeric_columns(df):
    return [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c]) and not df[c].isna().all()]


def ask(kind, cols, key, prefer_file=True, help=None):
    """A single value, or a column from the file. Returns a float, or (column, unit) for resolve()."""
    label, value, guess, step = PARAMS[kind]
    found = next((c for c in guess if c in cols), None)
    opts = ["Column in file", "Single value"] if cols else ["Single value"]
    fo2_found = next((c for c in PARAMS["logfo2"][2] if c in cols), None)
    fe3_found = next((c for c in PARAMS["Fe3"][2] if c in cols), None)
    if kind in ("dQFM", "logfo2") and HAS_THERMOBAR:
        opts = opts[:-1] + ([FROM_FO2] if kind == "dQFM" and cols else []) + [FROM_FE3, "Single value"]
    index = opts.index("Column in file") if (found and prefer_file) else \
        opts.index(FROM_FO2) if (FROM_FO2 in opts and fo2_found) else opts.index("Single value")
    src = st.radio(label, opts, index=index, horizontal=len(opts) < 3, key=f"{key}_src", help=help)
    if src == FROM_FE3:
        st.caption("log₁₀ fO₂ from Fe³⁺/ΣFe and the melt composition with Kress and Carmichael (1991), via "
                   "Thermobar's convert_fe_partition_to_fo2" + (" (ΔQFM relative to Frost 1991)" if kind == "dQFM" else "")
                   + ", at this Fe³⁺/ΣFe, temperature and pressure:")
        fe3 = ask("Fe3", cols, f"{key}_kcfe3")
        T = ask("T", cols, f"{key}_kcT")
        P = ask("P", cols, f"{key}_kcP")
        return ("fe3", fe3, T, P, kind)
    if src == FROM_FO2:
        fcol = st.selectbox("log₁₀ fO₂ column", cols, index=cols.index(fo2_found) if fo2_found else 0,
                            key=f"{key}_fo2col", help="log₁₀ fO₂ in bars, e.g. Lg(fO2) from Petrolog3.")
        st.caption("Converted to ΔQFM with Thermobar's convert_fo2_to_buffer, at this temperature and pressure:")
        T = ask("T", cols, f"{key}_fo2T")
        P = ask("P", cols, f"{key}_fo2P")
        buf = st.radio("QFM buffer", ["Frost (1991)", "O'Neill (1987)"], horizontal=True, key=f"{key}_fo2buf")
        return ("fo2", fcol, T, P, buf)
    if src == "Single value":
        if kind == "T":
            unit = st.radio("Units", ["°C", "K"], horizontal=True, key=f"{key}_unit", label_visibility="collapsed")
            v = st.number_input("Temperature", value=value if unit == "°C" else value + 273.15, step=step,
                                key=f"{key}_val_{unit}", label_visibility="collapsed")
            return v + 273.15 if unit == "°C" else v
        return st.number_input(label, value=value, step=step, key=f"{key}_val", label_visibility="collapsed",
                               format="%.4f" if step < 0.1 else "%.2f")
    col = st.selectbox("Column", cols, index=cols.index(found) if found in cols else 0, key=f"{key}_col",
                       label_visibility="collapsed")
    unit = None
    if kind == "T":
        is_k = col.upper().startswith("T_K") or "(K)" in col
        unit = st.radio("Units", ["K", "°C"], index=0 if is_k else 1, horizontal=True, key=f"{key}_unit_col")
    return col, unit


def resolve(value, df):
    """Turns ask()'s answer into a float or a Series aligned with df (temperatures in K)."""
    if not isinstance(value, tuple):
        return value
    if value[0] == "fo2":
        return fo2_to_dqfm(value, df)
    if value[0] == "fe3":
        return fe3_to_redox(value, df)
    col, unit = value
    s = pd.to_numeric(df[col], errors="coerce").reset_index(drop=True)
    return s + 273.15 if unit == "°C" else s


def fo2_to_dqfm(value, df):
    """ΔQFM from a log10 fO2 column, with Thermobar's convert_fo2_to_buffer (Frost 1991 or O'Neill 1987 QFM)."""
    import Thermobar as pt
    _, fcol, T, P, buf = value
    n = len(df)
    logfo2 = pd.to_numeric(df[fcol], errors="coerce").reset_index(drop=True)
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        out = pt.convert_fo2_to_buffer(fo2=10 ** logfo2, T_K=series(resolve(T, df), n), P_kbar=series(resolve(P, df), n))
    return out["DeltaQFM_Frost1991" if buf.startswith("Frost") else "DeltaQFM_ONeill1987"].reset_index(drop=True)


def fe3_to_redox(value, df):
    """log10 fO2 (or ΔQFM, Frost 1991) from Fe3+/ΣFe with Kress and Carmichael (1991), via Thermobar."""
    import Thermobar as pt
    _, fe3, T, P, kind = value
    n = len(df)
    liq = df.reset_index(drop=True)
    liq = liq[[c for c in liq.columns if c.endswith("_Liq") and c not in ("FeO_Liq", "Fe2O3_Liq", "Fe3Fet_Liq")]]
    liq = liq.apply(pd.to_numeric, errors="coerce")
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        out = pt.convert_fe_partition_to_fo2(liq_comps=liq, T_K=series(resolve(T, df), n),
                                             P_kbar=series(resolve(P, df), n), model="Kress1991",
                                             Fe3Fet_Liq=series(resolve(fe3, df), n))
    res = out["DeltaQFM_Frost1991"] if kind == "dQFM" else np.log10(out["fo2_calc"].astype(float))
    return pd.Series(np.asarray(res, dtype=float))


def series(x, n):
    """A Series of length n from a value or Series (some PySulfSat functions need Series)."""
    return x.reset_index(drop=True).astype(float) if isinstance(x, pd.Series) else pd.Series(np.full(n, float(x)))


def sample_col(df):
    return next((c for c in SAMPLE_COLS if c in df.columns and df[c].notna().any()), None)


# %% ----------------------------------------------------------------
# shared settings forms

def sulfide_form(cols, key):
    """Sulfide composition options shared by the SCSS models. Returns a function df -> keyword arguments."""
    how = st.radio("Sulfide composition", ["Fixed Fe/(Fe+Ni+Cu)", "Calculate (O'Neill 2021)",
                                           "Calculate (Smythe et al. 2017)", "Measured sulfide (Fe, Ni, Cu wt%)"],
                   key=f"{key}_how", help="O'Neill (2021): empirical, from FeOt, Ni, Cu and Fe³⁺/ΣFe of the melt. "
                                          "Smythe et al. (2017): solves for the sulfide in equilibrium with the "
                                          "melt's Ni and Cu, using Kiseeva and Wood (2015) partition coefficients.")
    if how.startswith("Fixed"):
        r = ask("FeFeNiCu", cols, f"{key}_ratio", help="E.g. 0.65, measured on sulfides in the sample suite.")
        return lambda df: dict(Fe_FeNiCu_Sulf=resolve(r, df))
    if how.startswith("Calculate"):
        ni, cu = ask("Ni", cols, f"{key}_ni"), ask("Cu", cols, f"{key}_cu")
        tag = "Calc_ONeill" if "O'Neill" in how else "Calc_Smythe"
        return lambda df: dict(Fe_FeNiCu_Sulf=tag, Ni_Liq=resolve(ni, df), Cu_Liq=resolve(cu, df))
    fe, ni, cu = ask("Fe_Sulf", cols, f"{key}_fes"), ask("Ni_Sulf", cols, f"{key}_nis"), ask("Cu_Sulf", cols, f"{key}_cus")
    return lambda df: dict(Fe_Sulf=resolve(fe, df), Ni_Sulf=resolve(ni, df), Cu_Sulf=resolve(cu, df))


def s6_form(cols, key, allow=("value",), kleinsasser=False):
    """S6+/ST choice. Returns (label, function(df, T_K, P_kbar) -> S6St Series or None, extra kwargs for total-S)."""
    opts = (["None (all S²⁻)"] if "none" in allow else []) + (["Measured or set value"] if "value" in allow else []) \
        + list(S6_MODELS) + (["Kleinsasser et al. (2022), dacite"] if kleinsasser else [])
    if not HAS_THERMOBAR:
        opts = [o for o in opts if not o.startswith("Gorojovsky")]
    choice = st.selectbox("S⁶⁺/ΣS from", opts, key=f"{key}_s6",
                          help="Gorojovsky and Wood (2026) needs Thermobar installed." if not HAS_THERMOBAR else None)
    if choice.startswith("None"):
        return choice, None, {}
    if choice.startswith("Measured"):
        v = ask("S6", cols, f"{key}_s6v")
        return choice, lambda df, T, P: series(resolve(v, df), len(df)), {}
    if choice.startswith("Kleinsasser"):
        d = ask("dQFM", cols, f"{key}_kdq")
        return choice, None, {"model": "Kleinsasser", "dQFM": d}
    spec = S6_MODELS[choice]
    redox = st.selectbox("Redox given as", spec["redox"], format_func=REDOX_LABEL.get, key=f"{key}_redox") \
        if len(spec["redox"]) > 1 else spec["redox"][0]
    r = ask(redox, cols, f"{key}_r_{redox}")  # separate keys, so a Fe3+ column choice isn't reused for ΔQFM

    used = {}

    def s6(df, T, P):
        n = len(df)
        rv = series(resolve(r, df), n)
        used[{"Fe3": "Fe3Fet_Liq_used", "logfo2": "logfo2_used", "dQFM": "DeltaQFM_used"}[redox]] = rv
        red = {"Fe3": dict(Fe3Fet_Liq=rv), "logfo2": dict(logfo2=rv), "dQFM": dict(deltaQFM=rv)}[redox]
        if choice.startswith("Nash"):
            out, _, _ = call(ss.calculate_S6St_Nash2019, T_K=series(T, n), **red)
            return series(out, n)
        if choice.startswith("Jugo"):
            out, _, _ = call(ss.calculate_S6St_Jugo2010_eq10, **red)
            return series(out, n)
        if choice.startswith("O'Neill"):
            out, _, _ = call(ss.calculate_OM2022_S6St, df=df, T_K=series(T, n), **red)
            return series(out["S6St_Liq"], n)
        if choice.startswith("Boulliung"):
            out, _, _ = call(ss.calculate_BW2022_OM2022_S6St, df=df, T_K=series(T, n), **red)
            return series(out["S6St_Liq"], n)
        out, _, _ = call(ss.calculate_GW2026_S6St, df=df, T_K=series(T, n), P_kbar=series(P, n), **red)
        return series(out[spec["col"]], n)
    return choice, s6, {"used": used}


# %% ----------------------------------------------------------------
# calculations: each settings form returns a job(df) -> (results, key columns, plotted column, warnings)

def form_scss(cols):
    name = st.selectbox("SCSS model", list(SCSS_MODELS), key="scss_model")
    spec = SCSS_MODELS[name]
    st.markdown("**Conditions**")
    T, P = ask("T", cols, "scss_T"), ask("P", cols, "scss_P")
    H2O = ask("H2O", cols, "scss_H2O", help="Used by this model; overrides any H2O_Liq column.") if spec["h2o"] else None
    Fe3 = ask("Fe3", cols, "scss_Fe3", help="Only Fe²⁺ enters the SCSS calculation.") if name != "Fortin et al. (2015)" else None
    lz = {}
    if name.startswith("Li and Zhang"):
        tr = st.radio("Temperature dependence", ["With the 1200 °C transition", "High-T form throughout",
                                                 "Low-T form throughout"], key="scss_lz",
                      help="The Li and Zhang (2022) model changes form at 1200 °C; using one form throughout "
                           "avoids the step.")
        lz = dict(T_K_transition=tr.startswith("With"), highT=tr.startswith("High"), lowT=tr.startswith("Low"))
    sulf = None
    if spec["sulf"]:
        st.markdown("**Sulfide**")
        sulf = sulfide_form(cols, "scss_sulf")
    st.markdown("**S⁶⁺ correction**")
    s6_label, s6, _ = s6_form(cols, "scss", allow=("none", "value"))

    def job(df):
        n = len(df)
        t, p = series(resolve(T, df), n), series(resolve(P, df), n)
        kw = dict(df=df, T_K=t, P_kbar=p, H2O_Liq=resolve(H2O, df) if H2O is not None else None,
                  Fe3Fet_Liq=resolve(Fe3, df) if Fe3 is not None else None, **lz)
        if sulf:
            kw.update(sulf(df))
        out, msgs, printed = call(spec["fn"], **kw)
        out = out.reset_index(drop=True)
        scss = out.pop(spec["out"])  # every model's result as SCSS2_ppm, first
        out = out.drop(columns=["SCSS2_ppm"], errors="ignore")
        out.insert(0, "SCSS2_ppm", scss)
        key = ["SCSS2_ppm"] + ([spec["err"]] if spec.get("err") in out else [])
        for c in ("Fe_FeNiCu_Sulf", "Fe_FeNiCu_Sulf_calc"):
            if c in out and c not in key:
                key.append(c)
        if s6 is not None:
            s6v = s6(df, t, p)
            tot, m2, _ = call(ss.calculate_SCSS_Total, SCSS=out["SCSS2_ppm"], S6St_Liq=s6v)
            out.insert(1, "S6St_Liq", s6v.values)
            out.insert(2, "SCSS_Tot_ppm", np.asarray(tot))
            key = ["SCSS2_ppm", "S6St_Liq", "SCSS_Tot_ppm"] + key[1:]
            msgs += m2
        return out, key, ("SCSS_Tot_ppm" if s6 is not None else "SCSS2_ppm"), msgs, printed
    return job, f"{name}" + (f"; S⁶⁺ correction: {s6_label}" if s6 is not None else "")


def form_scas(cols):
    name = st.selectbox("SCAS model", list(SCAS_MODELS), key="scas_model")
    st.markdown("**Conditions**")
    T = ask("T", cols, "scas_T")
    H2O = ask("H2O", cols, "scas_H2O", help="Overrides any H2O_Liq column.")
    st.markdown("**S²⁻ correction**")
    s6_label, s6, _ = s6_form(cols, "scas", allow=("none", "value"))
    P = ask("P", cols, "scas_P") if s6_label.startswith("Gorojovsky") else 1.0

    def job(df):
        n = len(df)
        t = series(resolve(T, df), n)
        out, msgs, printed = call(SCAS_MODELS[name]["fn"], df=df, T_K=t, H2O_Liq=resolve(H2O, df))
        out = out.reset_index(drop=True)
        key, plot = ["SCAS6_ppm"], "SCAS6_ppm"
        if s6 is not None:
            s6v = s6(df, t, series(resolve(P, df), n))
            tot, m2, _ = call(ss.calculate_SCAS_Total, SCAS=out["SCAS6_ppm"], S6St_Liq=s6v)
            out.insert(out.columns.get_loc("SCAS6_ppm") + 1, "S6St_Liq", s6v.values)
            out.insert(out.columns.get_loc("S6St_Liq") + 1, "SCAS_Tot_ppm", np.asarray(tot))
            key, plot, msgs = ["SCAS6_ppm", "S6St_Liq", "SCAS_Tot_ppm"], "SCAS_Tot_ppm", msgs + m2
        return out, key, plot, msgs, printed
    return job, name + (f"; S²⁻ correction: {s6_label}" if s6 is not None else "")


def form_total(cols):
    st.markdown("Combines an SCSS and an SCAS model with S⁶⁺/ΣS to give the S content at sulfide and sulfate "
                "saturation (Wieser and Gleeson 2023, following Jugo et al. 2010).")
    scss_name = st.selectbox("SCSS model", list(SCSS_MODELS), key="tot_scss")
    scas_name = st.selectbox("SCAS model", list(SCAS_MODELS), key="tot_scas")
    st.markdown("**Conditions**")
    T, P = ask("T", cols, "tot_T"), ask("P", cols, "tot_P")
    H2O = ask("H2O", cols, "tot_H2O", help="Overrides any H2O_Liq column.")
    Fe3 = ask("Fe3", cols, "tot_Fe3", help="Used by the SCSS model (only Fe²⁺ enters it).")
    spec = SCSS_MODELS[scss_name]
    sulf = None
    if spec["sulf"]:
        st.markdown("**Sulfide**")
        sulf = sulfide_form(cols, "tot_sulf")
    st.markdown("**S⁶⁺/ΣS**")
    s6_label, s6, extra = s6_form(cols, "tot", allow=("value",), kleinsasser=True)

    def job(df):
        n = len(df)
        t, p = series(resolve(T, df), n), series(resolve(P, df), n)
        h2o = resolve(H2O, df)
        kw = dict(df=df, T_K=t, P_kbar=p, H2O_Liq=h2o if spec["h2o"] else None,
                  Fe3Fet_Liq=resolve(Fe3, df) if scss_name != "Fortin et al. (2015)" else None)
        if sulf:
            kw.update(sulf(df))
        scss, m1, pr1 = call(spec["fn"], **kw)
        scas, m2, pr2 = call(SCAS_MODELS[scas_name]["fn"], df=df, T_K=t, H2O_Liq=h2o)
        SCSS, SCAS = scss[spec["out"]].reset_index(drop=True), scas["SCAS6_ppm"].reset_index(drop=True)
        if extra.get("model") == "Kleinsasser":
            tot, m3, pr3 = call(ss.calculate_S_Total_SCSS_SCAS, SCSS=SCSS, SCAS=SCAS, model="Kleinsasser",
                                deltaQFM=series(resolve(extra["dQFM"], df), n))
        else:
            tot, m3, pr3 = call(ss.calculate_S_Total_SCSS_SCAS, SCSS=SCSS, SCAS=SCAS, S6St_Liq=s6(df, t, p))
        out = pd.concat([tot.reset_index(drop=True), df.reset_index(drop=True)], axis=1)
        key = [c for c in ("Total_S_ppm", "S2_Tot_ppm", "S6_Tot_ppm", "S6St_Liq", "SCSS_2_ppm", "SCAS_6_ppm",
                           "SCSS_Tot", "SCAS_Tot") if c in out]
        return out, key, "Total_S_ppm", m1 + m2 + m3, pr1 + pr2 + pr3
    return job, f"SCSS: {scss_name}; SCAS: {scas_name}; S⁶⁺/ΣS: {s6_label}"


def form_s6(cols):
    st.markdown("**Model and conditions**")
    s6_label, s6, extra = s6_form(cols, "s6", allow=())
    needs_T = not s6_label.startswith("Jugo")
    T = ask("T", cols, "s6_T") if needs_T else 1473.15
    P = ask("P", cols, "s6_P") if s6_label.startswith("Gorojovsky") else 1.0

    def job(df):
        n = len(df)
        v = s6(df, series(resolve(T, df), n), series(resolve(P, df), n))
        out = df.reset_index(drop=True).copy()
        out.insert(0, "S6St_Liq", v.values)
        out.insert(1, "S2St_Liq", 1 - v.values)
        key = ["S6St_Liq", "S2St_Liq"]
        for i, (c, vals) in enumerate(extra.get("used", {}).items()):  # the redox state used, e.g. ΔQFM from fO2
            out.insert(2 + i, c, vals.values)
            key.append(c)
        return out, key, "S6St_Liq", [], ""
    return job, s6_label


def form_convert(cols):
    units = {"S (ppm)": "S_ppm", "S (wt%)": "S_wt", "SO₃ (wt%)": "SO3_wt", "SO₃ (ppm)": "SO3_ppm",
             "SO₂ (wt%)": "SO2_wt", "SO₂ (ppm)": "SO2_ppm", "SO₄ (wt%)": "SO4_wt", "SO₄ (ppm)": "SO4_ppm"}
    unit = st.selectbox("Your S is given as", list(units), key="cv_unit")
    if not cols:
        st.caption("Upload a file with a column of S contents.")
        return None, unit
    guess = next((c for c in cols if "S" in c and ("ppm" in c or "wt" in c or "SO3" in c)), cols[0])
    col = st.selectbox("Column", cols, index=cols.index(guess), key="cv_col")

    def job(df):
        out, msgs, printed = call(ss.convert_S_types, **{units[unit]: pd.to_numeric(df[col], errors="coerce").reset_index(drop=True)})
        out = out.reset_index(drop=True)
        return out, list(out.columns), None, msgs, printed
    return job, f"{col} as {unit}"


FORMS = {CALCS[0]: form_scss, CALCS[1]: form_scas, CALCS[2]: form_total, CALCS[3]: form_s6, CALCS[5]: form_convert}
ABOUT = {
    CALCS[0]: "Sulfur content at sulfide saturation (SCSS, as S²⁻). With an S⁶⁺ correction, SCSS_Tot_ppm is the "
              "total S (S²⁻ + S⁶⁺) at sulfide saturation.",
    CALCS[1]: "Sulfur content at anhydrite (sulfate) saturation (SCAS, as S⁶⁺). With an S²⁻ correction, SCAS_Tot_ppm "
              "is the total S at sulfate saturation.",
    CALCS[2]: "Total S, S²⁻ and S⁶⁺ at saturation, taking the lower of the sulfide and sulfate limits at the S⁶⁺/ΣS "
              "you choose.",
    CALCS[3]: "Proportion of S present as S⁶⁺, from Fe³⁺/ΣFe, fO₂ or ΔQFM.",
    CALCS[4]: "Your melt compositions plotted among the experiments each model was calibrated on.",
    CALCS[5]: "Converts S contents between S, SO₂, SO₃ and SO₄, in ppm or wt%.",
}


# %% ----------------------------------------------------------------
# results

@st.cache_data(show_spinner=False, max_entries=4)
def to_excel(results, df_in):
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as xw:
        results.to_excel(xw, sheet_name="Results", index=False)
        df_in.to_excel(xw, sheet_name="Input", index=False)
    return buf.getvalue()


def show_results(res, key, plot_col, msgs, printed, df_in, label):
    for m in msgs:
        st.warning(m)
    res = res.loc[:, ~res.columns.duplicated()].reset_index(drop=True)
    sample = sample_col(df_in)
    if sample is not None and sample not in res:
        res.insert(0, sample, df_in[sample].values)
    key = ([sample] if sample else []) + [c for c in key if c in res and c != sample]

    st.header("Summary")
    st.caption(label)
    c1, c2, _ = st.columns([1, 1, 2])
    c1.metric("Analyses", f"{len(res):,}")
    if plot_col and plot_col in res:
        v = pd.to_numeric(res[plot_col], errors="coerce")
        if v.notna().any():
            fmt = ".3g" if v.abs().median() < 10 else ".0f"
            c2.metric(f"Median {plot_col}", format(v.median(), fmt),
                      help=f"Range {format(v.min(), fmt)} to {format(v.max(), fmt)}")
            if plot_col == "SCAS_Tot_ppm" and v.median() > 1e5:
                st.info("SCAS_Tot is very large because S⁶⁺/ΣS is very low: SCAS (S⁶⁺ only) divided by a small S⁶⁺ "
                        "fraction. At these redox conditions the melt will not reach anhydrite saturation; sulfide "
                        "saturation limits the S content (see **Total S at sulfide and sulfate saturation**).")
            xs = [c for c in ("MgO_Liq", "FeOt_Liq", "SiO2_Liq", "T_K", "Fraction_melt") if c in res] + \
                 [c for c in numeric_columns(res) if c not in ("MgO_Liq", "FeOt_Liq", "SiO2_Liq", "T_K", "Fraction_melt")]
            x = st.selectbox("Plot against", xs, key="plot_x") if xs else None
            if x:
                ys = [c for c in key if c in numeric_columns(res) and (c.endswith("ppm") or c.startswith(("SCSS", "SCAS", "Total")))]
                ys = ys or [plot_col]
                long = res[[x] + ys].melt(x, var_name="Quantity", value_name="value")
                base = alt.Chart(long).encode(
                    x=alt.X(f"{x}:Q", title=x.replace("_", " "), scale=alt.Scale(zero=False)),
                    y=alt.Y("value:Q", title=("S (ppm)" if all(y.endswith(("ppm", "Tot")) for y in ys) else plot_col),
                            scale=alt.Scale(zero=False)),
                    color=alt.Color("Quantity:N", legend=alt.Legend(title=None, orient="top")),
                    tooltip=[alt.Tooltip(f"{x}:Q", format=".3g"), "Quantity", alt.Tooltip("value:Q", format=".4g")])
                chart = base.mark_line(point=True) if st.checkbox("Join points", value="Fraction_melt" in res,
                                                                  key="plot_join") else base.mark_circle(size=60)
                st.altair_chart(chart.properties(height=360).interactive(), width="stretch")
        n_bad = int(v.isna().sum())
        if n_bad:
            st.info(f"{n_bad:,} analyses returned no value (NaN): usually a missing input, such as an oxide, "
                    "temperature or redox value.")

    st.header("Results")
    fmt = {c: st.column_config.NumberColumn(format="%.4g") for c in key if pd.api.types.is_float_dtype(res[c])}
    st.dataframe(res[key], column_config=fmt, hide_index=True)
    d1, d2, _ = st.columns([1, 1, 3])
    d1.download_button("Download Excel", lambda: to_excel(res, df_in), file_name="PySulfSat_results.xlsx",
                       mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                       help="A Results sheet with every calculated column, plus your input.")
    d2.download_button("Download CSV", res.to_csv(index=False).encode(), file_name="PySulfSat_results.csv",
                       mime="text/csv")
    with st.expander("All columns"):
        st.dataframe(res, hide_index=True)
    printed = sorted({ln.strip() for ln in printed.splitlines() if ln.strip()})
    if printed:
        with st.expander("Messages from PySulfSat"):
            st.code("\n".join(printed), language=None)


def cali_page(df):
    """Your compositions among a model's calibration experiments (PySulfSat's return_cali_dataset)."""
    models = {f"SCSS: {k}": v["cali"] for k, v in SCSS_MODELS.items() if v.get("cali") and "eq. 12" not in k
              and v["cali"] != "O2021_SCSS"}  # O'Neill (2021): dataset not yet in PySulfSat
    models.update({f"SCAS: {k}": v["cali"] for k, v in SCAS_MODELS.items()})
    models.update({"S⁶⁺: O'Neill and Mavrogenes (2022), CS6": "O2022_CS6",
                   "S⁶⁺: Boulliung and Wood (2022), CS6": "BW2022_CS6"})
    with st.sidebar:
        st.header("Settings")
        label = st.selectbox("Calibration dataset", [k.replace(", eq. 11", "") for k in models], key="cali_model")
        model = {k.replace(", eq. 11", ""): v for k, v in models.items()}[label]
    cali, _, _ = call(ss.return_cali_dataset, model=model)
    numeric = lambda d: [c for c in d.columns if pd.api.types.is_numeric_dtype(d[c]) and not c.startswith("Unnamed")
                         and d[c].notna().any()]
    shared = [c for c in numeric(cali) if df is None or c in numeric(df)]
    if not shared:
        st.warning("Your data and this calibration dataset share no numeric columns to plot.")
        st.stop()
    order = ["SiO2_Liq", "MgO_Liq", "FeOt_Liq", "Na2O_Liq", "K2O_Liq", "CaO_Liq", "Al2O3_Liq", "T_K", "P_kbar", "H2O_Liq"]
    shared = [c for c in order if c in shared] + [c for c in shared if c not in order]
    with st.sidebar:
        x = st.selectbox("x axis", shared, index=0, key="cali_x")
        y = st.selectbox("y axis", shared, index=min(2, len(shared) - 1), key="cali_y")
    parts = [pd.DataFrame({"x": cali[x], "y": cali[y], "Dataset": f"{model} calibration",
                           "Label": cali["Sample"].astype(str) if "Sample" in cali else ""})]
    if df is not None:
        s = sample_col(df)
        parts.append(pd.DataFrame({"x": df[x], "y": df[y], "Dataset": "Your data",
                                   "Label": df[s].astype(str) if s else ""}))
    long = pd.concat(parts).dropna(subset=["x", "y"])
    dom = [f"{model} calibration"] + (["Your data"] if df is not None else [])
    chart = alt.Chart(long).mark_point(filled=True, size=70, opacity=0.75).encode(
        x=alt.X("x:Q", title=x.replace("_", " "), scale=alt.Scale(zero=False)),
        y=alt.Y("y:Q", title=y.replace("_", " "), scale=alt.Scale(zero=False)),
        color=alt.Color("Dataset:N", scale=alt.Scale(domain=dom, range=["#8c8c8c", "#d62728"][:len(dom)]),
                        legend=alt.Legend(title=None, orient="top")),
        shape=alt.Shape("Dataset:N", scale=alt.Scale(domain=dom, range=["circle", "triangle-up"][:len(dom)]),
                        legend=alt.Legend(title=None, orient="top")),
        tooltip=["Dataset", "Label", alt.Tooltip("x:Q", title=x, format=".4g"), alt.Tooltip("y:Q", title=y, format=".4g")])
    st.header("Calibration dataset")
    if df is not None:
        u = df[[x, y]].apply(pd.to_numeric, errors="coerce").dropna()
        inside = u[x].between(cali[x].min(), cali[x].max()) & u[y].between(cali[y].min(), cali[y].max())
        m1, m2, m3 = st.columns(3)
        m1.metric("Your analyses", f"{len(u):,}")
        m2.metric("Calibration experiments", f"{cali[[x, y]].dropna().shape[0]:,}")
        if len(u):
            m3.metric("Inside calibration range", f"{inside.mean():.0%}",
                      help="Share of your analyses whose x and y both fall within the calibration minimum to maximum.")
    st.altair_chart(chart.properties(height=460).interactive(), width="stretch")
    st.download_button("Download calibration dataset (CSV)", cali.to_csv(index=False).encode(),
                       file_name=f"PySulfSat_{model}_calibration.csv", mime="text/csv")
    with st.expander("Calibration dataset"):
        st.dataframe(cali, hide_index=True)


# %% ----------------------------------------------------------------
# page

with st.sidebar:
    calc = st.selectbox("Calculation", CALCS)
    st.header("Input data")
    src = st.radio("Input", ["Upload a file", "Example data"], horizontal=True, label_visibility="collapsed")
    data = name = sheet = fmt = None
    if src == "Example data":
        ex = st.selectbox("Example", list(EXAMPLE_FILES))
        fmt, rel, sheet = EXAMPLE_FILES[ex]
        path = EXAMPLES / rel
        if path.exists():
            data, name = path.read_bytes(), path.name
        else:
            st.error("The example files are not in this installation.")
    else:
        fmt = st.selectbox("File type", ["Spreadsheet", "Petrolog3 output", "MELTS .tbl file"],
                           help="Spreadsheet: one melt per row, oxides in wt% as SiO2_Liq or SiO2. "
                                "Petrolog3: the Excel output. MELTS: a liquid .tbl file.")
        types = {"Spreadsheet": ["csv", "xlsx", "xls"], "Petrolog3 output": ["xlsx", "xls"],
                 "MELTS .tbl file": ["tbl", "txt"]}[fmt]
        def _upload_changed():  # only fires when the user adds, replaces or removes the file
            f = st.session_state.get("upload_file")
            st.session_state["upload_kept"] = (f.getvalue(), f.name) if f is not None else None

        up = st.file_uploader("File", type=types, key="upload_file", on_change=_upload_changed)
        if up is not None:
            st.session_state["upload_kept"] = (up.getvalue(), up.name)
        kept = st.session_state.get("upload_kept")
        if kept is not None and Path(kept[1]).suffix.lower().lstrip(".") in types:
            # Keep using the file for the rest of the session even if the uploader comes back empty on a rerun
            data, name = kept
            if fmt == "Spreadsheet" and not name.lower().endswith(".csv"):  # Petrolog3: PySulfSat reads its own sheet
                sheets = pd.ExcelFile(io.BytesIO(data)).sheet_names
                sheet = st.selectbox("Sheet", sheets) if len(sheets) > 1 else sheets[0]

df_in, load_error = None, None
if data is not None:
    try:
        df_in = load_data(fmt, data, name, sheet)
    except Exception as e:
        load_error = f"Could not read {name}: {type(e).__name__}: {e}"
cols = [] if df_in is None else numeric_columns(df_in)

if calc != CALCS[4]:
    with st.sidebar:
        st.header("Settings")
        job, label = FORMS[calc](cols)

with st.sidebar:
    st.header("About")
    st.markdown("[Documentation](https://pysulfsat.readthedocs.io) · [GitHub](https://github.com/PennyWieser/PySulfSat) · "
                "[Paper](https://doi.org/10.30909/vol.06.01.107127)")
    st.caption(f"PySulfSat v{ss.__version__}, by Penny Wieser and coauthors. Questions or bugs? "
               "Open a [GitHub issue](https://github.com/PennyWieser/PySulfSat/issues).")

st.title("PySulfSat")
st.markdown("Sulfide and sulfate saturation, S⁶⁺/ΣS and total S in silicate melts, using published models. "
            "Choose the calculation, data and settings in the menu at left.")
st.info("**Please cite PySulfSat** (Wieser and Gleeson, 2023) and the paper for each model you use. "
        "See **How to cite** below.")
with st.expander("How to cite"):
    st.code(CITATION, language=None, wrap_lines=True)
    st.markdown("Also cite each model you use, e.g. Smythe et al. (2017) for the SCSS, Chowdhury and Dasgupta (2019) "
                "for the SCAS. The [documentation](https://pysulfsat.readthedocs.io) gives the full references.")
with st.expander("Input format"):
    st.markdown(
        "* **Spreadsheet**: one melt per row, oxides in wt% as `SiO2_Liq`, `FeOt_Liq`... or `SiO2`, `FeOt`... "
        "Iron as total FeO (`FeOt`). Optional columns: `H2O_Liq`, `Fe3Fet_Liq`, `Ni_Liq_ppm`, `Cu_Liq_ppm`, "
        "`T_K` or `T_C`, `P_kbar`, `Sample_ID`.\n"
        "* **Petrolog3** and **MELTS .tbl** files are read as exported, including temperature, pressure, H₂O and "
        "Fe³⁺/ΣFe.\n"
        "* Temperature, pressure and the other conditions can be a single value or a column of your file, chosen "
        "under **Settings**.")
st.markdown(f"**{calc}**: {ABOUT[calc]}")

if load_error:
    st.error(load_error)
    st.stop()
if calc == CALCS[4]:
    cali_page(df_in)
    st.stop()
if df_in is None:
    st.info("Upload a file in the menu at left, or choose **Example data**.")
    st.stop()
if df_in.empty:
    st.error("The file has no rows.")
    st.stop()
if len(df_in) > MAX_ROWS:
    st.error(f"{len(df_in):,} rows exceeds the {MAX_ROWS:,}-row web limit. Use `pip install PySulfSat` for larger datasets.")
    st.stop()
if job is None:
    st.stop()
try:
    with st.spinner(f"Calculating {len(df_in):,} rows..."):
        res, key, plot_col, msgs, printed = job(df_in)
except Exception as e:
    st.error(f"PySulfSat could not run this calculation: {type(e).__name__}: {e}")
    st.stop()
show_results(res, key, plot_col, msgs, printed, df_in, label)
