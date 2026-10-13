# PySulfSat web app

Browser front end for PySulfSat: upload melt compositions (a spreadsheet, a Petrolog3 output or a MELTS .tbl file)
or use the example data, choose a calculation and its settings, and get results with a plot and Excel/CSV downloads.
No Python needed for users.

| Calculation | Models |
|---|---|
| Sulfide saturation (SCSS) | Smythe et al. (2017), O'Neill (2021), O'Neill and Mavrogenes (2022), Li and Zhang (2022), Liu et al. (2021), Blanchard et al. (2021) eq. 11 and 12, Fortin et al. (2015). Sulfide composition fixed, calculated (O'Neill 2021 or Smythe et al. 2017, from melt Ni and Cu) or measured. Optional S6+ correction. |
| Sulfate saturation (SCAS) | Chowdhury and Dasgupta (2019), Zajacz and Tsay (2019), Masotta and Keppler (2015). Optional S2- correction. |
| Total S at sulfide and sulfate saturation | Any SCSS and SCAS model above, with S6+/ST from any model below or Kleinsasser et al. (2022). |
| S6+/ST | Nash et al. (2019), Jugo et al. (2010), O'Neill and Mavrogenes (2022), Boulliung and Wood (2022), Gorojovsky and Wood (2026) at 1 atm or with three pressure corrections. |
| Calibration datasets | Your compositions among each model's calibration experiments. |
| Convert S units | S, SO2, SO3, SO4 in ppm or wt%. |

Temperature, pressure, H2O, Fe3+/FeT, fO2 and the other inputs can each be a single value or a column of the file.
With the example Petrolog3 path, every model gives the same values as calling PySulfSat directly, as in the
example notebooks.

## Run locally

```
pip install -r WebApp/requirements.txt
streamlit run WebApp/pysulfsat_app.py
```

## Deploy (free)

1. Push to GitHub.
2. Go to https://share.streamlit.io, sign in with GitHub, and create an app.
3. Repository `PennyWieser/PySulfSat`, branch `main`, main file path `WebApp/pysulfsat_app.py` (case-sensitive),
   and set the App URL (e.g. `pysulfsat`) so it gets a readable address.

The example data come from `docs/Examples` in this repository, so deploy the whole repository. The app installs
PySulfSat from PyPI (pinned in `requirements.txt`); after a new release, bump the pin and push. Without Thermobar
installed, the Gorojovsky and Wood (2026) options are hidden.

## Not included yet

Monte Carlo error propagation, mantle melting, S isotope fractionation, and sulfide mass balance along a
fractionation path. MELTS .txt files are not offered: PySulfSat's `import_data(MELTS_txt=True)` uses
`pd.read_csv(delim_whitespace=True)`, which pandas has removed.
