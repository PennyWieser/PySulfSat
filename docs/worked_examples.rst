================
Worked examples
================

This page summarizes the different examples available. If you want a specific example, we can always include more.

Importing data
=====================

This :doc:`example <Examples\Data_Input\Importing_Spreadsheet>` shows how to load in data from a spreadsheet

This :doc:`example <Examples\Data_Input\Importing_Petrolog>` shows how to load in data from a Petrolog3 crystallization path

This :doc:`example <Examples\Data_Input\Importing_MeltsTBL>` shows how to load in data from different MELTS paths


Integrating PySulfSat and Petrolog3
==========================================

This :doc:`example <Examples\Intro_Example_Petrolog_FC\FC_Petrolog>` shows how to calculate the SCSS and SCAS for a petrolog3 fractionation path using a wide range of models.
It also shows how to perform corrections using different models for S6+/ST



Integrating PySulfSat and MELTS calculations
==========================================

This :doc:`example <Examples\Integrating_with_PetThermoTools\Single_FC_Model>` shows how to run a single MELTS crystallization model in PeTThermoTools
and then perform SCSS calculations. To run a FC calculation with more water, see this  :doc:`example <Examples\Integrating_with_PetThermoTools\Single_FC_Model_morewaterrich>`


This :doc:`example <Examples\Integrating_with_PetThermoTools\Polybaric_FC_Model>` shows how to run FC models at multiple pressures, and then calculate the SCSS


Modelling Mantle Melting
==========================

This :doc:`example <Examples\Mantle_Melting_Lee_Wieser\Simple_Melting_Cu_S_Ba>` shows how to model S, chalcophile and lithophile elements during mantle melting, following an adaptation of the model of Lee et al. (2012). We keep the S content of the melt fixed for simplicity. We show how to perform calculations at different S6/ST ratios

This :doc:`example <Examples\Mantle_Melting_Lee_Wieser\Complex_melt_changing_SCSS>` is a more complex example. It loads a Thermocalc melting path (e.g. major elements of mantle melts), and uses this to calculate a SCSS at each step in the model. It also accounts for changing silicate modes.

Calculating and correcting for S6
==========================
This :doc:`example <Examples\S6_S2_Corrections\S6_S2_Corrections_Nash_Jugo_Kleinsasser>` shows how to calculate SCSS2- and SCSSTot and SCAS6+ and SCASTot, and make the plot shown in the PySulfSat paper of total S vs. QFM. It also shows how to perfrom corrections for different s6+ models to a petrolog3 fractionation path.



This :doc:`example <Examples\S6_S2_Corrections\CalcS6ST_Muth>` shows how to calculate S6/St using different models, including how to propagate uncertainty using monte-carlo methods based on errors in input major element contents, for both S6/ST and SCSS

This :doc:`example <Examples\S6_S2_Corrections\CS6_S6ST_Correction>` shows how to calculate S6/St from lnC6, lnCS2 etc. using Oneill and Mavrogenes (2022) and Boulling and wood (2022) for comparison.

Icelandic case study
=======================
This This :doc:`example <Examples\Sulf_Evolution_During_FC\Sulfide_sat_magma_evolution_Iceland>` show how to model the SCSS following the method used in Liu et al. (2024) examining sulfide saturation at Holuhraun, Iceland. Calculated SCSS sulfide compositions are compared to measured compositions, and the amount of sulfide removed at each step is also calculated.


Other useful functions
==========================
This :doc:`example <Examples\Other_Useful_Functions\Calculating_KDs_Kiseeva>` shows how to use the KD model of Kiseeva and Wood (2015) and Brenan (2015).

This :doc:`example <Examples\Other_Useful_Functions\Converting_S_values>` shows how to convert between different S isotope notation, and convert SO3 and SO2 to S in ppm etc.

This :doc:`example <Examples\Other_Useful_Functions\Plotting_Cali_Datasets>` shows how to plot the calibration range of each model against your data.


This :doc:`example <Examples\S_isotope_Fractionation_Models\Frac_factors>` shows how to calculate sulfide-melt fractionation factors following Miyoshi et al. (1984), Fiege et al. (2015). We benchmark to Rezeau et al. (2023)



This :doc:`example <Examples/Liquid_Ol_Liq_Themometry/Olivine_Liquid_thermometry>` shows:

