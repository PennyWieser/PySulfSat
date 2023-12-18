================================================
Change Log
================================================
Documentation update only (Dec 18th, 2023)
================================
PyMELTScalc renamed to PetThermoTools - updated examples.

Version 1.0.4 (Nov 11th, 2023)
================================
Updated calibration datasets to read from csvs to get around issues with pandas 2 pkl compatability

Version 1.0.3 (Sept 7th, 2023)
================================
Samer Mashhour pointed out different versions of MELTS give different files from the .tbl in the example, providing a
file named _tbl.txt.
To accomadate these, we have added an option for the import_data function to have MELTS_txt=True
Thanks Samer!

Version 1.0.1 (April 30th, 2023)
================================
Added Sulfide isotope fractionation models following https://doi.org/10.1016/j.chemgeo.2023.121325

