==============================
Introduction and Citation
==============================

Welcome to PySulfSat- An Open-Source python3 tool for calculating sulfide and sulfate solubility in magmas.

This tool is published in Volcanica - https://doi.org/10.30909/vol.06.01.107127

Please make sure you cite PySulfSat if you use it, along with any of the underlying models you used.
e.g. 'We calculate the SCSS using the model of Smythe et al. (2017) implemented in the open-source Python3 tool PySulfSat (Wieser and Gleeson, 2023)'

==============================
Units
==============================

PySulfSat performs all calculations using  Kelvin for temperature and kbar for Pressure. Major element contents should be in wt%, all other units can be found in the doc strings of the relevant functions


==============================
Reporting bugs/issues with the code
==============================
No software is free of bugs, particularly when new features are being constantly added. We have extensively benchmarked PySulfSat to existing spreadsheets, and before the package is published on PyPI, automatic unit tests are run through GitHub in the attempt to catch problems introduced by changing Python dependencies/updates. However, if users spot any bugs, or wish to request features, they should submit an 'issue' on the GitHub page. Alternatively, they can email the author team. In both cases, please upload any files the code is using (e.g. excel, jupyter notebooks) so that I can run your code to see what the issue is!





