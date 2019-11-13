# ACCP

* Script plotWRF_2.py creates basic reflectivity, Doppler and scattering property fields from WRF output (Goddard 4ICE schemes).
* Script simAnalysis.py reads the output of the previous step, selects profiles with PIA(Ka)>10 and calculates reflectivity, Doppler and related variables with
instrument characteristics.
* Features and target variables are selected by script hbAnalysisClean.py
* Training and evaluation are carried out by script nnTest.py
