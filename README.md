# rf-envmod
Random forest models for predicting stream condition
A collection of R scripts for predicting local environmental conditions from StreamCat
variables. 

dataprep.R : downloads relevant predictor variables from StreamCat and interpolates air temperature and precipitation
from PRISM data. Site information that is required is location (lat and lon) and the associated COMID from NHDPlus. 
These data are loaded from the tab-delimited file example.dat.txt

envmod.R : Prediction local environment measurements for conductivity, dissolved phosphorus, substrate sand/fines, and 
total suspended solids for the sites specified in example.dat.txt. Should run without modification after 
running dataprep.R.
