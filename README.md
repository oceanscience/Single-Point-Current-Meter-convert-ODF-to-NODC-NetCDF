# Single-Point-Current-Meter-convert-ODF-to-NODC-NetCDF
Converts the ODF (Ocean Data Format) for moored Single Point Current Meter time series data to a NetCDF format that is compliant with IOOS, ACDD and CF standards. The ODF format is an historic format being used in the Department of Fisheries and Oceans Canada. Code is written in R.

INPUT: ODF files of moored moored Single Point Current Meter (MCM) data, GF3 Code Map_final_v2.csv, Shared_models_vessel_BIO.csv, OUTPUT: writes csv files saved in the processed folder; mooringprojects.csv, mooringprojectsh.csv, Shared_models_vessel_BIO.csv, a map of mooring locations and NetCDF file

Folder structure: working directory with code and some inputs data folder with ODF files and output

Code:
*  1ODFmctc_title_program_project.R - run this first to create mooringprojectsh.csv file that lists certain metadata. edit by hand before creating NetCDF 
*  2ODFmctd_to_netcdf2.R - run this next to create NetCDF.
**  ODFmctd_ships.R - this can be run if needed to list the ship names and instrument type, model and serial numbers to ensure the ships are listed in Shared_models_vessel_BIO.csv and there is correct instrument information

Notes:

ensure ship used for mission is included in file Shared_models_vessel_BIO.csv
edit the Shared_models_vessel_BIO.csv to include program, project and mission information
ensure conductivity has correct units and conductivity ration is converted to conductivity
ensure temperature scale is T90
ensure correct instrument model and serial numbers, and correct BODC code is used. refer to MicrocatSBE_models.txt file
creator_name: is who created the NetCDF
publisher_name: is Fisheries and Oceans Canada (DFO)

attached packages: oce, ocedata, stringr, glue, ncdf4, lubridate, dplyr, anytime, 
