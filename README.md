# Single-Point-Current-Meter-convert-ODF-to-NODC-NetCDF
Converts the ODF (Ocean Data Format) for moored Single Point Current Meter (MCM) time series data to a NetCDF format that is compliant with IOOS, ACDD and CF standards. The ODF format is an historic format being used in the Department of Fisheries and Oceans Canada. Code is written in R. The file MCM_netCDF_format.docx outlines the variable attributes and global attributes.

- INPUT: ODF files of moored moored Single Point Current Meter (MCM) data, GF3 Code Map_final_v2.csv, Shared_models_vessel_BIO.csv, 
- OUTPUT: writes csv files saved in the processed folder; mooringprojects.csv, mooringprojectsh.csv, mooringships.csv, a map of mooring locations mooring_map.jpeg and NetCDF file

Folder structure: working directory with code and some inputs data folder with ODF files and output

Code:
*  1_ODFmcm_title_program_project.R - run this first to create mooringprojectsh.csv file that lists certain metadata. edit by hand before creating NetCDF, see notes below
*  2_ODFmcm_ships.R - this can be run if needed to list the ship names, instrument type, model and serial numbers to ensure the ships are listed in Shared_models_vessel_BIO.csv and there is correct instrument information. Also creates a map to check mooring locations.
*  3_ODFmcm_to_netcdf5.R - run this next to create NetCDF.
   
Notes:

* ensure ship used for mission is included in file Shared_models_vessel_BIO.csv
* ensure conductivity has correct units and conductivity ration is converted to conductivity
* ensure temperature scale is T90
* ensure correct instrument model and serial numbers, and correct BODC code is used. refer to MicrocatSBE_models.txt file
* owner is Fisheries and Oceans Canada (DFO)
* originator is BIO	
* edit the mooringprojectsh.csv, see example file:
    * add columns: program, project, url and location_description
    * use the map created by the script ODFmcm_ships.R to determine the location,
    * the choices for program are Maritimes Region Historical offshore and Moored Time Series Data and Maritimes Region Historical coastal Moored Time Series Data, AZMP, AZOMP, RAPID  
    * The project and url can be left blank if there is nothing specified in the ODF
    * check the mission information

attached packages: oce, ocedata, stringr, glue, ncdf4, lubridate, dplyr, anytime, 
