# Welcome to OPTIO

Optio is a python-based tool that advises you on the optimal way to supply energy to commercial buildings using Solar PV systems coupled with batteries, maximising NPV, minimising LCOE, or maximising avoided GHG emissions.

***
OPTIO is divided as follows:
1. Energy Demand Module: demand
2. Weather Module: get_weather
3. Grid module: grid_parameters
4. PV and battery optimisation module: OPTIO

***
The following image explains the interlink between modules and inputs of the optimization module OPTIO: ![OPTIO Block Diagram](https://github.com/simondgarcia/optio/blob/main/optio%20block%20diagram.png)

#### demand modules:
1. **get_electric_profile:** take business type and annual energy demand as input to calculate an hourly Pandas Series of 8760 size of the energy demand

#### get_weather module has 3 main functions:
 1. **postcode_coordinates:** which takes a string of a postcode or address and makes an API call to OpenStreet to retrieve its coordinates as a 2-parameter list.
 2. **get_hourly_g_eff:** is the main function to retrieve hourly irradiance on plane data. It takes a list containing lat and lon s arguments for the coordinates and returns an 8760 hourly DataFrame containing the irradiance components and total. Optional arguments include start and end date, PV technology, tilt, azimuth, and irradiance database. Please refer to [PVGIS](https://re.jrc.ec.europa.eu/pvg_tools/en/) documentation for more information
 2. **get_weather_data:** Which takes a list containing lat and lon as coordinates for the location and returns a Pandas Dataframe containing hourly data of ambient temperature, relative humidity, GHi, DNI, DHI, Infrared, Wind Speed, Wind Direction, Air pressure. This is an API call to the TMY section of PVGIS [TMY](https://joint-research-centre.ec.europa.eu/pvgis-photovoltaic-geographical-information-system/getting-started-pvgis/api-non-interactive-service_en)
 
#### grid_parameters module has two main functions:
1. **get_grid_tariff**: this is a function that returns the hourly grid tariff based on peak and off-peak set-up values. The user inputs the tariff type as a string (constant or economy7). Tariffs have to be manually overwritten.
2. **get_carbon_intensity:** this function returns a 20-year list of annual carbon intensity factors for the UK grid based on UK API Carbon Intensity developed by National Grid


#### OPTIO:
This tool asks for user parameters about the project characteristics (business type), annual energy demand, location (preferably in the UK), and technical characteristics and economics to determine the optimal size and operation of electricity technologies.

##### OPTIO has four main functions:
1. **main:** the main code is where all the pre-defined variables are loaded. They can be overwritten by the user or modified. This includes location, annual demand, business type, PV costs, PV technical parameters, Emission factors, Battery costs, and Battery technical parameters, among others.
2. **run_model:** this function takes the loaded input data and performs a yearly calculation of the energy flows. It takes by argument PV capacity [kWp] and Battery capacity [kWh]. The return is a 6-size tuple containing (Lifecycle costs, results: a Dataframe with the financial and environmental results year by year, LCOE, output_df: a dataframe containing the hourly energy flows and weather info, NPV, total avoided emissions)
3. **iteration:** This function takes the 6 tuple result from run_mode and performs a grid search to generate a dataframe with the NPV, LCOE, and Avoided Emissions for each combination of PV and Battery set as lists at the beginning of the function. The output is a tuple size 3 containing the dataframes for NPV, LCOE, and avoided emissions respectively.
4. **results:** This function takes a 6 tuple size output from the run_model function and calculates key indicators of the configuration chosen. Return is a dictionary.


## Let's Start
