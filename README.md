# Welcome to OPTIO

Optio is a python-based tool that advice you on the optimal way to supply energy to commercial buildings using Solar PV systems coupled with batteries, maximising NPV, minimising LCOE, or maximsing avoided GHG emissions.

***
OPTIO is divided as follows:
1. Energy Demand Module: demand
2. Weather Module: get_weather
3. Grid module: grid_parameters
4. PV and battery optimisation module: OPTIO

***
The following image explains the interlink between modules and inputs of the optimization module OPTIO: !(optio block diagram.png)

#### demand modules:
1. **get_electric_profile:** take business type and annual energy demand as input to calculate an hourly Pandas Series of 8760 length of the energy demand

#### get_weather module has 3 main functions:
 1. **postcode_coordinates:** which takes a string of a postcode or address and makes an API call to OpenStreet to retrieve its coordinates as a 2 parameter list.
 2. **get_hourly_g_eff:** which is the main fuction to retrieve horuly irradiance on plane for the coordinates as input. Optional arguments include start and end date, PV technology, tilt, azimuth, and irradiance database. Please refer to [PVGIS](https://re.jrc.ec.europa.eu/pvg_tools/en/) documentation for more information
 2. **get_weather_data:** Which take a list containing lat and lon as coordinates of the location and returns a Pandas Dataframe containing hourly data of ambient temperature, realtive humidity, GHi, DNI, DHI, Infrared, Wind Speed, Wind Direction, Air pressure. This is an API call to TMY section of PVGIS [TMY](https://joint-research-centre.ec.europa.eu/pvgis-photovoltaic-geographical-information-system/getting-started-pvgis/api-non-interactive-service_en)
 
#### grid_parameters module has two main functions:
1. **get_grid_tariff** in which user can choose between constant and economy7. Tariffs have to be manually overwriten
2. **get_carbon_intensity:** this defines a 20 year list of annual carbon intensity factors for the UK grid based on UK API Carbon Intensity developed by National Grid


#### OPTIO:
This tool ask for user parameters about the project characteristics (business type), annual energy demand, location (preferibly in the UK), and technology characteristics and economics to determine the optimal size and operation of heat and electricity technologies.


## Let's Start
