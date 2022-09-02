# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 17:36:46 2022

@author: Simon D
"""

import pandas as pd
import requests
import matplotlib.pyplot as plt

# anvil.server.connect("server_O422MKDZYMEL65AZNAJPPAX6-O4NDYNEN6BTP2ZMD")

#----------------------------------------------
#I. LOCATION Geocoding and mapping
#----------------------------------------------
#   1. Get location coordinates
def postcode_coordinates(postcode):
    url_map = 'https://nominatim.openstreetmap.org/search/' + postcode +'?format=json'
    map_res = requests.get(url_map).json()
    coordinates=[map_res[0]['lat'],map_res[0]['lon']]
    return coordinates

#------------------------------------------------------------------------------------
#II. DATA PROVIDER: Choose data input provider PVGIS or customize Excel file
#------------------------------------------------------------------------------------
def get_weather_data(coordinates,data_provider='PVGIS'):
 
    #OPTION 1: API request call to PVGIS to download JSON info for TMY Weather and Irradiance
    if data_provider == 'PVGIS':
            
        # 1. Go to PVGIS and with location input by user
        url_initial='https://re.jrc.ec.europa.eu/api/v5_2/tmy?'
        url_params=[url_initial,'lat=',str(coordinates[0]),'&lon=',str(coordinates[1]),'&startyear=2005','&endyear=2020','&outputformat=json','&browser=0']
        url_weather="".join(url_params)
        
        # 2. Request call to url to get data in JSON
        weather_res = requests.get(url_weather)
        if weather_res:
            weatherdata=weather_res.json()
        
        # 3. Translate JSON to Dataframe and assign column names
        weather_df=pd.DataFrame.from_dict(weatherdata['outputs']['tmy_hourly'])
        weather_df.columns=["Time","Temp_amb","Rel_Humidity","GHI","DNI","DHI", "Infrared", "Wind_Speed","Wind_Direction","Air_Pressure"]
       
        # 4. Eliminate Time Column and Generate a new one to avoid having multiple years
        weather_df.drop('Time',axis=1)
        weather_df["Time"]=pd.date_range(start='2019-01-01 00:00', periods=8760, freq='1H')
        weather_resampled = weather_df.groupby([weather_df['Time'].dt.month, weather_df['Time'].dt.hour]).mean()
        #print(weather_df)
        return weather_df
    #-----------------------------------------------------------------------------------------
    #OPTION 2: Read excel/csv file with information
    elif data_provider == 'Excel':
    
        #Names of the columns to import
        
        names=["Time","Hour","Temp_amb","Rel_Humidity","GHI","DNI","DHI", "Infrared", "Wind_Speed","Wind_Direction","Air_Pressure","Complete Time"]
        table_weather = pd.DataFrame(columns = names)
        
        #Link to the file and datasheet name
        
        tmp=pd.read_excel("tmy_"+str(coordinates[0])+"_"+str(coordinates[1])+"_2005_2020.xlsx", sheet_name="tmy_51.546_-0.144_2005_2020",
                        names=names,
                        dtype={'Time': str,'Hour': str},
                        skiprows=16,
                        skipfooter=12
                         )
        
        # Merge date and hour and drop hour      
        
        table_weather = pd.concat([table_weather,tmp],ignore_index=True)
        table_weather['Time']=table_weather['Time']+" "+table_weather['Hour']
        
        # convert the column (it's a string) to datetime type
        my_datetime =pd.to_datetime(table_weather['Time'])
        table_weather['month']=my_datetime.dt.month
        table_weather['day']=my_datetime.dt.day
        table_weather['hour']=my_datetime.dt.hour
        table_weather.drop('Hour',axis=1,inplace=True)
        table_weather.drop('Time',axis=1,inplace=True)
        #print(month)
        #print(day)
        #print(hour)
        print(table_weather)
        
        #print(my_datetime)
        
        # create datetime index passing the datetime series
        datetime_index = pd.DatetimeIndex(my_datetime)
        table_weather2=table_weather.set_index(datetime_index)
        
        return table_weather
        #print(table_weather)
                    
    else:
        print('Please select PVGIS or Excel')
        #get_data='Excel'
def weather_df_to_list(df):
    
    year_tamb=[x for x in df.loc[:,"Temp_amb"]]
    year_wind=[x for x in df.loc[:,"Wind_Speed"]]
    
    return year_tamb , year_wind

def get_hourly_g_eff(coordinates,database='PVGIS-SARAH2',startyear='2019',endyear='2019',pv_tech='crystSi',tracking='0',slope='0',azimuth='0',optimize_slope='1',optimize_azimuth='1'):
    
    url_initial='https://re.jrc.ec.europa.eu/api/v5_2/seriescalc?'
    url_params=[url_initial,'lat=',str(coordinates[0]),'&lon=',str(coordinates[1]),'&startyear=',startyear,'&endyear=',endyear,'&pvtechchoice=',pv_tech,'&trackingtype=',tracking,'&angle=',slope,'&aspect=',azimuth,'&optimalinclination=',optimize_slope,'&optimalangles=',optimize_azimuth,'&components=1','&outputformat=json','&browser=0']
    url_irr="".join(url_params)
    # 2. Request call to url to get data in JSON
    irr_res = requests.get(url_irr)
    if irr_res:
        irr_data=irr_res.json()
    else:
        print('API call gone bad')
    
    # 3. Translate JSON to Dataframe and assign column names
    irr_df=pd.DataFrame.from_dict(irr_data['outputs']['hourly'])
    irr_df.columns=["Time","direct_poa","diffuse_poa","reflected_poa","sun_height","temp_amb","wind_speed", "reconstructed"]
   
    # 4. Eliminate Time, Wind Speed, Temperature and Generate a new date range to avoid having different years
    irr_df.drop('Time',axis=1,inplace=True)
    irr_df.drop('temp_amb',axis=1,inplace=True)
    irr_df.drop('wind_speed',axis=1,inplace=True)
    # irr_df["Time"]=pd.date_range(start='2019-01-01 00:00', end='2019-12-31 23:00', freq='1H')
    irr_df["total_poa"]=irr_df["diffuse_poa"]+irr_df["reflected_poa"]+irr_df["direct_poa"]
    irr_df.index=pd.date_range(start='2019-01-01 00:00', periods=8760, freq='1H')
    irr_df['Time']=pd.date_range(start='2019-01-01 00:00', periods=8760, freq='1H')
    
    #5. Changing from W/m2 to kW/m2
    irr_df.loc[:,"direct_poa"]=irr_df.loc[:,"direct_poa"]*0.001
    irr_df.loc[:,"diffuse_poa"]=irr_df.loc[:,"diffuse_poa"]*0.001
    irr_df.loc[:,"reflected_poa"]=irr_df.loc[:,"reflected_poa"]*0.001
    irr_df.loc[:,"total_poa"]=irr_df.loc[:,"total_poa"]*0.001
    #Average each hour into month 
    irr_resampled = irr_df.groupby([irr_df.index.month, irr_df.index.hour]).mean()
    
    return irr_df

def irr_df_to_list(df):
    g_eff= df["total_poa"]
    year_irr = [x for x in g_eff.loc[:]]
    return year_irr

#=============================================================================
def plot_results(df,start='2019-03-26 00:00',end='2019-03-27 00:00'):
    #Graphs---------------------------------------------------------------------------------------------------------------------------------------------
    #New DF for plotting
    plot_df=pd.DataFrame()
    plot_df['Time']=pd.date_range(start='2019-01-01 00:00', periods=8760, freq='1H')
    plot_df=df.loc[(df['Time'] >= start) & (df['Time'] < end)] #Start and end of input dataframe
    
    plot_df['Hour']=plot_df['Time'].dt.hour
    plot_df['Month']=plot_df['Time'].dt.month

    fig, ax1 = plt.subplots(figsize=(8, 8))

    # Change the style of plot
    plt.style.use('seaborn-darkgrid') #Seaborn darkgrid style for plot
     
    # Create a color palette
    palette = plt.get_cmap('Set1') #Set 1 palette of colors for plot plt
    
    ax1.set_ylabel("Irradiance [kW/m2]", fontsize=13)

    
    # # Plot multiple lines
    plt.plot(plot_df['Time'], plot_df['direct_poa'], marker='', color=palette(0), linewidth=0.5, alpha=0.9, label='Direct') #Plotting Variable 1
    plt.plot(plot_df['Time'],plot_df['diffuse_poa'], marker='', color=palette(1), linewidth=0.5, alpha=0.9, label='Diffuse') #Plotting Variable 2
    plt.plot(plot_df['Time'],plot_df['reflected_poa'], marker='', color=palette(2), linewidth=0.5, alpha=0.9, label='Albedo') #Plotting Variable 3
    
    # Add legend
    plt.legend(loc=2, ncol=2, fontsize=14)
     
    # Add titles
    title='Direct, Diffuse, and Reflected Irradiance for Central London'
    plt.title(title, loc='left', fontsize=16, fontweight=1, color='orange')
    plt.plot(dpi=500)
    #Save image as png
    plt.savefig("Direct, Diffuse, and Reflected Irradiance for Central London.png",dpi=500)