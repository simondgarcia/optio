# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 17:53:36 2022

@author: Simon D
"""

from demandlib import bdew
import pandas as pd
import matplotlib.pyplot as plt


def get_electric_profile(business_annual_demand: float,business='general'):
    '''
    Generate hourly electric profile of a commercial building (BETA)

    Parameters
    ----------
    business_annual_demand : float
        Total annual demand in kWh.
    business : str, optional
        Select the type of business. The default is 'general'.

    Returns
    -------
    dem_hourly : DataFrame
        Return dataframe of hourly electric load in kWh.

    '''
    business_type= {
        'general':'g0',
        'office':'g1',
        'bank':'g1',
        'school':'g1',
        'store':'g2',
        'supermarket':'g3',
        'hotel':'g6',
        'continous':'g3',
        'shop':'g4',
        'bakery':'g5',
        'club':'g6',
        'restaurant':'g6'}
    
    annual_demand={business_type[business.lower()]:business_annual_demand}

   
        #peak_demand*8760 #Assuming constant load throughout the year
    
    store_profile=bdew.ElecSlp(2019)
    dem=store_profile.get_profile(annual_demand,False)
    #Resample the data to have values for each hour
    dem_hourly = dem.resample("H").mean()
    #Convert index into DateTime
    dem_hourly.index = pd.to_datetime(dem_hourly.index)
    #Rename the Array
    dem_hourly.rename(columns={business_type[business]:business.capitalize()+' Demand'},inplace=True)
    
    #Average each hour into month 
    dem_hourly['Time']=pd.date_range(start='2019-01-01 00:00', periods=8760, freq='1H')
    dem_daily_profiles = dem_hourly.groupby([dem_hourly['Time'].dt.month, dem_hourly['Time'].dt.hour]).mean()
    dem_hourly['Hour']=dem_hourly.index.hour
    return dem_hourly


def plot_demand(dem: pd.DataFrame,start_date,end_date):
    '''
    Plot the demand as a series or DataFrame vector over the estipulated dates

    Parameters
    ----------
    dem : pd.DataFrame
        demand DataFrame or vector from function get_electric_profile or get_heating_profile.
    start_date : str
        Start date of the plot. Example: '2019-01-01 00:00' or '2019-01-01' must be in 2019.
    end_date : TYPE
        End date of the plot. Example: '2019-01-01 00:00' or '2019-01-01' must be in 2019..

    Returns
    -------
    Matplot plot of Energy [kWh] against time
    PNG image saved as "Timeseries of Demand [kWh].png".

    '''
    #New DF for plotting
    plot_df=dem.loc[(dem.index >= start_date) & (dem.index <= end_date)]

    # Change the style of plot
    plt.style.use('seaborn-darkgrid')
     
    # Create a color palette
    palette = plt.get_cmap('Set1')

    # Plot Hourly demand profile
    plt.plot(plot_df.index, plot_df[plot_df.columns[0]], marker='', color=palette(3), linewidth=1, alpha=0.9, label='Demand')
    # Add legend
    plt.legend(loc=2, ncol=2)
     
    # Add titles
    plt.title("Timeseries of Demand [kWh] ("+start_date+" to "+end_date+") ", loc='left', fontsize=12, fontweight=0, color='orange')
    plt.xlabel("Time")
    plt.ylabel("Energy [kWh]")
    plt.plot(dpi=500)
    plt.savefig("Timeseries of Demand [kWh].png",dpi=500)
