# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 21:55:06 2022

@author: Simon D
"""
import pandas as pd



# -----------------------------
# GRID PARAMATERS
#------------------------------

def get_grid_tariff(tariff_type,start='2019-01-01 00:00',end='2019-12-31 23:00'):
    '''
    Get the grid tariff either as a dual tariff day/night (economy7), or a constant

    Parameters
    ----------
    tariff_type : str
        Choose from 'constant',or 'economy7'.

    Returns
    -------
    grid_tariff

    '''
    # 1. Define a electricity price depending on the type chosen
    date=pd.date_range(start, end, freq='1H')
    grid_tariff=[]
    change=0.0 #% of increase or decrease in energy prices
    
    peak_tariff=0.30*(1+change) #£/kWh
    off_peak_tariff=0.25*(1+change) #£/kWh
    average_tariff=0.18 #£/kWh
    if tariff_type=='economy7':
        
        for i in range(0,8760):
            if int(date.hour[i]) < 6 or int(date.hour[i]) >21:
                grid_tariff.append(off_peak_tariff) #£/kWh
            else:
                grid_tariff.append(peak_tariff) #£/kWh
        #print('You selected Economy 7 grid tariff of 0.15 £/kWh during offpeak and 0.25 £/kWh during peak hours')
       
    else:
        grid_tariff=[average_tariff]*8760  #£/kWh

    return grid_tariff


def get_carbon_intensity(variable_carbon=False):
    '''
    Define a carbon intensity parameter: Constant throughout time or variable grabbed from API https://api.carbonintensity.org.uk/

    Parameters
    ----------
    variable_carbon : Boolean, optional
        Select true if carbon intensity is variable across time. The default is False.

    Returns
    -------
    20 year list with yearly average carbon intensity in gCO2/kWh for UK.

    '''

    if variable_carbon:
        carbon_intensity=[
            183.675,
            180.175,
            176.675,
            173.175,
            169.675,
            166.175,
            162.675,
            159.175,
            155.675,
            152.175,
            148.675,
            145.175,
            141.675,
            138.175,
            134.675,
            131.175,
            127.675,
            124.175,
            120.675,
            117.175
            ]    #gCO2/kWh
    else:
        carbon_intensity= [150]*20 #gCO2/kWh
    return carbon_intensity
    