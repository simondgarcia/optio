# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 16:19:17 2022

@author: Simon D
"""

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import numpy as np
import demand as dem
import get_weather as w
import pvlib
import matplotlib.pyplot as plt
import grid_parameters as grid
import time

# DATA INPUTS------------------------------------------------------------

#-------LOCATION--------------------------------------------------
   
test_postcode= 'bn6 9pq' #Project's Postcode ar address (If there is an error, check the address in OpenStreet.org)
test_coord = w.postcode_coordinates(test_postcode) #Function calling to turn a postcode into coordinates lat,lon

#Project Lifetime ---------------------------------------------------

N_lifetime = 20 #years

#------------------------------------------------------------------------------------------------------------------------
# PV Technical parameters
coef_temp = -0.0035  # [1/degrees celsius] #Temperature Coeficient PV panel
area_m = 1.65  # m2 #Module area
power_refm = 300  # W #Module Nominal Power
eff_m = 0.20  # % #Module efficiency
pv_deg=0.008 #Yearly PV degradation
dc_losses=0.05 #DC losses (wiring, connections, mismatch, nameplate)
ac_losses=0.03 #Wiring and Inverter loss
shadow_losses=0.04 #Shadow and enviromental losses
dc_ac_ratio = (1-0.01) #DC-AC Ratio loss (assuming DC AC Ratio of 1.2)

# Battery Technical Parameters
c_rate=2 #Battery c-rate discharge time of 2 hours
dod=0.8 #Depth of discharge of 85%
alpha_b=1 #Availability factor
SOC_min=0.15 #Min State of Charge of 15%
SOC_max=SOC_min+dod #Max State of Charge 95%
vol_b=1/80 #Volumetric Energy density of Li-ion battery [m3/kWh]

#Project constant parameters
area_total=5000 #m2 #Total Roof area
area_max = 750 #m2 #Available Roof Area for PV
vol_max = 75 #m3 #Available Volume for Battery

#Financial parameters
discount_rate=0.09
discount_factor=[1/((1+discount_rate)**x) for x in range(0,20)]
capex_bat=[1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0] #Capex occur in year 1 and 10 for Battery system (10 year life)
capex_pv=[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]  #Capex only occur in first year for PV system (20 year life)

#Cost parameters PV
c_max = 110000  # £
c_pv_mod = 250*1 # £/kWp
c_pv_inv = 120  # £/kWp
c_pv_BOS = 200 # £/ kWp 
c_pv_install = 250*1  # £/kWp
c_pv_capex=c_pv_install+c_pv_mod+c_pv_inv+c_pv_BOS
c_pv_opex = 10*1  # £/kWp-yr
# Cost parameters Battery Assuming Li-ion
red_fact=-0.0
c_b_cap = 180*(1-red_fact) #£/kWh
c_b_pcs = 220*(1-red_fact) #£/kW
c_b_BOS = 80*(1-red_fact) #£/kW
c_b_install = 80*(1-red_fact) #£/kWh
c_b_capex_pcs=c_b_BOS+c_b_pcs
c_b_capex_cap=c_b_cap+c_b_install

c_b_opex = 8*(1-red_fact)#£/kWh-yr
#Cost Parameters Grid
energy_inc=[1.03**y for y in range(0,20)]
c_elec = grid.get_grid_tariff('economy7')#Grid import tariff as a hourly list
c_exp=0.06 #£/kWh export Smart Export Guarantee
carbon_price=50*(1/1000) #£/kgCO2

#Emissions parameters Source:https://www.nrel.gov/docs/fy21osti/80580.pdf
emi_pv_up=28 #gCO2/kWh
emi_pv_on=10 #gCO2/kWh
emi_pv_down=5 #gCO2/kWh
emi_b_up=73000 #gCO2/kWh
emi_b_on=0 #gCO2/kWh
emi_b_down=0 #gCO2/kWh
emi_bat=[1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0] #Emissions occur in year 1 and 10 as they are embodied battery capacity
carbon_intensity=grid.get_carbon_intensity(variable_carbon=True) #gCO2/kWh

# Variable parameters
wdf=w.get_weather_data(test_coord) #Creation on weather dataframe asking PVGIS to get hourly data for the location test_coord
year_tamb,year_wind=w.weather_df_to_list(wdf) #Creation of temp_amb and wind_speed yearly vectors as lists from Dataframe wdf

irr_df= w.get_hourly_g_eff(test_coord) #function calling to get hourly irradiance for location test_coord and with parameters as default
year_irr=irr_df["diffuse_poa"]+irr_df["reflected_poa"]+irr_df["direct_poa"] #year effective irradiance vector summing the three components in kW/m2
year_irr_w= year_irr.loc[:]*1000 #year irradiance in W/m2

year_tpanel= pd.Series(pvlib.temperature.faiman(year_irr_w, year_tamb, year_wind)) #Panel temp series creation by calling function in pvlib with temp_amb and wind_speed as parameters

business_type='store' #Type of business
business_annual_demand=90000 #Annual electrcity demand in kWh
dem_series=dem.get_electric_profile(business_annual_demand, business_type)[business_type.capitalize()+' Demand'] #Electric Load profile generation calling demandlib with store defined as business_type and with annual demand. Generates a hourly dataframe/vector


def run_model(PV_cap, battery_energy):
    '''
    

    Parameters
    ----------
    PV_cap : TYPE
        DESCRIPTION.
    battery_energy : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    # Grab Currrent Time Before Running the Code
    start = time.time()
    #Time range

 # years #Technology and project lifetime
    t=range(0,8760)
    y=range(0,N_lifetime)
    
    c_capex= []
    c_opex=[]
    c_grid=[]
    c_export=[]
    c_carbon=[]
    net_flow=[]
    discounted_flow=[]
    net_flow_dif=[]
    discounted_flow_dif=[]
    grid_emissions=[]
    battery_emissions=[]
    pv_emissions=[]
    total_emissions=[]
    baseline_emissions=[]
    avoid_emissions=[]

    #Decision variables
    PV_cap=PV_cap
    N_panels=PV_cap/(power_refm/1000) #kWp
    e_b_max=battery_energy#battery_energy #kWh
    
    #Baseline

    baseline_grid=[-sum(dem_series[h]*c_elec[h] for h in t) *energy_inc[y] for y in range(0,20)] #Demand[h]*Energy Tariff[h]*Energy Increase Index for every hour in a year for every year in project lifetime
    baseline_carbon=[]
    baseline=[]
    for y in y:
        
        #PV losses   
        pv_losses = 1-((1-dc_losses)*(1-ac_losses)*(1-(pv_deg*(y+1)))*(1-shadow_losses))# % #Module losses
        pv_eff_losses = (1-pv_losses)  # % #System efficiency
        
        #Battery losses
        rte=0.86 #Round trip efficiency of 89%
        b_deg=0.005 #Battery degradation of 0.5% yearly
        b_losses = b_deg*(y+1)# % #Module losses
        rte -= b_losses  # % #System efficiency
        
        #CONSTRAINTS---------------------------------------------------------------------------------------------
        
        #Power to Energy ratio for Battery based on C-rate of battery
        k_b=e_b_max/c_rate
        
        #PV DC output at time t
        e_dc = N_panels*area_m*eff_m*year_irr*(1+coef_temp*(year_tpanel-25))
        
        #AC PV production at time t
        e_pv = e_dc*(pv_eff_losses*dc_ac_ratio)
        
        if y==0:
            pv_out=e_pv
        
        #PV surplus
        pv_surplus=np.maximum(e_pv-dem_series.values,0)
        
        e_b_c=[] #Energy Battery Charging
        e_b_d=[] #Energy Battery Discharging
        e_imp=[] #Energy import
        e_exp=[] #Energy export
        SOC=[] #Battery State of Charge  
        
        #Dispatch Strategy
        for h in t:
            if h !=0:
                if pv_surplus[h]>0:
                    e_b_d.append(0)
                    e_imp.append(0)
                    if SOC[h-1] == SOC_max:
                        e_b_c.append(0)
                    else:
                        e_b_c.append(min(pv_surplus[h],(SOC_max-SOC[h-1])*alpha_b*e_b_max,k_b))
                    e_exp.append(pv_surplus[h]-e_b_c[h])
        
                else:
                    if (SOC[h-1]-SOC_min)*rte*alpha_b*e_b_max > (dem_series[h]-e_pv[h]):
                        e_b_d.append(min(k_b,dem_series[h]-e_pv[h]))
                    else:
                        e_b_d.append((SOC[h-1]-SOC_min)*rte*e_b_max*alpha_b)
                    e_b_c.append(0)
                    e_exp.append(0)
                    e_imp.append(max(dem_series[h]-e_pv[h]-e_b_d[h],0))
                SOC.append(max(SOC_min,min(SOC_max,SOC[h-1]+(e_b_c[h]/e_b_max)-(e_b_d[h]/e_b_max))))
            else:
                SOC.append(SOC_min)
                e_b_c.append(0)
                e_b_d.append(0)
                e_imp.append(dem_series[0])
                e_exp.append(0)
                 
        
        #OBJECTIVE FUNCTIONS

        # Cost objective
        c_grid.append(sum(e_imp[h]*c_elec[h] for h in t)*energy_inc[y])
        c_export.append(sum(e_exp[h]*c_exp for h in t)*energy_inc[y])
        c_capex.append((PV_cap*c_pv_capex*capex_pv[y])+(e_b_max*c_b_capex_cap*capex_bat[y])+(k_b*c_b_capex_pcs*capex_bat[y]))  
        c_opex.append((PV_cap*c_pv_opex)+(e_b_max*c_b_opex))
        
        # Emissions objective
        pv_emissions.append((emi_pv_up+emi_pv_on+emi_pv_down)*sum(e_pv[h] for h in t))
        battery_emissions.append((emi_b_up+emi_b_on+emi_b_down)*e_b_max*emi_bat[y])
        grid_emissions.append(carbon_intensity[y]*sum(e_imp[h] for h in t))
        total_emissions.append((pv_emissions[y]+battery_emissions[y]+grid_emissions[y])/1000)
        baseline_emissions.append((carbon_intensity[y]*business_annual_demand)/1000)
        avoid_emissions.append(baseline_emissions[y]-total_emissions[y]) #Calculation of avoided emisions by substracting baseline emissions (carbon intesnity*demand) to total_emissions
        tree_eq=avoid_emissions[y]/25 #tree equivalence (25kg CO2/tree) according to https://www.encon.be/en/calculation-co2-offsetting-trees
        
        baseline_carbon.append(-baseline_emissions[y]*carbon_price)
        baseline.append(baseline_carbon[y]+baseline_grid[y])
        c_carbon.append(total_emissions[y]*carbon_price)
        
        #Discount Flows
        net_flow.append(-c_capex[y]-c_opex[y]-c_grid[y]-c_carbon[y]+c_export[y])
        # net_flow_array=np.array(net_flow)
        discounted_flow.append(net_flow[y]*discount_factor[y])
        net_flow_dif.append(net_flow[y]-baseline[y])
        discounted_flow_dif.append(net_flow_dif[y]*discount_factor[y])
        
    npv=sum(discounted_flow)
    npv_dif=sum(discounted_flow_dif)
    total_investment=sum((c_capex[y]+c_opex[y])*discount_factor[y] for y in range(0,20))
    
    
    
    #Total cost should be less that Max Budget
    # c_tot<=c_max
    
    # #OUTPUTS-------------------------------------------------------------------------------------------------------------------
    
    # Dataframe containing main data output
    output_df=pd.DataFrame()
    output_df['Time']=pd.date_range(start='2019-01-01 00:00', periods=8760, freq='1H')
    output_df['Demand']=pd.concat([dem_series], ignore_index=True) #Electricity demand for 20 years, assuming no increase
    output_df['Irradiance']=pd.concat([year_irr], ignore_index=True)
    output_df['Temp_panel']=pd.concat([year_tpanel], ignore_index=True)
    output_df['Temp_amb']=year_tamb
    output_df['PV Output AC']=pd.concat([e_pv], ignore_index=True)
    output_df['PV Surplus']=pd.concat([pv_surplus], ignore_index=True)
    output_df['SOC']=SOC
    output_df['E Charging [t]']=[e_b_c[h] for h in t]
    output_df['E Discharging [t]']=[e_b_d[h] for h in t]
    output_df['Grid Import']=[e_imp[h] for h in t]
    output_df['Grid Export']=[e_exp[h] for h in t]
    output_df['Grid Cost']=[e_imp[h]*c_elec[h] for h in t]
    output_df['Grid Revenue']=[e_exp[h]*c_exp for h in t]
    output_df['Battery Energy']=[SOC[h]*e_b_max for h in t]
    output_df['PV Out Y1']=pd.concat([pv_out],ignore_index=True)
    LCOE=total_investment/(0.001*sum(e_pv)*9.950114779) #total investment/sum(e_pv)*sum of discount factors

    results=pd.DataFrame(index=range(0,20))
    results['Capex']=c_capex
    results['Opex']=c_opex
    results['Import']=c_grid
    results['Export']=c_export
    results['Carbon Costs']=c_carbon
    results['Net Flow']=net_flow
    results['Discount Factor']=discount_factor
    results['Discounted Flow']=discounted_flow
    results['Discounted Flow Dif']=discounted_flow_dif
    results['Baseline Grid']=baseline_grid
    results['Baseline Carbon']=baseline_carbon
    results['Net Flow Dif']=net_flow_dif
    results['Total Emissions']=total_emissions
    results['Baseline Emissions']=baseline_emissions
    results['Grid emissions']=grid_emissions
    results['PV emissions']=pv_emissions
    results['Battery emissions']=battery_emissions
    results['Avoided Emissions']=avoid_emissions

    #-----------------------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------------------------------------------------------
    # Grab Currrent Time After Running the Code
    end = time.time()
    
    #Subtract Start Time from The End Time to calculate elapsed time running the code in minutes
    total_time = (end - start)
    print(f' \n Elapsed time solving the iteration: {total_time:,.2f} seconds')
    return npv,results,LCOE,output_df,npv_dif,sum(avoid_emissions)

#-------------------------------------------------------------------------------------------------------------------
def iteration():
        # Grab Currrent Time Before Running the Code
    start = time.time()
    panels= [0.1,20,29,30] #Fill this list with the desired bounds for PV capacity to run the iterations

    e_battery= [0.01,1,1.5,1.75,2,2.5,2.75,3,3.25,3.5,4,5]

    npv_dif_results = pd.DataFrame(index=panels, columns=e_battery)
    emi_results = pd.DataFrame(index=panels, columns=e_battery)
    lcoe_results = pd.DataFrame(index=panels, columns=e_battery)
    for n_panel in panels:
        for bat in e_battery:
            iteration=run_model(n_panel,bat)
            # npv = iteration[0]
            npv_dif=iteration[4]
            emi = iteration[5]
            lcoe = iteration[2]
            # npv_results.loc[n_panel, bat] = npv # store npv results in df
            npv_dif_results.loc[n_panel, bat] = npv_dif # store npv differential in df
            emi_results.loc[n_panel, bat] = emi #store emissions results in df
            lcoe_results.loc[n_panel, bat] = lcoe #store LCOE results in df
            

    # Grab Currrent Time After Running the Code
    end = time.time()
    #Subtract Start Time from The End Time to calculate elapsed time running the code in minutes
    total_time = (end - start)
    print(f'\n Elapsed time solving the model: {total_time:,.2f} seconds')
    
    return npv_dif_results,lcoe_results,emi_results


def results(variable):
    irradiance=sum(variable[3]['Irradiance'])
    grid_supply=sum(variable[3]['Grid Import'])
    grid_ratio=grid_supply/sum(variable[3]['Demand'])
    pv_supply=sum(variable[3]['PV Output AC'])
    pv_ratio=1-grid_ratio
    export=sum(variable[3]['Grid Export'])
    export_ratio=export/sum(variable[3]['Demand'])
    emissions_ratio=sum(variable[1]['Avoided Emissions'])/sum(variable[1]['Total Emissions'])
    emi_decrease=sum(variable[1]['Avoided Emissions'])/sum(variable[1]['Baseline Emissions'])
    return {'Total Irradiance':irradiance,
            'Grid Supply':grid_supply,
            'Grid Ratio':grid_ratio,
            'PV Supply':pv_supply,
            'PV Ratio':pv_ratio,
            'Export':export,
            'Export Ratio':export_ratio,
            'Emissions Ratio':emissions_ratio,
            'Emissions Decrease':emi_decrease}


def plot_results(df,N_panels, e_b_max,start='2019-03-26 00:00',end='2019-03-27 00:00'):
    '''
    Plot results is a function that takes the Output DataFrame from run_model function to plot PV Output, Import, Export and Battery Energy between a certain date range

    Parameters
    ----------
    df : PandasDataFrame
        Output Dataframe.
    N_panels : float
        Number of Panels.
    e_b_max : float
        Battery Capacity in kWh.
    start : str, optional
        Start Date. The default is '2019-03-26 00:00'.
    end : str, optional
        End Date. The default is '2019-03-27 00:00'.

    Returns
    -------
    Matplotlib plot for the specified date range and a locally saved PNG image.

    '''
    #Graphs---------------------------------------------------------------------------------------------------------------------------------------------
    #New DF for plotting
    plot_df=df.loc[(df['Time'] >= start) & (df['Time'] < end)] #Plot for the week 16 of 2019
    plot_df['Hour']=plot_df['Time'].dt.hour
    print(plot_df)
    # Change the style of plot
    plt.style.use('seaborn-darkgrid') #Seaborn darkgrid style for plot
     
    # Create a color palette
    palette = plt.get_cmap('Set1') #Set 1 palette of colors for plot plt
    
    # Plot multiple lines
    plt.plot(plot_df['Hour'], plot_df['PV Output AC'], marker='', color=palette(0), linewidth=1, alpha=0.9, label='PV Output AC') #Plotting PV output
    plt.plot(plot_df['Hour'], plot_df['Grid Import'], marker='', color=palette(1), linewidth=1, alpha=0.9, label='Import') #Plotting Grid Import
    plt.plot(plot_df['Hour'], plot_df['Demand'], marker='', color=palette(2), linewidth=2, alpha=0.9, label='Demand')   #Plotting Demand
    plt.plot(plot_df['Hour'], plot_df['Grid Export'], marker='', color=palette(3), linewidth=1, alpha=0.9, label='Export')   #Plotting Grid Export
    plt.plot(plot_df['Hour'], plot_df['Battery Energy'], marker='', color=palette(4), linewidth=1, alpha=0.9, label='Battery')   #Plotting Battery Energy
    
    # Add legend
    plt.legend(loc=2, ncol=2)
     
    # Add titles
    plt.suptitle("PV output, surplus, and Grid Supply Medium Size Hotel Central London", fontsize=11, fontweight=0, color='orange')
    title=f'{N_panels} panels and {e_b_max} kWh battery'
    plt.title(title, loc='left', fontsize=9, fontweight=0, color='orange')
    plt.xlabel("Hour")
    plt.ylabel("Energy [kWh]")
    plt.plot(dpi=500)
    #Save image as png
    plt.savefig("PV output, surplus, and Grid Supply Medium Size Hotel Central London 750 MWh.png",dpi=500)
    
#------------------------------------------------------------------------------------------------------------------------------------

        



#-------------------------------------------------------------------------------------------------------------------------
