#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: 
    Fatima Iqbal fatimaiqbal0313@gmail.com 
    Christopher (Sutton) Page 
    Lewis Stetson Rowles


systems included:
    influent
    sedeimentation Tank
    Screw Press
    Aerated Stabilization Basin
    Effluent 
    
    
references:
    
Bryant, C. Updating a model of pulp and paper wastewater treatment in a partial-mix aerated stabilization basin system, 
Water Sci Technol (2010) 62 (6): 1248–1255.
https://doi.org/10.2166/wst.2010.934

Bryant, C. A simple method for analysis of the performance of aerated wastewater lagoons,
Water Sci Technol (1995) 31 (12): 211–218.
https://doi.org/10.2166/wst.1995.0489

Brazil, B.L.; Summerfelt, S.T. Aerobic treatment of gravity thickening tank supernatant,
Aquacultural Engineering (2006) 34 (2): 92-102
https://doi.org/10.1016/j.aquaeng.2005.06.001

Metcalf & Eddy, Wastewater Engineering: Treatment and Resource Recovery,
McGraw-Hill 2013, 5th edition

Juneidi, S.J., Sorour, M.T. and Aly, S.A., 2022. Proposed systematic approach for
assessing different wastewater treatment plants alternatives: Case study of Aqaba city (South Jordan).
Alexandria Engineering Journal, 61(12), pp.12567-12580.
https://doi.org/10.1016/j.aej.2022.06.044

https://hachcompany.custhelp.com/app/answers/answer_view/a_id/1020001/~/how-can-i-convert-between-nh4-%28ammonium%29-to-nh3-%28ammonia%29%3F-
    

The following wastewater (WW) constituents are used as inputs into the treatment system:
IN_COD     #mg/L
IN_SBOD    #mg/L
IN_NH4     #mg/L
IN_PO4     #mg/L
IN_TSS     #mg/L

Flowrate into the system is set: 
IN_Q       #m3/d



list of parameters that users can change in the Dash

IN_COD     #mg/L
IN_SBOD    #mg/L
IN_NH4     #mg/L
IN_PO4     #mg/L
IN_Q       #m3/d
Power      #W
solids_percent  #%

"""

from scipy.interpolate import interp1d
from scipy.integrate import odeint
import numpy as np
import pandas as pd
import lhs
import math
import scipy
from scipy import stats
from setup import setup_data
import copy


# %% LCC/LCA modeling inputs

# outputs of interest
output_perc_mid = 50
output_perc_low = 5
output_perc_high = 95


# general parameters - import spreadsheet tabs as dataframes
general_assumptions = pd.read_excel(
    'assumptions_AerationBasin.xlsx', sheet_name='General', index_col='Parameter')
design_assumptions = pd.read_excel(
    'assumptions_AerationBasin.xlsx', sheet_name='Design', index_col='Parameter')
cost_assumptions = pd.read_excel(
    'assumptions_AerationBasin.xlsx', sheet_name='Cost', index_col='Parameter')
LCA_assumptions = pd.read_excel(
    'assumptions_AerationBasin.xlsx', sheet_name='LCA', index_col='Parameter')


# number of Monte Carlo runs
n_samples = int(general_assumptions.loc['n_samples', 'expected'])

# create empty datasets to eventually store data for sensitivity analysis (Spearman's coefficients)
correlation_distributions = np.full((n_samples, n_samples), np.nan)
correlation_parameters = np.full((n_samples, 1), np.nan)
correlation_parameters = correlation_parameters.tolist()


# loading in all of the data
result = setup_data([general_assumptions, design_assumptions, LCA_assumptions],
                    correlation_distributions, correlation_parameters, n_samples)

# creates variables for each of the variables in the excel file

for key in result.keys():
    exec(f'{key} = result[key]')







#%% Sedimentation


# Retention in the settled solids    
tf                 = t0 + t                                            #min

C_deg              = tot       * max_decay                             #mg

exp_one            = math_expo ** (-k_s*t0)         

exp_second         = math_expo ** (-k_s*t0) 

exp_calculation    = exp_one   - exp_second

Cavg               = C_deg/ k  * (tf - t0)  *  exp_calculation              #mg

loss               = C_deg     - Cavg                                       #mg


SED_IN_NH4         = IN_NH4                                                 #mg/L

SED_IN_PO4         = IN_PO4                                                 #mg/L

SED_IN_SBOD        = IN_SBOD                                                #mg/L 

SED_IN_Q           = IN_Q                                                   #m3/d  

SED_mass_dry       = SED_IN_Q *  solids_percent * specific_gravity * 1000   #kg/d

SED_Solids_TS      = SED_mass_dry  *  (TS_retention /100)                   #kg/d

SED_Liquid_TS      = SED_mass_dry  *  (100 - TS_retention) / 100            #kg/d

SED_IN_COD         = IN_COD * 1000 / 1000000                                #kg/m3

SED_Sol_COD        = split_COD * SED_IN_COD  * SED_IN_Q                     #kg/d
        
SED_Liq_COD        = SED_IN_COD * SED_IN_Q - SED_Sol_COD                    #kg/d        



 
#Calculate total mass of settled solids based on final solids content 

SED_Solids_i          = SED_Solids_TS / (final_solids_content / 100)   #kg/d

SED_mass_i            = SED_mass_dry / (final_solids_content/ 100)     #kg/d  

SED_Drained_water     = SED_mass_i - SED_Solids_i                      #kg/d

SED_liquid_mass_total = SED_Drained_water + SED_Liquid_TS              #kg/d

 

#Adjust solid water

SED_Sol_imass_H2O    = SED_Drained_water 	                            #kg/d                       

Settled_solids_a = SED_Solids_i / 1000                                  #m3/d
SED_Liquid_imass     =  IN_Q - Settled_solids_a      #kg/d

#Volume of Primary Sludge 
 
SED_Sludge_Volume    = SED_Solids_i / (1000 * solids_percent)             #m3/d
               




#%% Screw Press

solids_filtrate       = Filtrate_flowrate * Filtrate_Solids

#Solids capture

SP_Solids_Capture    = ((SED_Solids_i - solids_filtrate) / SED_Solids_i) * 100        #Percentage


#Sludge Cake characteristics
 
SP_Cake_Solids       =  SED_Solids_i * (SP_Solids_Capture / 100)                      #kg/d

#Volume of Cake Solids

SP_Cake_Volume       = SP_Cake_Solids / (specific_gravity * Sludge_cake_solids * 1000) #m3/d

#Centrate Characteristics

SP_Flowrate          = SED_Sludge_Volume - SP_Cake_Volume                              #m3/d

SP_Liq_TSS           = (SED_Solids_i) * (1 - SP_Solids_Capture / 100)                  #kg/d


cake_sol_COD            = SP_cake_COD * 1000 / 1000000                                 #kg

SP_sol_COD              = cake_sol_COD * SP_Cake_Volume                                #kg/d

SP_cake_sol_F_mass      = SP_Cake_Volume  * 1500                                       #kg/d

SP_cake_sol_mass        = SP_cake_sol_F_mass / (SP_Solids_Capture / 100)               #kg/d

SP_cake_water           = SP_cake_sol_mass - SP_cake_sol_F_mass                        #kg/d
                                    
liq_COD                 = liq_COD_IN * 1000 /1000000                                   #kg/d 

SP_liq_COD              = liq_COD * SP_Flowrate                                        #kg/d
 
SP_liq_imass_H2O        = SED_Liquid_imass  - SP_cake_water                            #kg/d                                  
SP_liq_imass_H2O_a      =  SP_liq_imass_H2O / 1000
SP_waste_TS             = SED_liquid_mass_total - SED_Liquid_imass                     #kg/d

SP_polymer_imass        = dewatering_polymer_dose * SP_waste_TS / 1000                 #kg/d

# %% Aerated Stabilization Basin

# process models for ASB


# Flow in Aerated Stablization Basin is set to the Flow volume from the sedimentation Tank


Q_in = IN_Q * 3785.4118

Q = Q_in + SED_Liquid_imass + SP_liq_imass_H2O_a 

TSS_mg = SED_Liquid_TS / 1000000

TSSi = TSS_mg / Q                                                                   #mg/L 




# STEP 1: TEMP SENSITIVITY- Calculate Arrhenius temperature sensitivity coefficient

# STEP 2: k.VSS- Adjust oxidation rate to cell temperature


K_t = K_20 * theta ** (T - 20)

Influent_nbVSS = Q * X_ot / V  # g/m3.d

Ks_SBODe = Ks + IN_SBOD  # mg/L

rsu = - k * X * IN_SBOD / Ks_SBODe  # g/m3.d

r_Xt_vss = -Y * (rsu) - kd * X + fd * kd * X + Influent_nbVSS  # g/m3.d

VSS_per_day = r_Xt_vss * Q  # g/d

oxidation_rate = r_Xt_vss * Oxygen_factor  # g/m3.d

Adj_Oxidation_rate = oxidation_rate * K_t  # g/m3.d

# STEP 3: Digest Factor-Temperature-adjust benthal feedback (not Arrhenius)

# Aerobic oxidation yield of 0.5 mg biomass per mg BOD utilized

aerobic_biomass_yield = 0.5

Biomass_Growth = aerobic_biomass_yield * (Influent_BOD - Effluent_BOD)  # mg/L

# Stoichiometric ratio between the uptake of phosphorus (P) and the growth of microbial biomass
Effluent_BOD = np.nan_to_num(Effluent_BOD, nan=0.0)  # Replace nan with 0, or choose an appropriate value

valid_indices = ~np.isnan(Effluent_BOD).flatten()  # Get indices where Effluent_BOD is not nan
Influent_BOD = Influent_BOD[valid_indices]
Effluent_BOD = Effluent_BOD[valid_indices]


P_Uptake_ratio = 2.2

Growth_P_Uptake = (P_Uptake_ratio / 115) * Biomass_Growth  # mgP/L

PO4_P_fb = Growth_P_Uptake - IN_PO4 + PO4_Pe  # mgP/L

TSS_solu_P = (PO4_P_fb * 115) / P_Uptake_ratio  # mg/L

Settled_solids = TSSi + Growth - TSSe  # mg/L

Digestion_eff_P = TSS_solu_P / Settled_solids * 100  # percentage

# Stoichiometric ratio between the uptake of nitrogen (N) and the growth of microbial biomass

N_Uptake_ratio = 14

Growth_N_Uptake = (N_Uptake_ratio / 115) * Biomass_Growth  # mgN/L

NH4_N_fb = Growth_N_Uptake - IN_NH4 + NH4_Ne  # mgN/L

TSS_solu_N = (NH4_N_fb * 115) / N_Uptake_ratio  # mg/L

Digestion_eff_N = (TSS_solu_N / Settled_solids) * 100  # percentage


# Calculation of Power


P_x = aeration / Volume  # W/m3

Power = P_x * Volume  # W

# STEP 4: MIXING- Calculate Mixing Intenstiy

Mixing_Intensity = Power / Volume  # W/m3

# STEP 5: SETTLING- Calculate percent of suspended solids that settle

TSS_sett = Settled_solids / (TSSi + Growth) * 100  # percentage

# STEP 6: PARTIAL MIX Calculate ration of cell mixing intensity to complete mix

Complete_Mixing = CM_Power / CM_Volume  # W/m3

Partial_Mix = Mixing_Intensity / Complete_Mixing  # fraction

# STEP 7: K1.VSS- Adjust baseline oxidation rate to cell partial-mix level

Adj_OR_Partial_mix = Partial_Mix * Adj_Oxidation_rate  # g/m3.d


# STEP 8: Cells in series- Select number of complete-mix cells to represent hydraulics

n = 1


# STEP 9: Calculate denominator of first-order rate equation.

#Rate = {1 + K_t * tn} ** n
Rate = (1 + K_t * tn) ** n

# STEP 10: Estimate benthal feedback of SBOD5

Sett_TSS = TSSi + Growth - TSSe


# Estimated soluble BOD feedback per mg of TSS settled

SBOD_feedback_ratio = 0.3

SBOD_fb = SBOD_feedback_ratio * Sett_TSS  # mg/L


# Calculate the BOD5 removed

# Calculate the BOD5 removed

K = 2.5 * (1.06 ** (T - 20))

Effluent_BOD5 = 1 / (1 + (K * Hydrau_Reten)) * Influent_BOD

BOD5_removal = Influent_BOD - Effluent_BOD  # mg/L

x_numerator = Y * (Influent_BOD - Effluent_BOD)

x_denomenator = 1 + (kd * Hydrau_Reten)

x = x_numerator / x_denomenator  # mg/L

Px = (x * Q) / 1000  # kg/d

Px_O2 = 1.42 * Px

O2_requirement = Q * (BOD5_removal / (f * 1000))-(Px_O2)  # kg/d


Power_req_bod = O2_requirement / (Coeff_Power)  # kW


# STEP 11:  LBOD to SBOD-Conversion of SBOD6-120 to SBOD5

Exponent_cal = (math_exp) ** -kL * t_n

LBOD_to_SBOD = Initial_Concentration * Exponent_cal  # mg/L

# STEP 12: Input SBOD- Calculate total reactant SBOD5

SBOD_input = IN_SBOD + SBOD_fb  # mg/L


# STEP 13: ΔSBOD- Calculate cell effluent SBOD5


Effluent_SBOD_denominator = Hydrau_Reten * (Y * k - kd) - 1  # mg/L.d


Effluent_SBOD = Ks * (1 + Hydrau_Reten * kd) / Effluent_SBOD_denominator

SBOD_removal = IN_SBOD - Effluent_SBOD  # mg/L or g/m3

SBOD_rem_g_d = SBOD_removal * Q  # g/d


SBOD_Yield = VSS_per_day / SBOD_rem_g_d  # g VSS / sBOD


SBOD_oxygen_used = (IN_SBOD * Q) - (Effluent_SBOD * Q)  # g O2/ g sBOD

Oxygen_per_unit_SBOD = SBOD_oxygen_used / SBOD_rem_g_d

# STEP 14- P Supply-Calculate available phosphorous from all supply sources


# Load parameters from Excel file
excel_file_path = 'assumptions_AerationBasin.xlsx'
parameters_df = pd.read_excel(excel_file_path, sheet_name='Design')

# Extract parameter values
parameter_values = parameters_df.set_index(
    'Parameter').loc[:, 'expected'].to_dict()

# Extract specific parameters needed for the simulation
IN_SBOD = parameter_values['IN_SBOD']
SBODe = parameter_values['SBODe']
IN_PO4 = parameter_values['IN_PO4']
PO4_Pe = parameter_values['PO4_Pe']
TSSi = parameter_values['TSSi']
TSSe = parameter_values['TSSe']
Growth = parameter_values['Growth']
T = parameter_values['T']

# Define the function representing the differential equation


def dSBOD5_dt(SBOD5, t, k20, OPO4_interp, Kopo4, T):
    OPO4 = OPO4_interp(t)
    return -k20 * (OPO4 / (Kopo4 + OPO4)) * (1.05 ** (T - 20)) * SBOD5


Sett_TSS = TSSi + Growth - TSSe

SBOD_feedback_ratio = 0.3

SBOD_fb = SBOD_feedback_ratio * Sett_TSS  # mg/L

SBOD_input = IN_SBOD + SBOD_fb


# Set initial conditions and parameters
initial_SBOD5 = SBOD_input
k20 = 0.02
Kopo4 = 0.05

# Function to calculate OPO4 at each time step


def calculate_OPO4(t, IN_SBOD, SBODe, IN_PO4, PO4_Pe):
    aerobic_biomass_yield = 0.5
    Biomass_Growth = aerobic_biomass_yield * (IN_SBOD - SBODe)
    P_Uptake_ratio = 2.2
    Growth_P_Uptake = (P_Uptake_ratio / 115) * Biomass_Growth
    PO4_P_fb = Growth_P_Uptake - IN_PO4 + PO4_Pe
    return IN_PO4 + PO4_P_fb


#  Time points for integration
# Adjust the time range and number of points as needed
time_points = np.linspace(0, 10, 100)

# Evaluate calculate_OPO4 at each time point
OPO4_values = [calculate_OPO4(t, IN_SBOD, SBODe, IN_PO4, PO4_Pe)
               for t in time_points]

# Create an interpolation function for OPO4
OPO4_interp = interp1d(time_points, OPO4_values,
                       kind='linear', fill_value='extrapolate')

# Solve the differential equation using odeint
result = odeint(dSBOD5_dt, initial_SBOD5, time_points,
                args=(k20, OPO4_interp, Kopo4, T))


# To get the value of SBOD5 at time t = 5 (for example)
desired_time = 26.4

# Find the index in the time_points array that is closest to the desired_time
index_at_desired_time = np.abs(time_points - desired_time).argmin()

# Get the corresponding value of SBOD5 from the result array
sbod5_P_at_desired_time = result[index_at_desired_time]

# Calculate OPO4 at the desired time
opo4_at_desired_time = calculate_OPO4(desired_time, IN_SBOD, SBODe, IN_PO4, PO4_Pe)

# Calculate the remaining amount of OPO4
remaining_opo4 = IN_PO4 - opo4_at_desired_time   

print(f"SBOD5 at time {desired_time}: {sbod5_P_at_desired_time[0]}")
print(f"Remaining OPO4 at time {desired_time}: {remaining_opo4}")

def calculate_OPO4_remaining(SBOD5_at_t, initial_SBOD5, IN_PO4, PO4_Pe):
    aerobic_biomass_yield = 0.5
    Biomass_Growth = aerobic_biomass_yield * (initial_SBOD5 - SBOD5_at_t)
    P_Uptake_ratio = 2.2
    Growth_P_Uptake = (P_Uptake_ratio / 115) * Biomass_Growth
    PO4_P_fb = Growth_P_Uptake - IN_PO4 + PO4_Pe
    return IN_PO4 - Growth_P_Uptake

# Calculate SBOD5 and OPO4 at the desired time
sbod5_at_desired_time = result[index_at_desired_time][0]
opo4_remaining_at_desired_time = calculate_OPO4_remaining(sbod5_at_desired_time, initial_SBOD5, IN_PO4, PO4_Pe)


print(f"SBOD5 at time {desired_time}: {sbod5_at_desired_time}")
print(f"Remaining OPO4 at time {desired_time}: {opo4_remaining_at_desired_time}")


# Aerobic oxidation yield of 0.5 mg biomass per mg BOD utilized

aerobic_biomass_yield = 0.5

Biomass_Growth_isr = aerobic_biomass_yield * (Influent_BOD - Effluent_BOD5)  # mg/L

# Stoichiometric ratio between the uptake of phosphorus (P) and the growth of microbial biomass

P_Uptake_ratio = 2.2

Growth_P_Uptake_isr = (P_Uptake_ratio / 115) * Biomass_Growth_isr  # mgP/L

PO4_P_fb_isr = Growth_P_Uptake_isr - (IN_PO4 + remaining_opo4) # mgP/L



# STEP 15: PΔSBOD- Calculation of SBOD5 supported by phosphorous


PO4_E = kpo4 + IN_PO4

SBOD_O_P = (PO4_P / PO4_E)

SBOD_P = K20 * X_O * SBOD_O_P * Coeff_30_to_20  # mg/L

# STEP 16: P-limitΔSBOD- select higher effluent SBOD5, if phosphorous limited

PO4_E_limi = kpo4 + Limi_PO4_P

SBOD_O_limi_P = (Limi_PO4_P / PO4_E_limi)

SBOD_limi_P = K20 * X_O * SBOD_O_limi_P * Coeff_30_to_20


SBOD_removal_no_P = SBOD_removal - SBOD_P


SBOD_removal_limi_P = SBOD_removal_no_P + SBOD_limi_P

Effluent_SBOD_P_limi = Input_SBOD - SBOD_removal_limi_P  # mg/L

# STEP 17: ΔSBOD no P- calculate SBOD5 removal after phosphorus exhausted


SBOD_removal_no_P = SBOD_removal - SBOD_P  # mg/L

# STEP 18- AEROBIC ΔSBOD- Calculate overall aerobic cell effluent SBOD5

SBOD_P_1 = IN_SBOD - sbod5_P_at_desired_time

Aerobic_SBOD = SBOD_removal - SBOD_P_1  # mg/L


# STEP 19- Ox supply- Calculate available oxygen supply

Power_req_sbod = Power - Power_req_bod

DO_supply = Power_req_sbod * (Coeff_Power)  # kg/d


# Define the system of ODEs

def model(y, t):
    S, X, DO = y
    dSdt = -(Us * S / (Ks + S)) * D0 / (Ko + DO)
    dXdt = a * dSdt - d * X
    dDOdt = -a_prime * dSdt - d_prime * X + KLa * (DOs - DO)
    return [dSdt, dXdt, dDOdt]


# Set parameter values
Us = 0.1
Ks = 100
D0 = 29310
Ko = 0.1
a = 0.70
d = 0.002
a_prime = 0.34
d_prime = 0.0008
KLa = 2
DOs = 6

# Load parameters from Excel file
excel_file_path = 'assumptions_AerationBasin.xlsx'
parameters_df = pd.read_excel(excel_file_path, sheet_name='Design')

# Extract parameter values
parameter_values = parameters_df.set_index(
    'Parameter').loc[:, 'expected'].to_dict()

# Extract specific parameters needed for the simulation
Coeff_Power = parameter_values['Coeff_Power']
aeration = parameter_values['aeration']
Volume = parameter_values['Volume']
Influent_BOD = parameter_values['Influent_BOD']
Effluent_BOD = parameter_values['Effluent_BOD']
Y = parameter_values['Y']
kd = parameter_values['kd']
Hydrau_Reten = parameter_values['Hydrau_Reten']
f = parameter_values['f']


P_x = aeration / Volume  # W/m3

Power = P_x * Volume


BOD5_removal = Influent_BOD - Effluent_BOD # mg/L

x_numerator = Y * (Influent_BOD - Effluent_BOD)

x_denomenator = 1 + (kd * Hydrau_Reten)

x = x_numerator / x_denomenator  # mg/L

Px = (x * Q) / 1000  # kg/d

Px_O2 = 1.42 * Px

O2_requirement = Q * (BOD5_removal / (f * 1000))-(Px_O2)  # kg/d


Power_req_bod = O2_requirement / (Coeff_Power)

Power_req_sbod = Power - Power_req_bod

DO_supply = Power_req_sbod * (Coeff_Power)

DO_supply = DO_supply[0]

# Initial conditions
inital_SBOD = sbod5_P_at_desired_time

initial_conditions = [sbod5_P_at_desired_time, 60, DO_supply]

# Time points for integration
t = np.linspace(0, 20, 100)  # Replace with your time range

# Solve the differential equations using odeint
solution = odeint(model, initial_conditions, t)

# Print the results
for i in range(len(t)):
    print(f"At time {t[i]}:")
    print(f"BOD concentration (S) = {solution[i, 0]}")
    print(f"X concentration = {solution[i, 1]}")
    print(f"DO concentration = {solution[i, 2]}")
    print()


# Solve the differential equations using odeint
solution = odeint(model, initial_conditions, t)

# Extract BOD concentration at a specific time
specific_time = 26.4  # Replace with the desired time
# Find the index closest to the specified time
index = np.abs(t - specific_time).argmin()

# BOD concentration at the specified time
sbod_concentration_at_specific_time = solution[index, 0]
DO_concentration_at_specific_time = solution[index, 1]

print(f"At time {specific_time}:")
print(f"BOD concentration (S) = {sbod_concentration_at_specific_time}")
print(f"DO concentration) = {DO_concentration_at_specific_time}")

# STEP 20- DOΔSBOD- Calculate SBOD5 removal supported by oxygen supply


DO_SBOD_a = 1.3 + DO_supply  # kg/d

DO_SBOD_b = DO_supply / DO_SBOD_a

DO_SBOD = Oxygen_per_unit_SBOD * X_O * DO_SBOD_b  # mg/L

Combine_SBOD_removal = DO_SBOD + SBOD_limi_P  # mg/L

Intr_SBOD_removal = IN_SBOD - Combine_SBOD_removal  # mg/L

Intr_SBOD_removal_123 = IN_SBOD - sbod_concentration_at_specific_time

# Step 21: Aerobic ΔSBOD- Select higher effluent aerobic SBOD5, if oxygen limited

DO_supply_limi = Power_limi * Coeff_Power  # kg/d

DO_b = 1.3 + DO_supply_limi

DO_c = DO_supply_limi / DO_b


DO_SBOD_limi = 0.25 * X_O * DO_c  # mg/L


SBOD_rem_no_DO = IN_SBOD - DO_SBOD

SBOD_removal_limi_DO = SBOD_rem_no_DO + DO_SBOD_limi

Effluent_SBOD_DO_limi = IN_SBOD - SBOD_removal_limi_DO  # mg/L

# STEP 22: ΔSBOD no P0-Recalculate SBOD5 removal after phosphorus exhausted

SBOD_rem_no_P = SBOD_removal_limi_DO - SBOD_P  # mg/L

# Step 23: Aerobic ΔSBOD Recalculate overall aerobic cell effluent SBOD5

Overall_Aerobic_SBOD = IN_SBOD - SBOD_rem_no_P  # mg/L

# Step 24: Anoxic ΔSBOD Calculate SBOD5 removal after oxygen exhausted


DO_x = 1.3 + DO_supply_limi_x

DO_no_O = DO_supply_limi_x / DO_x


SBOD_no_O = 0.082 * X_O * DO_no_O

Anoxic_SBOD_eff = Aerobic_SBOD - SBOD_no_O  # mg/L


# Step 25:Total ΔSBOD- Calculate aerobic-plus-anoxic cell effluent SBOD5

Total_SBOD = Overall_Aerobic_SBOD + Anoxic_SBOD_eff  # mg/L

# Sep 26: Aerobic growth-Calculate new aerobic biomass growth

Aerobic_growth = IN_SBOD - Overall_Aerobic_SBOD

Aerobic_Biomass_growth = aerobic_biomass_yield * Aerobic_growth  # mg/L

# Step 27: : Anoxic growth-Calculate new anoxic biomass growth

Anoxic_growth = IN_SBOD - Anoxic_SBOD_eff  # mg/L

Anoxic_Biomass_growth = 0.3 * Anoxic_growth  # mg/L

# Step 28: Uptake N-Calculate nitrogen uptake by biomass growth

Overall_Biomass_growth = Aerobic_Biomass_growth + Anoxic_Biomass_growth

Uptake_N = N_Uptake_ratio / 115 * Overall_Biomass_growth  # mg/L

# Step 29: Uptake P-Calculate phosphorus uptake by biomass growth

Uptake_P = P_Uptake_ratio / 115 * Overall_Biomass_growth  # mg/L

b_t = b * t_h

a_b = a + b_t

TSS_removal_percent = t_h / a_b

TSS_removal = TSS_removal_percent / 100 * TSSi

Effluent_TSS = TSSi - TSS_removal


# Step 30: Settled TSS-Calculate new plus inlet suspended solids that settle

Overall_Settled_TSS = TSSi + Overall_Biomass_growth - Effluent_TSS  # mg/L


# Calculation of Effluent NH4

DO_mg = DO_supply * 1000000

Volume_liter = Volume * 1000


DO = DO_mg / Volume_liter  # mg/L

Temp_Corr_factor = 2.718 ** (0.098 * (T - 15))

pH_Corr_factor = 1 - 0.833 * (7.2 - pH)


Max_growth = Max_speci_growth * Temp_Corr_factor * \
    DO / Ko2 + DO * pH_Corr_factor  # d-1

k_o = Max_growth / Y_o  # d-1

Min_Resi_time = 1 / Y_o * k_o - kd_o  # d

Design_Resi_time = SF * (Min_Resi_time)  # d

Substrate_Utilization = (1 / Design_Resi_time + kd_o) * 1 / Y_o  # d-1

Kn = - 10 ** (0.051 * T - 1.158)  # mg/L

N = 1 - k_o / Substrate_Utilization

Effluent_NH4 = Kn / N  # mg/L

Ammonia_N = IN_NH4 * 0.9441


# Step 31: Benthal N-Calculate nitrogen feedback from settled biomass solids

NH4_N_feedback = Uptake_N - IN_NH4 + Effluent_NH4  # mg/L

# Step 32: Benthal P-Calculate phosphorus feedback from settled biomass solids

PO4_P_feedback = Uptake_P - IN_PO4 + PO4_Pe  # mg/L

# Step 33:  Nitrogen fixation- Calculate nitrogen fixation required to meet nitrogen demand

Nitrogen_Demand = IN_NH4 - Effluent_NH4

Nitrogen_Fixation = Nitrogen_Demand - Uptake_N  # mg/L

# Step 34:  Feedback check- Compare SBOD5 feedback with starting estimate- step 10

Overall_SBOD_fb = SBOD_feedback_ratio * Overall_Settled_TSS  # mg/L




# Calculation of effleunt NH3

import numpy as np

# Define the values
pKw = 14.0  # Example value at 25°C
pKb = 4.75  # Example value for ammonia at 25°C

# Assume pH is an array of pH values
pH = 7 # Example array of pH values

# Calculate the exponent for each pH value
exponent = pKw - pKb - pH

# Calculate 10^(exponent) for each pH value
fraction1 = 10 ** exponent

# Print the result
print("10^(pKw - pKb - pH):", fraction1)

fraction2 = 1 / (1 + fraction1)

kN1 = 2.71828 ** (1.57 * (pH - 8.5))

kN2 = 2.71828 ** (0.13 * (T - 20))

kN = kN1 * kN2

A = V / 3

AQ = A / Q

effluentN = AQ * kN * fraction2

effluentN1 = 1 / (1 + effluentN)

NH3 = IN_NH3 

effluentN2 = IN_NH3 * effluentN1


# Calculation of feedback


biomass_growth = 0.5 * (Influent_BOD - Effluent_BOD)

biomass_growth_array = np.full_like(IN_NH4, biomass_growth)

Uptake_N1 = (14 / 115 ) * biomass_growth_array 

NH4_N_feedback_1 = Uptake_N1 - (IN_NH4 + effluentN2)


# Step 35: ULTIMATE OXYGEN DEMAND

Effluent_sbod_cal = IN_SBOD - Intr_SBOD_removal_123

BOD = Effluent_BOD+ Effluent_sbod_cal

TKN = (effluentN2 / 1.215) + Organic_nitrogen  # mg/L

# TKN removal 

f_pH = 2.71828 ** (0.2 * (pH - 7))

AQ_k_fpH = AQ * K_t * f_pH

TKN_removal = TKN / (1 + AQ_k_fpH)

NBOD = 4.6 * TKN_removal  # mg/L


CBOD5 = BOD - NBOD

NH3 = Effluent_NH4 * 0.9441  # mg/L


Flow_rate_2 = Q * 0.000409  # ft3/s

Flow_rate_3 = Flow_rate_2 * 0.538171


UOD_2 = ((cBOD5_Multiplier * CBOD5) + 4.57 * effluentN2 ) * Flow_rate_3 * 8.34


# Aerated Stabilization Basin 2



# Define the function representing the differential equation

def dSBOD5_dt1(SBOD51, t1, k201, OPO41, Kopo41, T1):
    return -k201 * (OPO41 / (Kopo41 + OPO41)) * (1.05 ** (T - 10)) * SBOD51

# Load parameters from Excel file
excel_file_path = 'assumptions_AerationBasin.xlsx'
parameters_df = pd.read_excel(excel_file_path, sheet_name='Design')

# Extract parameter values
parameter_values = parameters_df.set_index(
    'Parameter').loc[:, 'expected'].to_dict()

# Extract specific parameters needed for the simulation
T1 = parameter_values['T1']

# Set initial conditions and parameters
initial_SBOD51 = sbod_concentration_at_specific_time
k201 = 0.02
OPO41 = opo4_remaining_at_desired_time
Kopo41 = 0.05


# Set the time points for integration
# Adjust the time range and number of points as needed
time_points1 = np.linspace(0, 10, 100)

# Solve the differential equation using odeint
result1 = odeint(dSBOD5_dt1, initial_SBOD51, time_points1,
                 args=(k201, OPO41, Kopo41, T1))


# Define your differential equation function for OPO4
def dOPO4_dt(OPO4, t, k201, OPO41, Kopo41, T1):
    dOPO4_dt_value = k201 * OPO4 * (OPO41 / (OPO41 + Kopo41))**(T1 / 12)
    return dOPO4_dt_value


result_OPO4 = odeint(dOPO4_dt, OPO41, time_points1, args=(k201, OPO41, Kopo41, T1))
# The result variable now contains the integrated values of SBOD5 over time

# Suppose you want to get the value of SBOD5 at time t = 5 (for example)
desired_time1 = 26.4

# Find the index in the time_points array that is closest to the desired_time
index_at_desired_time1 = np.abs(time_points1 - desired_time1).argmin()

# Get the corresponding value of SBOD5 from the result array
sbod5_P_at_desired_time1 = result1[index_at_desired_time1]

opo4_at_desired_time1 = result_OPO4[index_at_desired_time1]

# Calculate the remaining amount of OPO4 after consumption for SBOD removal
opo4_left_after_consumption = opo4_at_desired_time1 - OPO41

print(f"SBOD5 at time {desired_time1}: {sbod5_P_at_desired_time1[0]}")
print(f"Remaining OPO4 at time {desired_time1}: {opo4_left_after_consumption}")

import numpy as np
import pandas as pd
from scipy.integrate import odeint

# Define the function representing the differential equation for SBOD5

def dSBOD5_dt1(SBOD51, t1, k201, OPO41, Kopo41, T1):
    return -k201 * (OPO41 / (Kopo41 + OPO41)) * (1.05 ** (T1 - 10)) * SBOD51

# Define the function representing the differential equation for OPO4

def dOPO4_dt(OPO4, t, SBOD51, k201, OPO41, Kopo41, T1):
    # Assuming the OPO4 is consumed in proportion to the SBOD5 degradation
    return -k201 * OPO4 * (SBOD51 / (Kopo41 + SBOD51)) * (1.05 ** (T1 - 10))

# Load parameters from Excel file
excel_file_path = 'assumptions_AerationBasin.xlsx'
parameters_df = pd.read_excel(excel_file_path, sheet_name='Design')

# Extract parameter values
parameter_values = parameters_df.set_index('Parameter').loc[:, 'expected'].to_dict()

# Extract specific parameters needed for the simulation
T1 = parameter_values['T1']

# Set initial conditions and parameters
initial_SBOD51 = sbod5_P_at_desired_time  # SBOD5 at the start of this calculation
k201 = 0.02
OPO41 = opo4_remaining_at_desired_time  # Initial OPO4 at the start of this calculation
Kopo41 = 0.05

# Set the time points for integration
# Adjust the time range and number of points as needed
time_points1 = np.linspace(0, 10, 100)

# Solve the differential equation for SBOD5 using odeint
result1 = odeint(dSBOD5_dt1, initial_SBOD51, time_points1, args=(k201, OPO41, Kopo41, T1))

# Solve the differential equation for OPO4 using odeint
# Use the SBOD5 values from result1 in the OPO4 equation
result_OPO4 = odeint(dOPO4_dt, OPO41, time_points1, args=(initial_SBOD51, k201, OPO41, Kopo41, T1))

# Suppose you want to get the value of SBOD5 and OPO4 at a specific time, say t = 26.4
desired_time1 = 26.4

# Find the index in the time_points array that is closest to the desired_time
index_at_desired_time1 = np.abs(time_points1 - desired_time1).argmin()

# Get the corresponding values of SBOD5 and OPO4 from the result arrays
sbod5_P_at_desired_time1 = result1[index_at_desired_time1]
opo4_at_desired_time1 = result_OPO4[index_at_desired_time1]

# Calculate the remaining amount of OPO4 after consumption for SBOD removal
# Here we are using the OPO4 value directly from the result_OPO4 array
# There's no need to subtract OPO41 as it was already adjusted in the differential equation
opo4_left_after_consumption1 = opo4_at_desired_time1

print(f"SBOD5 at time {desired_time1}: {sbod5_P_at_desired_time1[0]}")
print(f"Remaining OPO4 at time {desired_time1}: {opo4_left_after_consumption[0]}")




Growth_P_Uptake = (P_Uptake_ratio / 115) * biomass_growth  # mgP/L

PO4_P_fb1 = Growth_P_Uptake -  (opo4_remaining_at_desired_time +  opo4_left_after_consumption) # mgP/L



Calculated_Effluent_PO4 = (opo4_left_after_consumption) + (PO4_P_fb1)


Px_1B = aeration_1B / Volume_1B

Power_1B = Px_1B * Volume_1B

K1 = 0.062 * (1.042 ** T)

Effluent_BOD5_1B = 1 / (1 + (K1 *Hydrau_Reten)) * Effluent_BOD

BOD5_removal_1B = BOD5_removal - Effluent_BOD_1B    #mg/L


x_numerator_1B = Y * (BOD5_removal - Effluent_BOD_1B)

x_denomenator_1B = 1 + (kd * Hydrau_Reten)

x_1B = x_numerator_1B / x_denomenator_1B  # mg/L

Px_1B = (x_1B * Q) / 1000  # kg/d

Px_O2_1B = 1.42 * Px_1B

O2_requirement_1B = Q * (BOD5_removal_1B / (f * 1000))-(Px_O2_1B)  # kg/d


Power_req_bod_1B = O2_requirement_1B / (Coeff_Power)  # kW


# STEP 19- Ox supply- Calculate available oxygen supply

Power_req_sbod_1B = Power_1B - Power_req_bod_1B

DO_supply_1B = Power_req_sbod_1B * (Coeff_Power)  # kg/d


# STEP 13: ΔSBOD- Calculate cell effluent SBOD5

SBOD_removal_1B = Effluent_SBOD - SBOD_1B  # mg/L or g/m3

SBOD_rem_g_d_1B = SBOD_removal_1B * Q  # g/d


SBOD_Yield_1B = VSS_per_day / SBOD_rem_g_d_1B  # g VSS / sBOD


SBOD_oxygen_used_1B = (Effluent_SBOD * Q) - \
    (SBOD_removal_1B * Q)  # g O2/ g sBOD

Oxygen_per_unit_SBOD_1B = SBOD_oxygen_used_1B / SBOD_rem_g_d_1B

SBOD_P_1B = Effluent_sbod_cal - sbod5_P_at_desired_time1

Aerobic_SBOD_1b = SBOD_removal_1B - SBOD_P_1B

# STEP 20- DOΔSBOD- Calculate SBOD5 removal supported by oxygen supply


# Define the system of ODEs

def model1(y1, t1):
    S1, X1, DO1 = y1
    dSdt1 = -(Us1 * S1 / (Ks1 + S1)) * D01 / (Ko1 + DO1)
    dXdt1 = a1 * dSdt1 - d1 * X1
    dDOdt1 = -a_prime1 * dSdt1 - d_prime1 * X1 + KLa1 * (DOs1 - DO1)
    return [dSdt1, dXdt1, dDOdt1]


# Set parameter values
Us1 = 0.1  # maximum rate of sbod removal
Ks1 = 100  # saturation coefficient for sBOD
D01 = 12634  # concentration of dissolved oxygen
Ko1 = 0.1  # saturation coefficient for DO
a1 = 0.70   # yield coefficient
d1 = 0.002   # yield cofficient
a_prime1 = 0.34  # oxygen-use coefficient
d_prime1 = 0.0008  # oxygen-use coefficient
KLa1 = 2  # overall-oxygen transfer rate
DOs1 = 6  # saturation concentration of dissolved oxygen

# Load parameters from Excel file
excel_file_path = 'assumptions_AerationBasin.xlsx'
parameters_df = pd.read_excel(excel_file_path, sheet_name='Design')

# Extract parameter values
parameter_values = parameters_df.set_index(
    'Parameter').loc[:, 'expected'].to_dict()

# Extract specific parameters needed for the simulation
Coeff_Power = parameter_values['Coeff_Power']
aeration_1B = parameter_values['aeration_1B']
Volume_1B = parameter_values['Volume_1B']
Effluent_BOD_1B = parameter_values['Effluent_BOD_1B']
Y = parameter_values['Y']
kd = parameter_values['kd']
Hydrau_Reten = parameter_values['Hydrau_Reten']
f = parameter_values['f']


K1 = 2.5 * (1.06 ** (T1 - 20))   

Effluent_BOD5_1B = 1 / (1 + (K1 *Hydrau_Reten)) * Effluent_BOD5

BOD5_removal_1B = BOD5_removal - Effluent_BOD_1B    #mg/L


x_numerator_1B = Y * (BOD5_removal - Effluent_BOD_1B)

x_denomenator_1B = 1 + (kd * Hydrau_Reten)

x_1B = x_numerator_1B / x_denomenator_1B  # mg/L

Px_1B = (x_1B * Q) / 1000  # kg/d

Px_O2_1B = 1.42 * Px_1B

O2_requirement_1B = Q * (BOD5_removal_1B / (f * 1000))-(Px_O2_1B)  # kg/d


Power_req_bod_1B = O2_requirement_1B / (Coeff_Power)  # kW


Px_1B = aeration_1B / Volume_1B

Power_1B = Px_1B * Volume_1B

Power_req_sbod_1B = Power_1B - Power_req_bod_1B

DO_supply_1B = Power_req_sbod_1B * (Coeff_Power)

# Initial conditions
initial_SBOD1 = sbod5_P_at_desired_time1
# Replace with your initial values
initial_conditions1 = [initial_SBOD1, 20, 12643]

# Time points for integration
t1 = np.linspace(0, 20, 100)  # Replace with your time range

# Solve the differential equations using odeint
solution1 = odeint(model1, initial_conditions1, t1)

# Print the results
for i in range(len(t1)):
    print(f"At time {t1[i]}:")
    print(f"BOD concentration (S) = {solution1[i, 0]}")
    print(f"X concentration = {solution1[i, 1]}")
    print(f"DO concentration = {solution1[i, 2]}")
    print()

# Solve the differential equations using odeint
solution1 = odeint(model1, initial_conditions1, t1)

# Extract BOD concentration at a specific time
specific_time1 = 26.4  # Replace with the desired time
# Find the index closest to the specified time
index1 = np.abs(t1 - specific_time1).argmin()

# BOD concentration at the specified time
sbod_concentration_at_specific_time1 = solution1[index1, 0]

aerobic_biomass_yield = 0.5

Biomass_Growth1 = aerobic_biomass_yield * (Effluent_BOD - Effluent_BOD_1B)  # mg/L



# Calculation of feedback1


biomass_growth1 = 0.5 * ( Effluent_BOD5 - Effluent_BOD5_1B )

biomass_growth_array1 = np.full_like(IN_NH4, biomass_growth1)

Uptake_N12 = (14 / 115 ) * biomass_growth1

NH4_N_feedback_12 = Uptake_N1 - (effluentN2 + NH4_N_feedback_1) 

# Define the values
pKw = 14.0  # Example value at 25°C
pKb = 4.75  # Example value for ammonia at 25°C

# Assume pH is an array of pH values
pH1 = 7 # Example array of pH values

# Calculate the exponent for each pH value
exponent1 = pKw - pKb - pH1

# Calculate 10^(exponent) for each pH value
fraction1 = 10 ** exponent1

# Print the result
print("10^(pKw - pKb - pH):", fraction1)

fraction2_1 = 1 / (1 + fraction1)

kN1_1 = 2.71828 ** (1.57 * (pH1 - 8.5))

kN2_1 = 2.71828 ** (0.13 * (T1 - 20))

kN_1 = kN1_1 * kN2_1

A = V / 4

AQ = A / Q

effluentN1_1 = AQ * kN_1 * fraction2_1

effluentN_12 = 1 / (1 + effluentN1_1)

NH3_i = effluentN2 + NH4_N_feedback_12

effluentN2_1 =  NH3_i  * effluentN_12


# Stoichiometric ratio between the uptake of phosphorus (P) and the growth of microbial biomass


P_Uptake_ratio = 2.2

Growth_P_Uptake1 = (P_Uptake_ratio / 115) * Biomass_Growth1  # mgP/L

PO4_P_fb1 = Growth_P_Uptake - opo4_remaining_at_desired_time + opo4_left_after_consumption1  # mgP/L

Effluent_PO4 = opo4_left_after_consumption1 + PO4_P_fb1

Effluent_sbod_cal1 = Effluent_sbod_cal - sbod_concentration_at_specific_time1

BOD1 =  Effluent_BOD_1B +Effluent_sbod_cal1

Ammonia_N1 = IN_NH41 * 0.9441

TKN1 = (effluentN2_1 / 1.215) + Organic_nitrogen1  # mg/L

# TKN removal 

f_pH1 = 2.71828 ** (0.2 * (pH1 - 7))

AQ_k_fpH1 = AQ * K_t * f_pH1

TKN_removal1 = TKN_removal / (1 + AQ_k_fpH1)

NBOD1 = 4.6 * TKN_removal1

CBOD51 = BOD1 - NBOD1

UOD_3 = ((cBOD5_Multiplier * CBOD51) + 4.57 * effluentN2_1 ) * Flow_rate_3 * 8.34






import pandas as pd
import matplotlib.pyplot as plt

# Load the data
data = {
    'Time': ["2020-05-06 00:00:00", "2020-05-11 00:00:00", "2020-05-18 00:00:00", "2020-05-25 00:00:00", "2020-05-27 00:00:00", "2020-06-08 00:00:00", "2020-06-10 00:00:00", "2020-06-17 00:00:00", "2020-06-22 00:00:00", "2020-07-27 00:00:00", "2020-08-03 00:00:00", "2020-08-05 00:00:00", "2020-08-12 00:00:00", "2020-08-17 00:00:00", "2020-08-26 00:00:00", "2020-09-23 00:00:00", "2020-10-12 00:00:00", "2020-10-14 00:00:00", "2020-10-19 00:00:00", "2020-10-21 00:00:00", "2020-10-26 00:00:00", "2020-10-28 00:00:00", "2020-11-04 00:00:00", "2020-11-09 00:00:00", "2020-11-22 00:00:00", "2020-11-30 00:00:00", "2020-12-09 00:00:00", "2020-12-27 00:00:00", "2021-01-27 00:00:00", "2021-02-03 00:00:00", "2021-02-15 00:00:00", "2021-02-17 00:00:00", "2021-03-03 00:00:00", "2021-03-08 00:00:00", "2021-03-29 00:00:00", "2021-04-14 00:00:00", "2021-04-19 00:00:00", "2021-05-03 00:00:00", "2021-05-05 00:00:00", "2021-05-17 00:00:00", "2021-05-31 00:00:00", "2021-06-02 00:00:00", "2021-06-09 00:00:00", "2021-06-30 00:00:00", "2021-07-05 00:00:00", "2021-07-07 00:00:00", "2021-07-12 00:00:00", "2021-07-14 00:00:00", "2021-08-11 00:00:00", "2022-08-03 00:00:00", "2022-08-10 00:00:00", "2022-08-24 00:00:00", "2022-08-29 00:00:00", "2022-09-05 00:00:00", "2022-09-12 00:00:00", "2022-09-14 00:00:00", "2022-09-21 00:00:00", "2022-09-26 00:00:00", "2022-09-28 00:00:00", "2020-05-13 00:00:00", "2020-06-15 00:00:00", "2020-08-31 00:00:00", "2020-09-09 00:00:00", "2020-09-21 00:00:00", "2020-10-07 00:00:00", "2020-11-04 00:00:00", "2020-11-11 00:00:00", "2020-11-18 00:00:00", "2020-12-02 00:00:00", "2020-12-07 00:00:00", "2021-01-06 00:00:00", "2021-01-18 00:00:00", "2021-02-08 00:00:00", "2021-02-22 00:00:00", "2021-03-01 00:00:00", "2021-03-10 00:00:00", "2021-03-22 00:00:00", "2021-03-24 00:00:00", "2021-03-31 00:00:00", "2021-04-28 00:00:00"],
    'Actual': [21, 8, 16, 23, 36, 26, 16.03, 13.36, 14.48, 13, 14, 28, 12.2, 10, 14, 18, 15.5, 16.5, 14.5, 15, 16, 12.5, 14, 19, 20, 29, 19, 39, 30, 20, 21, 17, 23, 37, 27, 30, 31, 36, 34, 31.5, 35, 34.5, 18, 20, 20, 30, 14.5, 14, 17, 17, 10, 16, 24, 25, 16, 10, 11, 11, 23, 13, 11, 12, 16, 12, 14, 16, 21, 18, 13, 19, 29, 21, 23, 14, 20, 14, 20, 14, 27, 17],
    'Process model': [23.9, 6.32, 16.744, 20.78, 31, 25.77, 18.8, 13.64, 15, 10.2, 11.38, 28.1, 12.2, 19.2, 13.89, 13.04, 12.7, 18, 13.35, 13.84, 13.53, 15, 7, 18.8, 18.1, 24.1, 21.8, 35.6, 26.4, 16.1, 26, 19.7, 24.7, 31, 18.8, 23.8, 32.5, 38.5, 34.25, 31.7, 46, 37.69, 11.38, 22.7, 15.9, 30.5, 15.51, 15.5, 13.5, 16, 7, 16.5, 27.4, 25.6, 20, 8.3, 13, 12.3, 25.6, 7, 6.2, 19, 10.8, 4.57, 7, 25.62, 31.9, 28.51, 47, 34.9, 44.5, 42, 34.6, 3.35, 11.6, 24.7, 7.51, 24, 54.3, 28],
    'Temperature correct': [23.9, 6.32, 16.744, 20.78, 31, 25.77, 18.8, 13.64, 15, 10.2, 11.38, 28.1, 12.2, 19.2, 13.89, 13.04, 12.7, 18, 13.35, 13.84, 13.53, 15, 7, 18.8, 18.1, 24.1, 21.8, 35.6, 26.4, 16.1, 26, 19.7, 24.7, 31, 18.8, 23.8, 32.5, 38.5, 34.25, 31.7, 46, 37.69, 11.38, 22.7, 15.9, 30.5, 15.51, 15.5, 13.5, 16, 7, 16.5, 27.4, 25.6, 20, 8.3, 13, 12.3, 25.6, 10, 9.12, 11.5, 14.6, 10.7, 14.4, 15.9, 24, 20, 11.6, 18.06, 32, 21, 18, 9, 22.4, 16.5, 16, 13.97, 19.47, 18.8]
}
# Check lengths of each column
for key, value in data.items():
    print(f"Length of {key}: {len(value)}")


# Convert to DataFrame
df = pd.DataFrame(data)

# Convert 'Time' column to datetime
df['Time'] = pd.to_datetime(df['Time'])

# Set 'Time' column as index
df.set_index('Time', inplace=True)

# Plot the data
df.plot(figsize=(14, 7))

plt.xlabel('Time')
plt.ylabel('Values')
plt.grid(True)
plt.show()


#power
power= O2_requirement / Hydrau_Reten * 24 / aeration_efficiency #kWh


#ASB

#Electricity_Consumption
Electricity_Consumption= power * Hydrau_Reten * 24 #kWh

#electricity_annual_costs
electricity_annual_costs= Electricity_cost * 365 * Electricity_Consumption #$/year


