"""
Radiator heat disipation trough conduction and convection

Approach: Assumed that each row will get equal amound of flow of water at same temperature, compute for 1 row to then multiply by the number of rows.
          Divide a row into cells per fins, to compute Q for the first cell to then compute the next cell with inlet water temp being previous outlet water temp.

Assumtion: constant heat of surface of the tube across each cells 
           constant heat of surface across the depth of radiator
           correction factor for delta temp log mean of cross flow = 1
           take average conduction resitance for fins 
"""
import math 

#Inputs
flr_w = 0.005 #m^3/s          Flow rate of water
Vm_a = 100 #m/s               Mean Speed of air
T_in_w = 80 #c                Temperature initial of water 
T_in_a = 15 #c                Temperature initial of air
width = 0.5 #m                Width of the radiator
hight = 0.4 #m                Hight of radiator
depth = 0.05 #m               Depth of the radiator
no_fins = 100 #               Number of fins in one stage
hight_fins = 0.01 #           Hight of a fins
hight_tubes = 0.005 #         Hight of the water tubes 
thickness_tubes = 0.001 #     Thickness of wall of the tubes
thickness_fins = 0.001 #      Thickness of the fins
Cp_w = 4070 #J/kg*K           Specific heat capacity of water  *at 80c, 1bar
Cp_a = 1006 #J/Kg*K           Specific heat capacity of air  *at 15c, 1bar
p_w = 971 #Kg/m^3             Density of water  *at 80c, 1bar
p_a = 1.225 #Kg/m^3           Density of air  *at 15c at 1bar
K_w = 0.67002 #W/m*K          Thermoconductivity of water  *at 80c, 1bar
K_a = 0.025219 #W/m*K         Thermoconductivity of air  *at 15c, 1bar
K_m = 237 #W/m*K              Thermoconductivity of aluminum
mu_w = 0.00034599 #N*s/m^2    Dynamic viscosity of water  *at 80c, 1bar
mu_a = 0.000014657 #N*s/m^2   Dynamic viscosity of air *at 15c, 1bar
laminar_flow_w = False #      Initial flow condition for water True or false
laminar_flow_a = True #       Initial flow condition for air True or false
Ts_percent_error = 0.1 #%     Temperature surface percent error for empirical loop
Q_sum = 0 #                   Initialization of Q per row
T_out_a_sum = 0 #             Initialization of Temp out air average

#Script
def Convection_heat_transfer_coefficient(Ac, P, p, Vm, Cp, mu, K, x, ratio, laminar_flow):
    D = (4*Ac)/P #hydraulic diameter
    Re = (p*Vm*D)/mu #Reynolds Number
    Pr = (Cp*mu)/K #Prandlt Number
    Gz = (D*Re*Pr)/x #Graetz Number
    if laminar_flow == True: #If laminar
        Nu = 0.382365*ratio+2.67011 #Nusselt Number
    else:                    #If turbulent
        Nu = 0.023*(Re**(4/5))*(Pr**0.3) #Nusselt Number
    Nu_ = Nu+(0.0668*Gz)/(1+0.04*Gz**(2/3)) #Average Nusselt Number
    h_ = (Nu_*K)/D #Average heat transfer coefficient
    return h_ 
 
no_stages = hight/(hight_fins+hight_tubes) #Number of stages
Ac_w = (hight_tubes-2*thickness_tubes)*(depth-2*thickness_tubes) #Cross section area of the inside of the water tube
Vm_w = flr_w/(no_stages*Ac_w)

for i in range(no_fins): #Comput Q for each cell with T_in_w being the previous T_out_w
    T_out_w = T_in_w-((1/5)*(T_in_w-T_in_a)) #Initialization of ampirical values for otlet temperature of water 
    T_out_a = T_in_a+((1/3)*(T_in_w-T_in_a)) #Initialization of ampirical values for otlet temperature of air
    T_out_w_og = T_in_w*2  #Initialization of otlet temperature of water
    T_out_a_og = T_in_a/2 #Initialization of otlet temperature of air
    while ((abs(T_out_w-T_out_w_og))/T_out_w)*100 > Ts_percent_error or ((abs(T_out_a-T_out_a_og))/T_out_a)*100 > Ts_percent_error: #Empirical computation for T_out until it reaches the wanted percent error
         
        x_w = width/no_fins #Lenght the water has to travel in one cell 
        Ac_w = (hight_tubes-2*thickness_tubes)*(depth-2*thickness_tubes) #Cross section area of the inside of the water tube
        Ac2_w = Ac_w/2 #Half the cross section area of the inside of the water tube
        ratio_w = depth/(hight_tubes-2*thickness_tubes) #Ratio of width to hight of the inside of the water tube
        P_w = (2*(hight_tubes-2*thickness_tubes))+(2*(depth-2*thickness_tubes)) #Cross section perimiter of the inside of the water tube
        P2_w = P_w/2 #Half the cross section perimiter of the inside of the water tube
        A_w = (x_w*(depth-2*thickness_tubes)) #Surface area of the bottom face inside the water tube
        h_w = Convection_heat_transfer_coefficient(Ac_w, P_w, p_w, Vm_w, Cp_w, mu_w, K_w, x_w, ratio_w, laminar_flow_w) #Average heat transfer coefficient for water
        R_conv_w = 1/(h_w*A_w) #Resistance for convection of water to the tube
    
        L_base = thickness_tubes #Thickness of the tube's wall
        R_cond = L_base/(K_m*x_w*depth) #Resitance for conduction throught the tube's wall 
    
        P_a = (2*hight_fins)+(2*((width/no_fins)-thickness_fins)) #Cross section perimiter of one cell copartment for air
        P_fin = hight_fins #Cross section perimiter of a fin
        Ac_a = hight_fins*(width/no_fins) #Cross section area for one cell
        Ac2_a = Ac_a/2 #Half of the cross section area for one cell
        x_a = depth #Lenght the air has to travel
        ratio_a = ((width/no_fins)-thickness_fins)/hight_fins #Ratio of width to hight of one cell compartment for air
        h_a = Convection_heat_transfer_coefficient(Ac_a, P_a, p_a, Vm_a, Cp_a, mu_a, K_a, x_a, ratio_a, laminar_flow_a) #Average heat transfer coefficient for water
        A_fins = hight_fins*depth #Surface area of one half a fine for both sides
        R_cond_fins = ((hight_fins/4)/(K_m*thickness_fins*depth)) #Average resistance for conduction of half a fin
        R_conv_fins = 1/(h_a*A_fins) #Resitance for convection of half a fin
        R_fins = R_cond_fins+R_conv_fins #Resitance total of half a fin
        A_a = (x_w*depth)-(thickness_fins*depth) #Surface area of the base in contact with air (not including fin)
        R_conv_base = 1/(h_a*A_a) #Resitance for convection of the base 
        R_conv_a = ((1/R_conv_base)+(1/(R_fins)))**(-1) #Resistance total from metal to air

        Delta_T1 = T_in_w-T_out_a #Change in temperature between water inlet temp to air outlet temp
        Delta_T2 = T_out_w-T_in_a #Change in temperature between water outlet temp to air inlet temp
        Delta_T_lm = (Delta_T1-Delta_T2)/(math.log(Delta_T1/Delta_T2)) #Delta Temperature log mean for one cell between water and air at a cross flow 

        R_tot = R_conv_w+R_cond+R_conv_a #Total resitance of one cell from water to air

        Q = (Delta_T_lm*0.5)/R_tot #Heat transfer rate for one cell
        
        T_out_w_og = T_out_w #Initialization outlet water temperature original
        T_out_a_og = T_out_a #Initialization outlet air temperature original
        time_w = x_w/Vm_w #Time it take water to travel across one cell
        volume_w = ((hight_tubes/2)-thickness_tubes)*x_w*(depth-(2*thickness_tubes)) #Volume of water in one cell at any given momment
        mass_w = volume_w*p_w #mass of water in one cell
        T_out_w = (-Q*time_w/(mass_w*Cp_w))+T_in_w #New found empiricall value for the outlet temperature of water
        time_a = x_a/Vm_a #Time it take air to travel across one cell or depth of radiator
        volume_a = (hight_fins/2)*(x_w-thickness_fins)*(depth) #Volume of air in one cell at any given momment
        mass_a = volume_a*p_a #mass of air in one cell
        T_out_a = (Q*time_a/(mass_a*Cp_a))+T_in_a #New found empiricall value for the outlet temperature of air
        
    T_in_w = T_out_w #Intlet temperature of water now being the previous outlet water temperature
    T_out_a_sum += T_out_a #Sum of oultlet temperature of water
    Q_sum += Q #Total heat transfer rate for one row
    
Q_tot = Q_sum*no_stages #Total heat transfer rate for the radiator
T_out_a_avr = T_out_a_sum/no_fins #Average outlet water remperature

print(f"The total thermal power dissipated trough this radiator is: {round(Q_tot)} W")
print(f"The output temperature of water is: {round(T_out_w)} c")
print(f"The output temperature of air is: {round(T_out_a_avr)} c")


