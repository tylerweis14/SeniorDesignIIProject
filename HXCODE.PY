# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 20:57:11 2020

@author: Jorge
"""


### CODE FOR SHELL-AND-TUBE HX WITH TWO TUBE PASSES/ONE SHELL PASS

import numpy as np


CHOICE= 2 ### 1== Shell & Tube HX and 2== Printed Circuit Heat Exchanger
 
T_h_in= 923 #ALL TEMPERATURES LISTED IN KELVIN, K 
T_c_in= 573
T_c_out= 873


#FLUID PROPERTIES -- FLiBe (shellside hot)

rho_h= 2518-0.406*(T_h_in+T_c_in)/2 #density in kg/m^3
mu_h= 0.000116*np.exp(3755/((T_h_in+T_c_in)/2)) #viscosity in Pa*s
c_ph= 2415.78  #specific heat in J/kg*K
k_h= 0.629697+0.0005*((T_h_in+T_c_in)/2) #thermal conductivity in W/mK
Pr_h= 13.525

#FLUID PROPERTIES -- Solar Salt (tubeside cold)

rho_c= 1804
mu_c= 0.00169
c_pc= 1520 
k_c= 0.530 
Pr_c= 4.85



#TUBE PROPERTIES 

d_o= 0.02 #outer tube diameter in m 
t_w= 0.001 #tube wall thickness 
d_i= d_o-2*t_w #inner tube diameter 

#GUESSES 

U= 100
U_guess= 200 #Overall HT Coefficient in W/m^2*K
v_tube_guess= 1.5  #Tube velocity in m/s 
#Energy Balance 

mdot_h= 1500 #mass flow rate in kg/s 
mdot_c= 1550 
Qdot= mdot_c*c_pc*(T_c_out-T_c_in)
T_h_out= T_h_in-mdot_c*c_pc*(T_c_out-T_c_in)/(mdot_h*c_ph)

#Log Mean Temperature Difference 

DeltaT_Mean= ((T_h_in-T_c_out)-(T_h_out-T_c_in))/(np.log((T_h_in-T_c_out)/(T_h_out-T_c_in)))
if mdot_h*c_ph>mdot_c*c_pc:
    P= mdot_c*c_pc/(mdot_h*c_ph)
else: P= mdot_h*c_ph/(mdot_c*c_pc)          #P is Capacity Ratio 
R= Qdot/(np.min([mdot_h*c_ph,mdot_c*c_pc])*(T_h_in-T_c_in)) #R is Effectiveness 
Fnum= np.sqrt((R**2)+1)*np.log((1-P)/(1-P*R))
Fden= (R-1)*np.log(((2/P)-1-R+np.sqrt((R**2)+1))/(2/P)-1-R-np.sqrt((R**2)+1))
F= Fnum/Fden
if np.isnan(F) is True :              #Determining if F is a real number 
    DeltaT_Mean= DeltaT_Mean*F
else: DeltaT_Mean= DeltaT_Mean*1 


if CHOICE==1: 

    #TEMA Correlations 
    
    L_tp= 1.25*d_o   #Tube Pitch in meters
    L_tb= 0.0008     #Tube-to-baffle diametral clearance in meters 
    L_b= 0.005       #Tube sheet thickness in meters 
    SS= 0.2          #No. sealing strips per crossflow row
    B_c= 0.25        #Baffle-cut as decimal of shell diameter
    
    #ITERATIONS TO GET U_GUESS=U
    
    Ui= np.arange(200, 10100, 100)
    i= 0
    
    while U/U_guess<0.99 or U/U_guess>1.01 :
         
        U_guess= Ui[i] 
            
    
    #Number of Tubes 
            
        N_p= 2  #Number of Shell Passes 
        N_t= mdot_c*N_p*4/(np.pi*(d_i**2)*v_tube_guess*rho_c)
        
        #Length of Tube 
        
        A_req= Qdot/(U_guess*DeltaT_Mean)  #Heat Transfer Area
        L_tube= A_req/(np.pi*d_o*N_t)
        
        #Rest of TEMA Correlations 
        
        D_ctl= np.sqrt(N_t*0.86*L_tp**2/0.78) #Centerline Tube Limit for 30 degree layout
        D_s= 0.995*(D_ctl+12/1000+d_o)        #Shell Diameter 
        D_otl= 0.995*D_s-12/1000              #Outer Tube Limit 
        L_bc= 0.3*D_s                         #Baffle Spacing
        L_sb= (3.1+0.004*D_s)/1000            #Shell-to-Baffle Diametral Clearance 
        L_bb= (12+5*D_s)/1000                 #Bundle-to-Shell Diametral Clearance 
        
        #Number of Baffles 
        
        t_b= 0.01  #Baffle thickness in meters 
        N= (L_tube/(L_bc+t_b))-1 #Number of Baffles 
        
        #Wall Interferance HT Coefficient 
        
        k_w= 14.9   #Thermal Conductivity of Material AISI 304 Stainless Steel
        U_wall= k_w/t_w    #HTC of Wall in W/Km^2 
        
        #DIMENSIONLESS PARAMETERS 
        
        S_m= L_bc*(D_s-D_otl+((D_otl-d_o)/L_tp)*(L_tp-d_o)) #minimal crossflow area at bundle centerline 
        v_max= mdot_h/(S_m*rho_h) #MAX velocity on shellside 
        
        Re= rho_h*v_max*d_o/mu_h  #Reynolds Number
        Pr= c_ph*mu_h/k_h         #Prandtl Number 
        
        if Re<=300:     #FACTORS FOR Nu for in line banks
            a=0.742 
        if Re<=300: 
            m= 0.431     
        if 300<Re<=200000:
            a= 0.211
        if 300<Re<=200000:
            m= 0.651 
        if 200000<Re<=2000000:
            a= 0.116 
        if 200000<Re<=2000000:
            m= 0.7    
        
        Nu= a*(Re**m)*Pr**0.34 
        
        #Shellside HT Coefficient - IDEAL
        
        U_shellside_ideal= Nu*k_h/d_o 
        
        #Correction Factors for Shellside HTC
        
        L_c= B_c*D_s  #Baffle Length 
        F_c= (1/np.pi)*(np.pi+(2*((D_s-2*L_c)/D_otl)*np.sin(2*np.arccos((D_s-2*L_c)/D_otl)))-np.arccos((D_s-2*L_c)/D_otl))
        
        J_c= 0.55+0.72*F_c  #CONFIGURATION CORRECTION FACTOR 
        
        S_sb= D_s*L_sb*(np.pi-np.arccos(1-(2*L_c/D_s)))/2 #Shell to Baffle Leakage Area
        S_tb= np.pi*d_o*L_tb*N_t*(1+F_c)/4                #Tube to Baffle Leakage Area
        r_lm= (S_sb+S_tb)/S_m 
        r_s= S_sb/(S_sb+S_tb) 
        x= -0.15*(1+r_s)+0.8 
        
        J_L= (0.44/(1-r_s))+(1-0.44*(1-r_s))*np.exp(-2.2*r_lm)  #LEAKAGE Correction Factor 
        
        F_bp= (D_s-D_otl)*L_bc/S_m  #Bypass Area Fraction
        N_c= D_s*(1-(2*L_c/D_s))/(0.866*L_tp) #Crossflow Rows for Triangular Layout 
        N_ss= N_c*SS 
        
        if SS>0.5:       #Bypass Correction Factor 
            J_B= 1 
        else: J_B= np.exp(-1.35*F_bp*(1-(2*r_s))**(1/3))
        
        #Shellside HT Coefficient 
        
        U_shellside= U_shellside_ideal*J_c*J_L*J_B 
        
        #PRESSURE DROP CALCULATION 
        
        if Re<=2300 : 
            K_f= 0.795+0.247*1000/Re+0.335*10000/Re**2-0.155*10000/Re**3+0.241*10000/Re**4
        else: K_f= 0.245+0.339*10000/Re-0.984*10000000/Re**2+0.133*10**11/Re**3-0.599*10**13/Re**4
        
        DeltaP_CI= N_c*K_f*0.5*rho_h*v_max**2 #Ideal Crossflow Pressure Drop
        
        N_cw= 0.8*L_c/L_tp #Number Effective Crossflow Rows 
        S_w= (D_s**2/4)*(np.arccos(1-2*L_c/D_s)-(1-2*L_c/D_s)*np.sqrt(1-(1-2*L_c/D_s)**2))-N_t*(1-F_c)*np.pi*d_o**2/8
        #S_w is the window flow area
        
        DeltaP_WI= (2+0.6*N_cw)*mdot_h**2/(2*S_m*S_w*rho_h) #IDEAL Window Zone Pressure Drop 
        
        S_b= L_bb*L_bc 
        if SS>0.5:
            R_B= 1 
        else: R_B= np.exp(-3.7*S_b*(1-r_s**(1/3))/S_m)
        
        R_L= np.exp(-1.33*(1+r_s))*r_lm**x 
        
        DeltaP_h= (((N-1)*DeltaP_CI*R_B+N*DeltaP_WI)*R_L+2*DeltaP_CI*R_B*(1+N_cw/N_c))/100000
        #PRESSURE DROP IN PASCALS or N/m^2 
        
        #Tubeside HT Coefficient 
        
        A_tube= 0.25*np.pi*d_i**2  #Tube Area in square meters 
        A_front=  A_tube*N_t/N_p   #Frontal Area in Square meters 
        v_tube= mdot_c/(A_front*N_t)  #velocity in tube in m/s
        Mdot= mdot_c/A_front    #Mass flux tube side in kg/s/m^2
        
        #Dimensionless Parameters 
        
        Re_tube= d_i*Mdot/mu_c  #Reynolds Number 
        Pr_tube= c_pc*mu_c/k_c  #Prandtl Number
        
        f= (0.79*np.log(Re_tube-1.64))**-2  #Friction Factor 
        
        if Re_tube<=2300: #LAMINAR
            Nu_tube= 4.36  #NUSSELT Number 
        if Re_tube>2300:   #TURBULENT 
            Nu_tube= (f/8)*(Re_tube-1000)*Pr_tube/(1+12.7*(f/8)**0.5*(Pr_tube**(2/3)-1))
        
        #Tubeside HT Coefficient 
        
        U_tubeside= Nu_tube*k_c/d_i #Tubeside HTC in W/m^2K
        
        #OVERALL HEAT TRANSFER COEFFICIENT
        
        U= ((1/U_shellside)+(d_o/d_i)/U_tubeside+(d_o*np.log(d_o/d_i)/(2*U_wall)))**-1  
        
        #Tubeside PRESSURE DROP 
        
        if Re_tube<=2300: 
            f_0= 16/Re_tube 
        if Re_tube>2300: 
            f_0= 0.079*Re_tube**-0.25 
        
        DeltaP_c= 4*f_0*L_tube*N_p*Mdot**2/(100000*2*d_i*rho_c)
        
        
        i= i+1 
            
    
    
    
    print('The overall heat transfer coefficient is', U,'W/m^2K')
    print('The shellside pressure drop is',DeltaP_h,'bar')
    print('The tubeside pressure drop is',DeltaP_c,'bar')
    print('The effectiveness is',R)

if CHOICE==2: 

    #PRINTED CIRCUIT HEAT EXCHANGER 
    
    
    #cHANNEL PROPERTIES 
    
    d_pche= 0.003 #channel diameter in meters 
    v_pche_guess= 1.5   #channel velocity in m/s 
    t_p= 0.00317  #ratio between channel diameter and plate thickness 
    p_pche= 1.1*d_pche #pitch between two horizontal channels in m 
    t_wall= 0.00317 #wall thickness in meters
    
    #FRONTAL AREA 
    
    Vdot_h= mdot_h/rho_h   #Volumetric flow rate for hot fluid m^3/s
    Vdot_c= mdot_c/rho_c   #Volumetric flow rate for cold fluid 
    A_front_h= Vdot_h/v_pche_guess  #Frontal Area of Hot Fluid in m^2
    A_front_c= Vdot_c/v_pche_guess  #Frontal Area of Cold Fluid 
    
    if A_front_h>A_front_c: #Conditional statement to choose higher frontal area 
        A_front_pche= A_front_h 
    else: A_front_pche= A_front_c 
    
    #New velocities
    
    v_pche_h= Vdot_h/A_front_pche 
    v_pche_c= Vdot_c/A_front_pche 
    
    #SIZING 
    
    A_channel= np.pi*(d_pche**2)/8 #Cross-sectional area of channel in m^2
    
    N_channels= 2*A_front_pche/A_channel #Total number of channels 
    
    Lx= 0.5*p_pche*N_channels  #Dimension in x-direction in meters
    Ly= 0.5*t_p*N_channels     #Dimenion in y-direction 
    
    #Dimensionless Parameters 
    
    P_w= np.pi*d_pche*0.5+d_pche  #Wetted Perimeter of channel in m
    D_h= 4*A_channel/P_w          #Hydraulic Diameter in m
    
    Re_h= (2*mdot_h/N_channels)*D_h/(A_channel*mu_h)
    Re_c= (2*mdot_c/N_channels)*D_h/(A_channel*mu_c)
    
    if Re_h<=2300: 
        Nu_h= 4.36 #Fully-Developed Laminar flow and uniform heat flux 
    elif Re_h>2300: 
        Nu_h= 0.023*(Re_h**0.8)*(Pr_h**0.3) #Fully-Developed Turbulent Flow cooling 
    
    if Re_c<=2300: 
        Nu_c= 4.36 
    else: Nu_c= 0.023*(Re_c**0.8)*(Pr_c**0.4) #Heating 
        
    
    #Heat Transfer Coefficients 
    
    U_h= Nu_h*k_h/D_h #Hot-side HTC in W/m^2K
    U_c= Nu_c*k_c/D_h #Cold-side HTC 
    
    k_pche= 23.9 #ALLOY 617
    f_bar =0.4499*((d_pche/t_p)**2)+(0.3403)*(d_pche/t_p)+1  #correction factor for wal interference
    U_w= f_bar*(k_pche/(t_wall))
    
    
    U= ((1/U_h)+(1/U_c)+(1/U_w))**-1 #OVERALL HTC
    
    #Pressure Drop 
    
    A_req= Qdot/(DeltaT_Mean*U) #HT Area 
    L_pche= A_req/((d_pche+np.pi*d_pche/2)*0.5*N_channels) #Length of HX 
    
    if Re_h<=2300:  #Friction Factors for each fluid 
        f_ih= 15.78/Re_h 
    if Re_h>2300: 
        f_ih= 0.478*Re_h**(-0.26)
    
    if Re_c<=2300: 
        f_ic= 15.78/Re_c
    if Re_c>2300: 
        f_ic= 0.478*Re_c**(-0.26)
        
    DeltaP_h= 2*f_ih*L_pche*rho_h*(v_pche_h**2)/(D_h*100000)  #pressure drop on hot-side in Pa
    DeltaP_c= 2*f_ic*L_pche*rho_c*(v_pche_c**2)/(D_h*100000)  #pressure drop on hot-side in Pa
    
    print('PCHE- The overall heat transfer coefficient is', U,'W/m^2K')
    print('PCHE- The hot-side pressure drop is',DeltaP_h,'bar')
    print('PCHE- The cold-side pressure drop is',DeltaP_c,'bar')

# Calculating the estimation of Equipment Cost
#  
#    C_b= np.exp(8.202+0.0156*np.log(A_req) + (0.0681*(np.log(A_req))**2)) #Base cost
#    F_d= np.exp(-0.903 + 0.906*np.log(A_req))  #Fixed Head Cost Factor
#    F_p= 0.895 + 0.4981*np.log(A_req) #Pressure Cost Factor
#    F_m= 1.4144 +  0.23296*np.log(A_req) #Material Cost Factor (SS316)
#    cost_in_2020= 740.0/567.7 # ratio between index of 2020 and index of 2015
#    cost_tot= C_b*F_d*F_p*F_m
#    
#    print ('the total area is ', A_req)

    
    #calculating efficiency
    
    
C_c=c_pc*mdot_c 
C_h=c_ph*mdot_h
C_min=min(C_c,C_h)
Cr=C_h/C_c
NTU= U*A_req/C_min
Fa=NTU*np.sqrt(1+Cr**2)/2 #This fin analogy for one shell pass only
eff=np.tanh(Fa)/Fa # efficiency

#pumping power calculation
P_=(1/eff)*((mdot_c*DeltaP_c/rho_c)+(mdot_h*DeltaP_h/rho_h))

# Calculating cost 

tau= 7920 # hours operating
K_el= 0.1319 #USD/kWh as of April 2020
cost_elec=K_el*tau*(P_/1000) #electricity cost


print('The efficiency is ',eff)
print('The pumping power required is ',P_,'W')
print ('The cost electricity is ',cost_elec, 'USD')
       
