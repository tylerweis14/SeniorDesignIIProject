import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import ode
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def dydt(t,y, params):
    x, y, z_f, z_c=y
    alpha,lamb,beta,c_pf,c_pc,m_f,m_c,W_ce,T_cine,a_f,n_e,alpha_f,alpha_c,h=params

    T_fe=T_cine+(1/(2*W_ce*c_pc)+(1/h))*a_f*n_e #equillibrium of fuel temp
    T_ce=T_cine+(a_f*n_e/(2*W_ce*c_pc)) #equillibrium of coolant temp
    u=(T_cine-T_cine)/T_cine
    w=(1300-W_ce)/W_ce
    Power = 1 #percentage
    p_c=(Power-(x))*10
    p=p_c+alpha_c*T_ce*z_c+alpha_f*T_fe*z_f
    
    
    dydt1 = -(beta*x/alpha)+(beta*y/alpha)+(p/alpha)+(p*x/alpha)
    dydt2 = (x-y)*lamb
    dydt3 = ((a_f*n_e*x)/(m_f*c_pf*T_fe))-(h*z_f/(m_f*c_pf))+(h*T_ce*z_c/(m_f*c_pf*T_fe))
    dydt4 = (h*T_fe*z_f/(m_c*c_pc*T_ce))-((2*c_pc*W_ce+h)*z_c/(m_c*c_pc))+((2*W_ce*T_cine*u)/(m_c*T_ce))
    -(2*W_ce*w*(T_ce-T_cine)/(m_c*T_ce))-(2*W_ce*w*z_c/m_c)+(2*W_ce*T_cine*u*w/(m_c*T_ce))
    
    derivs=[dydt1, dydt2, dydt3, dydt4]

    return derivs

def tempoutput(params2):
    c_pc,W_ce,T_cine,a_f,h,finaltempchange=params2
    T_fe=T_cine+(1/(2*W_ce*c_pc)+(1/h))*a_f*n_e #equillibrium of fuel temp
    T_ce=T_cine+(a_f*n_e/(2*W_ce*c_pc)) #equillibrium of coolant temp
    Tout = (T_fe-T_ce)/(Rf*W_ce*c_pc) + T_ce + finaltempchange
    power = W_ce*c_pc*(Tout-T_cine)
    return Tout, power/1e7



finaltemps = []
T_cine = 600

for j in range(10):
    
    alpha=0.001
    lamb=0.1
    beta=7.5*10**-3
    c_pf=717 #specific heat of graphite moderator
    c_pc=2414.7 #specific heat of FliBE
    m_f=470000*(1.5/1000) #mass of u235 in 470,000 pellets
    m_c=90830.8 #mass of coolant
    W_ce=1500 #mass flow rate
#    T_cine = newtemp
#    T_cine=600 #Temperature in
    a_f=7.0e6
    n_e=200.0
    alpha_f=-3.8e-5 #Change in reactivity based on temp of fuel
    alpha_c=-1.8e-5 #Change in reactivity based on temp for moderator
    h=4700*1940 #heat transfer coefficient and total area of fuel
    Rf = .0005 #fouling factor
    
    
    params=[alpha,lamb,beta,c_pf,c_pc,m_f,m_c,W_ce,T_cine,a_f,n_e,alpha_f,alpha_c,h]
    x0=0.0 #starting neutron pop
    y0=0.0 #starting precursors
    z_f0=1.0 #starting fuel temp
    z_c0=1.0 #starting moderator temp
    y0=[x0,y0,z_f0,z_c0]
    t0=0
    A=[]
    
    # Solver
    r = ode(dydt).set_integrator('dopri5', method='nsteps')
    r.set_initial_value(y0, t0).set_f_params(params)
    
    t1 =3600.0
    dt = 0.1
    T=[]
    while r.successful() and r.t < t1:
        r.integrate(r.t+dt)
        T=np.append(T,r.t)
        A=np.append(A,r.y)
    #print np.size(A)
    B= A.reshape(np.size(T),4)
    
    finaltempchange = sum(B[:,3])
    
    
        
    params2 = [c_pc,W_ce,T_cine,a_f,h,finaltempchange]    
    
    
    finaltemps.append(tempoutput(params2)[0])
    
    T_h_in= finaltemps[j] #ALL TEMPERATURES LISTED IN KELVIN, K 
    T_c_in= 300
    T_c_out= 600
    
    
    #FLUID PROPERTIES -- FLiBe (shellside hot)
    
    rho_h= 2518-0.406*(T_h_in+T_c_in)/2 #density in kg/m^3
    mu_h= 0.000116*np.exp(3755/((T_h_in+T_c_in)/2)) #viscosity in Pa*s
    c_ph= 2415.78  #specific heat in J/kg*K
    k_h= 0.629697+0.0005*((T_h_in+T_c_in)/2) #thermal conductivity in W/mK
    
    #FLUID PROPERTIES -- Solar Salt (tubeside cold)
    
    rho_c= 1804
    mu_c= 0.00169
    c_pc= 1520 
    k_c= 0.530 
    
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
    T_cine = T_h_out
  
plt.plot(finaltemps,[1,2,3,4,5,6,7,8,9,10])
plt.xlabel('Temperature')
plt.ylabel('Cycles')
plt.title('Temperature Vs Cycles')