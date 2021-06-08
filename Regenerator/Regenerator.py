## By Noah Wilson
## 2/23/20

import CoolProp.CoolProp as cp
from scipy.optimize import fsolve,minimize
import numpy as np
import matplotlib.pyplot as plt

class AirState:
    def __init__(self,**ipts):

        props,vals = zip(*ipts.items()) # Separate tuples for propertes and values
        propval = [props[0],vals[0],props[1],vals[1]] # Single list used for PropsSI input
        self.T = cp.PropsSI('T',*propval,'Air') # Temperature in K
        self.p = cp.PropsSI('P',*propval,'Air') # Pressure in Pa
        self.h = cp.PropsSI('HMASS',*propval,'Air') # Specific enthalpy in J/kg
        self.s = cp.PropsSI('SMASS',*propval,'Air') # Specific entropy in J/kg-K
        self.k = cp.PropsSI('CONDUCTIVITY',*propval,'Air') # Thermal conductivity in W/m-K
        self.Cp = cp.PropsSI('CPMASS',*propval,'Air') # Isobaric heat capacity in J/kg-K
        self.mu = cp.PropsSI('VISCOSITY',*propval,'Air') # Dynamic viscosity in Pa-s
        self.rho = cp.PropsSI('DMASS',*propval,'Air') # Density in kg/m^3
        self.Pr = cp.PropsSI('PRANDTL',*propval,'Air') # Density in kg/m^3

class RotaryRegenerator:

    def channels(self,State,mdot,Dh,Av,As,Af,L):

        Re = mdot*Dh/(State.mu*0.5*Av) # Reynolds number
        if Re >= 2300:
            f = (1.58*np.log(Re) - 3.28)**(-2) # Friction factor for turbulent flow
            Nu = (.5*f*(Re-1000)*State.Pr)/(1+12.7*((State.Pr**(2/3))-1)*np.sqrt(.5*f)) # Nusselt number for turbulent flow
        else:
            f = 16/Re # Friction factor for laminar flow
            Nu = 3.11 # Nusselt number for laminar flow in triangular tubes (From Incropera and Dewitt 2007)
        H = State.k*Nu/Dh # Convective heat transfter coefficient
        NTU =  H*As/(mdot*State.Cp) # Number of transfer units
        dP = (f*L*mdot**2)/(2*State.rho*Dh*(.5*Af)**2) # Pressure drop across channels
        return H,NTU,dP

    def effectiveness(self,C_h,C_c,NTU_h,NTU_c,M,C_r,Omega):
        C_min, C_max = min(C_h,C_c),max(C_h,C_c) # Max and min heat capacity rates in W/K
        NTU_0 = (((1/(C_h*NTU_h))+(1/(C_c*NTU_c)))**-1)/C_min # Combined hot flow and cold flow NTU
        Cr_st = M*C_r*2*Omega/C_min # Normalized heat capacity rate for matrix material
        eps_r = (NTU_0/(1+NTU_0))*(1-(1/(9*Cr_st**1.93))) # Effectiveness of the wheel
        if eps_r <= 0: # eps_r has a chance of being negative at very small Omega
            return 0 # Return zero effectiveness if calculated effectiveness is less than zero
        else:
            return eps_r

    def T_out(self,mh,mc,Hot,Cold,epsilon,M,C_r,Omega):

        ## Equations for to find the temperature of the hot and cold side of
        ## the wheel assuming the output gas and the wheel are in thermal equilibrium
        ## on both sides.
        def Wheel_Temps(ipts):
            Twh,Twc = ipts[0],ipts[1]
            Q_max = M*C_r*2*Omega*(Twh - Twc)
            Out = [0,0]
            Out[0] = Hot.T - Twh - (Q_max/(mh*Hot.Cp))
            Out[1] = Cold.T - Twc + (Q_max/(mc*Cold.Cp))
            return Out
        Twh,Twc = fsolve(Wheel_Temps,[0,0]) # Equilibrium temperatures of the wheel in K
        Q_max = M*C_r*2*Omega*(Twh - Twc) # Maximum heat transfer rate to or from the wheel in W
        Q = epsilon*Q_max # Actual heat rate to or from the wheel in W
        T_c_out = Hot.T - Q/(mh*Hot.Cp) # Output temperature on hot side
        T_h_out = Cold.T + Q/(mc*Cold.Cp) # Output temperature on cold side
        return T_h_out,T_c_out

    ## Inputs:
    ##      Hot - State of incoming air from hot side
    ##      Cold - State of incoming air from cold side
    ##      m_dot - Mass flow rate, assumed to be equal on hot and cold side in kg/s
    ##      a - Half base of one channel in m
    ##      t - Thickness of material in m
    ##      D - Diameter of wheel in m
    ##      L - Width of wheel in m
    ##      rho_r - Density of wheel material in kg/m^3
    ##      C_r - Specific heat of wheel material in J/kg-K
    ##      Omega - Rotational rate of wheel in 1/s
    def __init__(self,Cold,Hot,m_dot,a,t,D,L,rho_r,C_r,Omega):

        A_ch = np.sqrt(3)*(a**2) # Area of a single channel in m^2
        Per = 6*a # Wetted perimeter of each channel in m
        D_h = 4*A_ch/Per # Hydraulic diameter 4*Area of channel/Wetted perimeter in m
        A_fr = (np.pi*D**2)/4 # Frontal area exposed to flow in m^2
        N_ch = A_fr/(A_ch+3*a*t) # Number of channels in the regenerator
        A_void = N_ch*A_ch # Void area of front of wheel in m^2
        A_surf = Per*L*N_ch # Total Heat transfer surface area in m^2
        M_r = rho_r*(A_fr - A_void)*L # Mass of the matrix material in kg
        mleak = A_void*L*Omega*(Cold.rho - Hot.rho) # Rate of mass leakage between the two streams in kg/s
        m_dot_h,m_dot_c = m_dot,m_dot+mleak # Mass flow rate for hot and cold streams respectively in kg/s
        Ch,Cc = m_dot_h*Hot.Cp, m_dot_c*Cold.Cp # Heat capcity rates for hot flow and cold flow in W/K
        self.Hh,self.NTUh,self.dPh = self.channels(Hot,m_dot_h,D_h,A_void,A_surf,A_fr,L) # Hot to cold stream
        self.Hc,self.NTUc,self.dPc = self.channels(Cold,m_dot_c,D_h,A_void,A_surf,A_fr,L) # Hot to cold stream
        self.epsilon = self.effectiveness(Ch,Cc,self.NTUh,self.NTUc,M_r,C_r,Omega) # Effectiveness of the wheel
        self.Th_out,self.Tc_out = self.T_out(m_dot_h,m_dot_c,Hot,Cold,self.epsilon,M_r,C_r,Omega) # Output gas temperatures in K

def Cycle(Phigh,a,t,L,omega):

    ## Optimizable properties
    ## Phigh - Higher pressure after compressor in Pa
    ## a - Half of one channel base in m
    ## t - Thickness of channel material in m
    ## omega - Rotational rate of the regenerator in 1/s

    ## Cycle properties
    Wout = 10e6 # Work output in W (10 MW)
    Th = 1500 # Temp after combustion in K
    Tc = 300 # Inlet temp in K
    Plow = 0.1e6 # Lower pressure at inlet and outlet in Pa (1 bar)
    eta_comp = 0.85 # Compressor efficiency
    eta_turb = 0.9 # Turbine efficiency

    ## Regenerator properties
    D = 2 # Diameter of wheel in m
    rhor = 8135 # Density of material in kg/m^3
    Cr = 537.6 # Heat capacity of material in J/kg-K

    dPh,dPc = 0,0 # Intial values for pressure drop across regenerator
    old_eta,eta = 100 ,10 # Values to initialize while loop

    ## While loop that iterates the system incorporating pressure drop from the
    ## regenerator until the calculated efficeincy of the cycle does not change from
    ## the previous iteration by more than 0.01%
    while abs(1-(eta/old_eta)) > 1e-4:

        S1 = AirState(T=Tc,P=Plow-dPh) # State 1 : Inlet to compressor

        h2s = cp.PropsSI('H','P',Phigh,'S',S1.s,'Air') # Enthalpy for ideal compressor
        h2 = S1.h - (S1.h - h2s)/eta_comp # Real compressor enthalpy
        S2 = AirState(P=Phigh,H=h2) # State 2 : Output from compressor, cold input to regenerator

        S4 = AirState(P=Phigh-dPc,T=Th) # State 4 : Post-combustion, turbine intake

        h5s = cp.PropsSI('H','P',Plow,'S',S4.s,'Air') # Enthalpy for ideal turbine
        h5 = S4.h - eta_turb*(S4.h - h5s) # Real turbine enthalpy
        S5 = AirState(P=Plow,H=h5) #State 5 : Turbine exhaust, warm input to regenerator

        mdot = Wout/(S4.h - S5.h) # mass flow rate in kg/s

        Regen = RotaryRegenerator(S2,S5,mdot,a,t,D,L,rhor,Cr,omega)
        dPh,dPc = Regen.dPh,Regen.dPc # Pressure drop across each flow of the regenerator

        S3 = AirState(P=Phigh,T=Regen.Th_out) # State 3 : Warm output from regenerator, combustion input
        S6 = AirState(P=Plow,T=Regen.Tc_out) # State 6 : Cold output from regenerator

        Win = mdot*(S2.h - S1.h) # Work input to compressor
        Qin = mdot*(S4.h - S3.h) # Heat input from combustor

        old_eta = eta
        eta = (Wout - Win)/Qin # Thermal efficiency of the cycle
    Q_woRegen = mdot*(S4.h-S2.h) # Heat input without regenerator
    eta_woRegen = (Wout - Win)/Q_woRegen # Efficiency without regenerator
    return eta,eta_woRegen

Opt_PLW = lambda ipts : -Cycle(Phigh = ipts[0]*1e6,a = 0.001,t = 0.0001,L=ipts[1],omega = ipts[2])[0]
P,L,W = minimize(Opt_PLW,[1,0.5,0.25]).x
eta,eta_woRegen = Cycle(Phigh = P*1e6,a = 0.001,t = 0.0001,L=L,omega = W)

print(f'\nEfficiency : {100*round(eta,3)} %')
print(f'Efficeincy without regenerator : {100*round(eta_woRegen,3)} %')
print('Carnot efficiency : 80 %')
print('Turbine output pressure : 0.1 MPa')
print(f'Compressor output pressure : {round(P,2)} MPa')
print('High temp : 1500 K')
print('Low temp : 300 K')
print('Wheel diameter : 2 m')
print(f'Wheel width : {round(L,2)} m')
print(f'Rotational rate : {round(W,2)} 1/s')
print('Half base width of channels : 1 mm')
print('Thickness of wheel material : 0.1 mm\n')
