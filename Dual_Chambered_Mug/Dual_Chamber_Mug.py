## Code by Noah Wilson
## 5/12/21
## This code simulates the physical and thermal properties of the
## dual-chambered travel mug concept developed by Abigail Sickbert,
## Robert Smith, and myslef.

import numpy as  np
import matplotlib.pyplot as plt

## General thermodynamic properties
sig = 5.67e-8 # Stephan-Boltzmann Constant W/m^2-K^4
C = 4.186e3 # Specific heat of water in J/kg-K
rho = 1000 # Density of water in kg/m^3
Tamb = 293.15 # Ambient tamperatrure (20 C) in K
T_target = 333.15 # Target drinking temperature (60 C ~ 140 F) in K
rhost = 7930 # Density of 18/8 steel in kg/m^3
kst = 16.2 # Thermal conductivity of 18/8 stainless steel in W/m-K

class bottle:
    def __init__(self,Cap,t,d,H,h,F):

        self.Cap = Cap
        self.t = t # Steel thickness
        self.tc = 4*t # Cooling chamber wall thickness
        self.d = d
        self.H = H
        self.h = h
        self.h_lid = 0.025 # Lid thickness in m
        self.t_lid = 0.005 # Width of lid threading
        self.h_sep = 0.005 # Thickness of separator b/t chambers
        self.t_sep = 0.005 # Width of threading to connect cooling chamber to insulating chamber

        ## Thermodynamic properties of the bottle
        self.ep = 0.02 # Emissivity of polished steel
        self.rhomatH = rhost
        self.k_matH = kst
        self.rhomatC = 857 # Density of Polyethelene
        self.k_matC = 0.4  # Thermal conductivity of polyproelene
        self.rhopoly = 920 # Density of polypropelene in kg/m^3
        self.k_poly = 0.11 # Thermal conductivity of polypropelene in W/m-K
        self.F = F # Fraction of the cooling chamber filled
        self.T_H = 373.15 # Temperature of water in the insulating chamber in K (Initial is boiling)
        self.T_C = int(self.F>0)*self.T_H + int(self.F==0)*Tamb # Temperature of the water in cooling chamber in K (Initial is equal to water in insulating chamber or ambient temperature if empty)
        self.T_out = Tamb

    @property
    def r(self):
        ## Inner radius of inner wall of the insulating chamber in m
        return np.sqrt(self.Cap/(np.pi*(self.H - self.t_sep)))

    @property
    def R(self):
        ## Inner radius of outer wall of the insulating chamber in m
        return self.r + self.t + self.d

    @property
    def A(self):
        ## Outer surface area of inner wall of the insulating chamber
        return 2*np.pi*(self.r + self.t)*(self.H - self.d) +  np.pi*(self.r + self.t)**2

    @property
    def Cap_C(self):
        ## Capacity of the cooling chamber in m^3
        return (np.pi*self.R**2)*(self.h-self.t_lid)

    @property
    def m_H(self):
        ## Mass of water being kept hot in kg
        return rho*(self.Cap - self.F*self.Cap_C)

    @property
    def m_C(self):
        ## Mass of water being cooled in kg
        return self.F*rho*self.Cap_C

    @property
    def CM(self):
        ## Finds height of the center of mass

        R = self.R
        r = self.r

        ## Insulating chamber
        m_base = self.rhomatH*self.t*np.pi*R**2 # Mass of base of the bottle
        m_outer = self.rhomatH*np.pi*((R + self.t)**2 - R**2)*self.H # Mass of outer wall
        m_rim = self.rhomatH*np.pi*(R**2 - (r +self.t)**2)*self.t # mass of rim at top of insulating chamber
        m_inner = self.rhomatH*np.pi*((r + self.t)**2 - r**2)*(self.H - self.t - self.d) # Mass of inner wall
        m_bot = self.rhomatH*self.t*np.pi*r**2 # Mass of the botoom of the inside of the insulating chamber
        m_hot = m_base + m_outer + m_rim + m_inner + m_bot # Total mass of the insulating chamber
        CM_hot = (m_base*(self.t/2) + m_outer*(self.H/2) + m_rim*(self.H - self.t/2) + \
        m_inner*((self.H + self.t + self.d)/2) + m_bot*(3*self.t/2 + self.d))/m_hot # Center of mass of the insulating chamber

        ## Polypropelene separator
        m_hsep = self.rhopoly*self.h_sep*np.pi*R**2 # Mass of non-threaded part of separator
        m_tsep = self.rhopoly*self.t_sep*np.pi*r**2 # Mass of threaded part of separator
        m_sep = m_hsep + m_tsep # Total mass of the separator
        CM_sep = (m_hsep*(self.H + self.h_sep/2) + m_tsep*(self.H - self.t_sep/2))/m_sep # CM of the separator

        ## Cooling chamber
        m_cool = self.rhomatC*self.h*np.pi*(R**2 - (R - self.tc)**2) # Mass of steel in the cooling chamber
        CM_cool = self.H + self.h_sep + self.h/2 # Center of mass of the cooling chamber

        m_H,m_C = self.m_H,self.m_C
        m_h2o = m_H + m_C # Total mass of water
        CM_h2o = (m_C*(self.H + self.h_sep + self.F*self.h/2) + m_H*(2*self.t + self.d + m_H/(2*rho*np.pi*r**2)))/m_h2o # Center of mass of water

        ## Lid
        m_hlid = self.rhopoly*self.h_lid*np.pi*R**2 # Mass of main part of the Lid
        m_tlid = self.rhopoly*self.t_lid*np.pi*(R - self.tc)**2 # Mass of threaded part of the lid
        m_lid = 0.75*(m_hlid + m_tlid) # Total mass of the lid
        CM_lid = (m_hlid*((self.h_lid/2) + self.H + self.h_sep + self.h) + m_tlid*(self.H + self.h_sep + self.h - (self.t_lid/2)))/m_lid

        m_tot = m_hot + m_sep + m_cool + m_h2o + m_lid
        return (m_hot*CM_hot + m_sep*CM_sep + m_cool*CM_cool + m_h2o*CM_h2o + m_lid*CM_lid)/m_tot # Total center of mass
    @property
    def mtot(self):
        ## Finds height of the center of mass

        R = self.R
        r = self.r

        ## Insulating chamber
        m_base = self.rhomatH*self.t*np.pi*R**2 # Mass of base of the bottle
        m_outer = self.rhomatH*np.pi*((R + self.t)**2 - R**2)*self.H # Mass of outer wall
        m_rim = self.rhomatH*np.pi*(R**2 - (r +self.t)**2)*self.t # mass of rim at top of insulating chamber
        m_inner = self.rhomatH*np.pi*((r + self.t)**2 - r**2)*(self.H - self.t - self.d) # Mass of inner wall
        m_bot = self.rhomatH*self.t*np.pi*r**2 # Mass of the botoom of the inside of the insulating chamber
        m_hot = m_base + m_outer + m_rim + m_inner + m_bot # Total mass of the insulating chamber

        ## Polypropelene separator
        m_hsep = self.rhopoly*self.h_sep*np.pi*R**2 # Mass of non-threaded part of separator
        m_tsep = self.rhopoly*self.t_sep*np.pi*r**2 # Mass of threaded part of separator
        m_sep = m_hsep + m_tsep # Total mass of the separator separator

        ## Cooling chamber
        m_cool = self.rhomatC*self.h*np.pi*(R**2 - (R - self.tc)**2) # Mass of steel in the cooling chamber
        ## Lid
        m_hlid = self.rhopoly*self.h_lid*np.pi*R**2 # Mass of main part of the Lid
        m_tlid = self.rhopoly*self.t_lid*np.pi*(R - self.tc)**2 # Mass of threaded part of the lid
        m_lid = 0.75*(m_hlid + m_tlid) # Total mass of the lid

        m_tot = m_hot + m_sep + m_cool + m_lid
        return np.array([m_hot, m_sep,m_cool,m_lid])

    @property
    def tip(self):
        # Tipping angle of the bottle in degrees
        return 180*np.arctan(self.R/self.CM)/np.pi

    @property
    def Q_rad_w2w(self):
        ## Radiative heat transfer between the walls of the insulating chamber in W
        return sig*self.ep*self.A*((self.T_H**4) - (Tamb**4))

    @property
    def Q_cond_w2w(self):
        ## Conductive heat transfer across the connection betwen the insulating walls in W.
        ## Treats the connecting between the walls as a flat ring.

        return 1/(2*np.pi*self.t*self.k_matH/np.log(self.R/self.r))

    @property
    def Q_conv_out(self):
        ## Free convective cooling from the vertical  walls of the cooling chamber

        Pr = 0.708 # Prandtl number of air at STP

        v = 1.821e-05# Kinematic viscosity of air at STP in Pa-s
        B = 1/Tamb # Volumetric thermal expansion for ideal gas
        g = 9.81 # Acceleration due to gravity
        Ts = self.T_out # Surface temperature
        L = self.H # Vertical wall length
        Gr = (g*B*(Ts - Tamb)*L**3)/(v**2) # Grashof number

        Ra = Gr*Pr # Rayleigh number

        Nu = 0.68 + (0.67*Ra**0.25)/((1 + (0.492/Pr)**(9/16))**(4/9)) # Nusselt number for a vertical plate

        k = 0.0259 # Thermal conductivity of air at STP in W/m-K

        h = Nu*k/L # Convective heat transfer coefficient

        A = 2*np.pi*self.R*self.h

        Rtot = (1/(h*A))
        return Rtot

    @property
    def Q_cond_sep(self):
        ## Conductive heat transfer from water through the separator
        ## Only applies if the insulating chamber is full

        Asep = np.pi*self.r**2
        Lsep = self.h_sep + self.t_sep
        return int(self.F==0)/(self.k_poly*(Asep/Lsep))

    @property
    def Q_cond_air(self):
        ## Heat conducted through the air then the separator in the insulating chamber if it is paritally full
        k = 0.033 # Thermal conductivity of air at 100 C and 1 atm
        A = np.pi*self.r**2
        L = A*(self.Cap - (self.m_H/rho))
        Asep = np.pi*self.r**2
        Lsep = self.h_sep + self.t_sep
        Rtot = (L/(k*A)) + Lsep/(self.k_poly*Asep)
        return int(self.F>0)*Rtot


    @property
    def Q_wconv(self):
        ## Free convective cooling from the vertical  walls of the cooling chamber

        Pr = 0.708 # Prandtl number of air at STP

        v = 1.821e-05# Kinematic viscosity of air at STP in Pa-s
        B = 1/Tamb # Volumetric thermal expansion for ideal gas
        g = 9.81 # Acceleration due to gravity
        Ts = self.T_C # Surface temperature
        L = self.h # Vertical wall length
        Gr = (g*B*(Ts - Tamb)*L**3)/(v**2) # Grashof number

        Ra = Gr*Pr # Rayleigh number

        Nu = 0.68 + (0.67*Ra**0.25)/((1 + (0.492/Pr)**(9/16))**(4/9)) # Nusselt number for a vertical plate

        k = 0.0259 # Thermal conductivity of air at STP in W/m-K

        h = Nu*k/L # Convective heat transfer coefficient

        A = 2*np.pi*self.R*self.h

        Rtot = (1/(h*A)) + self.tc/(self.k_matC*A)

        return Rtot

    @property
    def Q_lconv(self):
        ## Free convective cooling from the top of the lid

        Pr = 0.708 # Prandtl number of air at STP

        v = 1.821e-05# Kinematic viscosity of air at STP in Pa-s
        B = 1/Tamb # Volumetric thermal expansion for ideal gas
        g = 9.81 # Acceleration due to gravity
        Ts = self.T_C # Surface temperature
        L = self.R/2 # Characteristic length of the lid
        Gr = (g*B*(Ts - Tamb)*L**3)/(v**2) # Grashof number

        Ra = Gr*Pr # Rayleigh number

        Nu = 0.54*Ra**(0.25)

        k = 0.0259 # Thermal conductivity of air at STP in W/m-K

        h = Nu*k/L # Convective heat transfer coefficient

        A = np.pi*(2*L)**2

        Rtot = (1/(h*A)) + (self.h_lid + self.t_lid)/(0.75*self.k_poly*A)

        return Rtot

    @property
    def Q_tot_H(self):
        return self.Q_rad_w2w + (self.T_H - Tamb)/(self.Q_cond_w2w + self.Q_cond_sep + self.Q_cond_air + self.Q_conv_out)
        # return self.Q_rad_w2w + self.Q_cond_w2w + self.Q_cond_sep + self.Q_cond_air + self.Q_conv_out

    @property
    def Q_tot_C(self):
        ## Total heat transfer from the cooling chamber
        return (self.T_C - Tamb)*(self.Q_wconv + self.Q_lconv) - (self.T_H - self.T_C)/self.Q_cond_air

    @property
    def dTdt_C(self):
        if self.m_C > 0:
            return self.Q_tot_C/(self.m_C * C) # Rate of temperature change of the cooling water
        else:
            return 0

    @property
    def dTdt_H(self):
        return (self.Q_tot_H)/(self.m_H*C) # Rate of temperature change of insulated water in K/s

    @property
    def dTdt_st(self):
        Cst = 500 # Heat capacity of steel in J/kg-K
        return (self.Q_tot_H)/(self.mtot[0]*Cst) # Rate of temperature change of the outer wall

def cool_time(Height,F):

    # Capacity = 0.000946353 # Bottle capacity in m^3 (~32 oz)
    Capacity = 6e-4 #  capacity in m^3 (~20 oz)
    Thickness = 0.0005 # Wall thickness in m
    Gap = 0.0025 # Insulating gap in m
    H = Height # Height of the insulating chamber in m
    h = 0.05 # Height of the cooling chamber in m

    bot = bottle(Capacity,Thickness,Gap,H,h,F)
    dt = 60 # Timestep in seconds
    t = 0 # Initialize time
    times = []
    TempH = []
    TempC = []

    while bot.T_H > 300:

        times.append(t/3600)
        TempH.append(bot.T_H-273)
        TempC.append(bot.T_C-273)
        t += dt # Increment timestep
        dT_C = dt*bot.dTdt_C # Change in cooling chamber temperature
        dT_H = dt*bot.dTdt_H # Change in hot chamber temperature
        dT_out = dt*bot.dTdt_st # Change in outer wall temperature
        bot.T_C -= dT_C
        bot.T_H -= dT_H
        bot.T_out += dT_out

    TempH = np.array(TempH)
    TempC = np.array(TempC)
    time2coolH = times[np.argmax(TempH<=60)] # Find the time at which the hot chamber reached 60 C
    time2coolC = times[np.argmax(TempC<=60)] # Find the time at which the hot chamber reached 60 C

    return time2coolH,time2coolC,bot

def main():
    th,tc,Bot = cool_time(0.15,1)
    #
    print(f'Insulating chamber height : {100*Bot.H} cm')
    print(f'Cooling chamber height : {100*Bot.h} cm')
    print(f'Outer bottle diameter : {200*(Bot.R+Bot.t)} cm')
    print(f'Inner diameter of insulating chamber : {200*Bot.r} cm')
    print(f'Thickness of steel : {100*Bot.t} cm')
    print(f'Insulating gap : {100*Bot.d} cm')
    print(f'Total height with cap and separator : {100*(Bot.h+Bot.H+Bot.h_sep+Bot.h_lid)} cm')
    print(f'Mass of bottle : {Bot.mtot} kg')
    print(f'Insulating capacity : {33814*Bot.Cap} oz')
    print(f'Cooling capacity : {33814*Bot.Cap_C} oz')
    print(fr'Time for insulating chamber to reach $60^\circ C$ when full: {th} hrs')
    print(fr'Time for cooling chamber to reach $60^\circ C$ when full: {tc} hrs')
    print(f'Tipping point when cooling chamber is full : {Bot.tip} Deg')


    Heights = np.linspace(0.1,0.3,100) # Range of heights
    Fs = 0*np.ones(100) # Set cooling chamber to empty
    timesHemp = []
    timesHfull = []
    timesC = []
    Capc = []
    Ds = []
    mass = []
    tipping_empty = []
    tipping_full = []

    for i in range(100):
        timeH,timeC,Bot = cool_time(Heights[i],Fs[i])
        timesHemp.append(timeH)
        Capc.append(33814*Bot.Cap_C)
        Ds.append(200*Bot.R)
        mass.append(1000*Bot.mtot)
        tipping_empty.append(Bot.tip)

    Fs = 1*np.ones(100)# Set cooling chamber to Full

    for i in range(100):
        timeH,timeC,Bot = cool_time(Heights[i],Fs[i])
        timesHfull.append(timeH)
        timesC.append(timeC)
        tipping_full.append(Bot.tip)

    fig1,axs1 = plt.subplots()
    mass = np.array(mass)
    axs1.plot(Heights,mass[:,0])
    axs1.plot(Heights,mass[:,2])
    axs1.plot(Heights,mass[:,1])
    axs1.plot(Heights,mass[:,3])
    axs1.legend(['Hot Chamber','Cooling Chamber','Separator','Lid'])
    axs1.set_xlabel('Bottle Height (cm)')
    axs1.set_ylabel('Mass (g)')
    axs1.set_title('Mass of Individual Components')

    fig2,axs2 = plt.subplots(2,1,sharex=True,figsize=(6,6))
    axs2[0].plot(100*Heights,timesHemp,100*Heights,timesHfull)
    axs2[1].plot(100*Heights,timesC)
    axs2[0].set_ylabel('Time (hrs)')
    axs2[0].set_title('Hot Chamber')
    axs2[0].legend(['Cooling Chamber Empty','Cooling Chamber Full'])
    axs2[1].set_xlabel('Bottle Height (cm)')
    axs2[1].set_title('Cooling Chamber')
    axs2[1].set_ylabel('Time (hrs)')
    plt.suptitle(r'Time to cool to $60^\circ$ C')

    fig3,axs3 = plt.subplots(2,2,sharex=True,figsize=(8,5))

    axs3[0,0].plot(100*Heights,mass)
    axs3[1,0].plot(100*Heights,Ds,'C1')
    axs3[0,1].plot(100*Heights,Capc,'C2')
    axs3[1,1].plot(100*Heights,tipping_empty,'C3',100*Heights,tipping_full,'C4')

    axs3[0,0].legend([r'Time to cool to $60^\circ$C'])
    axs3[0,0].set_ylabel('Time (hrs)')
    axs3[0,0].legend([r'Total bottle mass sans water'])
    axs3[0,0].set_ylabel('Mass (g)')
    axs3[1,0].legend(['Bottle diameters'])
    axs3[1,0].set_ylabel('D (cm)')
    axs3[1,0].set_xlabel('Bottle Height (cm)')
    axs3[0,1].set_ylabel('Capacity (oz)')
    axs3[0,1].legend(['Cooling Chamber Capacity'])
    axs3[0,1].yaxis.tick_right()
    axs3[0,1].yaxis.set_label_position("right")
    axs3[1,1].set_ylabel('Tip Angle (Deg)')
    axs3[1,1].set_xlabel('Bottle Height (cm)')
    axs3[1,1].legend(['Empty Cooling Chamber','Full Cooling Chamber'])
    axs3[1,1].yaxis.tick_right()
    axs3[1,1].yaxis.set_label_position("right")
    fig.subplots_adjust(wspace=.1)
    fig.suptitle('Various Parameters 20 oz Bottle')
    plt.show()

if __name__ == '__main__':
    main()
