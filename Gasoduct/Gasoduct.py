## Thermal Design Assignment 3
## By Noah Wilson
## 3/13/21

## IMPORTANT ##
## To run this code you must have the file "Compressor_Map.npy" in the same directory.  Otherwise
## you will need to change the filepath in line 37.

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.optimize import minimize_scalar
import time

class Gas_Properties:

    Pamb = 101353. # Ambient (minimum) pressure in Pa (14.7 psi)
    P = Pamb # Actual pressure in Pa (default isth e ambient pressure)
    T = 288.17 # Methane temp in K (518.7 R)
    mdot = 31.7515 # Mass flow rate in kg/s (70 lb/s)
    MW = 0.01604 # Molar mass of CH4 in kg/mol
    Rk = 1.07168 # Gas constant to isentropic coefficient ratio for methane to air at ambient conditions

    @property
    def rho(self):
        A,B,C,D,E,F,G = [-3.32545752e-38 ,3.06145459e-31, -1.14567539e-24,
                          2.27771128e-18, -2.77950237e-12, 5.31425467e-6, 1.76901594e-1] # Polynomial fit parameters for density
        return A*self.P**6 + B*self.P**5 + C*self.P**4 + D*self.P**3 + E*self.P**2 + F*self.P + G # Density in kg/m^3

    @property
    def mu(self):
        A,B,C = [ 4.93143680e-8, 1.10473563e-6, -8.41671669e-6] # Logarithmic fit parameters for viscosity
        return A*np.log(self.P)**2 + B*np.log(self.P) + C # Dynamic viscosity in Pa-s

Gas = Gas_Properties()

class Compressor_Properties:

    def __init__(self,comp_ratio=10):

        self.comp_ratio = comp_ratio # Compression ratio P_out/P_in (Default 10)
        self.Fixed_Cost = 5e6 # Fixed cost of one compressor in $

        X,Y,Effs = np.load('Compressor_Map.npy') # Loads compressor map data
        self.interp = LinearNDInterpolator(list(zip(X,Y)),Effs) # Interpolates compressor map data


    @property
    def Max_eff(self): # Finds optimum efficiency and "dimensionless" mass flow rate

        def Func(dimm):
            E = -float(self.interp(dimm,self.comp_ratio)) # Efficiency value from compressor map
            return {True:E,False:0}[E==E]

        return minimize_scalar(Func,method='Bounded',bounds=(30,120))

    @property
    def D_ratio(self): # Optimum ratio of compressor diameter to pipe diameter
        dimm = self.Max_eff.x # "Dimensionless" mass flow rate in lb/s
        # return round(np.sqrt(70/(self.comp_ratio*np.sqrt(Gas.Rk)*dimm)),3)
        return round(np.sqrt(70/(np.sqrt(Gas.Rk)*dimm)),3)

    @property
    def efficiency(self): # Optimum efficeincy of compressor
        eff = -round(self.Max_eff.fun,4)
        if eff <= 0:
            raise ValueError(f'Compression Ratio is out of operating range.')
        else:
            return eff

    @property
    def Power(self):
        A,B,C = [23499.80265741, 135474.18610104, 10487.99995938] # Logarithmic fit parameters for enthalpy change as a function of compression ratio
        dH = A*np.log(self.comp_ratio)**2 + B*np.log(self.comp_ratio) + C # Enthalpy change in J/kg
        return 0.001*Gas.mdot*dH/self.efficiency # Power requried for compressor in kW

    @property
    def yearly_cost(self):
        energy = self.Power * 8760 # Yearly energy required in kWh for one compressor
        elec_cost = 0.1 # Cost of electricity in $/kWh
        return elec_cost*energy

class Pipe_Properties(Gas_Properties):

    Materials = {'PVC'     :      [0,322.92], # [Roughness (m), Cost factor ($/m^2)]
                 'Steel'   :[4.57e-5,215.28], # [Roughness (m), Cost factor ($/m^2)]
                 'Concrete':[9.14e-4,107.64]} # [Roughness (m), Cost factor ($/m^2)]

    def __init__(self,material,D=1):

        self.material = material
        self.D = D # Diameter of pipe in m (Default is 1)
        self.L = 1.448e6 # Total length in m (900 miles)
        self.NComp = 1 # Number of compressors along the pipe
        self.ep = self.Materials[material][0] # Roughness in m
        self.cost_factor = self.Materials[material][1] # Cost factor in $/m^2
        self.Total_Cost = 0

    @property
    def Cost(self):
        return self.L * self.D * self.cost_factor # Pipe cost in $

    @property
    def head_loss(self):
        g = 9.81 # Gravitaional accleration (m/s^2)
        u = 4*self.mdot/(np.pi*self.rho*self.D**2) # Flow rate in the pipe (m/s)
        Re = self.rho*u*self.D/self.mu # Reynolds number
        f = .25*np.log10(((self.ep/self.D)/3.7) + (5.74/(Re**0.9)))**-2 # Darcy friction factor
        return (f/self.D)*((u**2)/(2*g)) # Head loss per unit length of the pipe (Pa/m)

Comp = Compressor_Properties() # Initialize compressor

def opt_comp_ratio(Cipt,Pipe): # Function used to optimize compresion ratio

    Comp.comp_ratio=Cipt # Initialize the compressor

    def opt_pipe_D(Dipt): # Function used to optimize pipe diameter for each compression ratio

        Pipe.D = Dipt # Set pipe diameter
        Pipe.P = Comp.comp_ratio * Pipe.Pamb # Set pipe pressure after first compression
        Pipe.Total_Cost = Pipe.Cost + Pipe.NComp*(15*Comp.yearly_cost + 5.0e6)
        dx = 100 # Distance step along the entire gasoduct in m
        L = 1.448e6 # Length of entire gasoduct in m
        Distance = np.arange(0,L,dx)
        ## The following loop subtracts the frictonal pressure loss for a section of pipe
        ## length dx and calculates the new head loss with that pressure and repeats.  When
        ## the pressure in the pipe reaches the ambient pressure the loop ends and the number of
        ## compresors needed are calculated by dividing the total length of pipe by the distance
        ## at which the pressure drops back to ambient pressure.
        for x in Distance:
            Pipe.P -= dx*Pipe.head_loss # Reduce pressure by the fricitonal loss
            if Pipe.P <= Pipe.Pamb: # Checks if pipe pressure is below ambient
                Pipe.NComp = int(L/x) + 1 # Calculates number of compressors needed
                Pipe.Total_Cost = Pipe.Cost + Pipe.NComp*(15*Comp.yearly_cost + 5.0e6)
                break
        return Pipe.Total_Cost

    minimize_scalar(opt_pipe_D,method='bounded',bounds=(.1,20)) # Optimizes diameter for each compression ratio
    return Pipe.Total_Cost

PVC = Pipe_Properties('PVC') # Initialize PVC pipe
Steel = Pipe_Properties('Steel') # Initialize steel pipe
Concrete = Pipe_Properties('Concrete') # Initialize concrete pipe

def main(Pipe):

    minimize_scalar(opt_comp_ratio,args=(Pipe),method='bounded',bounds=(5,25)) # Optimize compression ratio
    print(f'\nMaterial              : {Pipe.material}')
    print(f'Comp Ratio            : {Comp.comp_ratio}')
    print(f'Comp eff              : {Comp.efficiency}')
    print(f'D Ratio               : {Comp.D_ratio}')
    print(f'Pipe Diameter         : {Pipe.D}')
    print(f'Number of Compressors : {Pipe.NComp}')
    print(f'Total Cost            : {Pipe.Total_Cost}')
    print(f'Elec Cost             : {15*Comp.yearly_cost}')

if __name__ == '__main__':
    
    start = time.time()
    main(PVC) # Run the optimizatino for PVC
    main(Steel) # Run the optimization for Steel
    main(Concrete) # Run the optimization for concrete
    print(f'\n Computation time (sec) : {time.time() - start}')
