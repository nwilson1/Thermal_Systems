# Thermal_Systems

## Codes Included

### Regenerator.py

This code optimizes the variuos dimensions of a rotary regenerator implemented as part of a Brayton cycle.  The code can be run from the command line and requires no user input, but some parameters, such as the operating temperatures, can be modified in the code itself.  The CoolProp library is required to run this code.

Refer to the accompanying PDF for a more detailed explanation of the relevant physics and the overall numerical methodology.  It requires a minute or so to run the optimization, but all of the relevant parameters are shown in the output.  A sample output is shown below. 

~~~~
Efficiency : 41.4 %
Efficeincy without regenerator : 37.9 %
Carnot efficiency : 80 %
Turbine output pressure : 0.1 MPa
Compressor output pressure : 1.42 MPa
High temp : 1500 K
Low temp : 300 K
Wheel diameter : 2 m
Wheel width : 0.31 m
Rotational rate : 0.17 1/s
Half base width of channels : 1 mm
Thickness of wheel material : 0.1 mm
~~~~
### Gasoduct.py

This code estimates the minimum cost of a hypothetical 900 mile long methane pipeline.  The code examines three different pipe materials (PVC, concrete, steel) and determines which of these would requiere the least cost.  The analysis was based on the amount of frictional pressure loss each pipe would experience.  The more frictional loss, the more compressors are required to maintain pressure in the pipeline.  So the main parameters begin optimized were the diameter of the pipe, the compression ratio of the compressors, the number of compressors used.  

The code can be run from the command line and requires no user input, but some parameters can be modified in the code itself.  The numpy array file "Compressor_Map.npy" must be in the same directory as Gasoduct.py for the simulation to run.  The CoolProp library is also required to run this code.

Refer to the accompanying PDF for a more detailed explanation of the relevant physics and the overall numerical methodology.  A sample output is shown below. 

~~~~
Material              : PVC
Comp Ratio            : 13.251017855547635
Comp eff              : 0.85
D Ratio               : 0.899
Pipe Diameter         : 0.5222854982035631
Number of Compressors : 3
Total Cost            : 1021218262.501945
Elec Cost             : 254001249.1340859

Material              : Steel
Comp Ratio            : 19.37979950556771
Comp eff              : 0.8409
D Ratio               : 0.81
Pipe Diameter         : 0.45557909439321465
Number of Compressors : 2
Total Cost            : 765805717.5330757
Elec Cost             : 306895061.93927467

Material              : Concrete
Comp Ratio            : 10.341722956569452
Comp eff              : 0.84
D Ratio               : 0.956
Pipe Diameter         : 0.7040284389170407
Number of Compressors : 2
Total Cost            : 571950232.5104828
Elec Cost             : 226109222.5317595
~~~~
