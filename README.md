# AliAerospace Rocket Turbine Studio

<img width="1479" height="330" alt="image" src="https://github.com/user-attachments/assets/b94ab483-3567-4b7f-992e-ffa9dca7e37a" />





https://github.com/user-attachments/assets/e464f51e-7cfa-41ea-b1ed-9f3542678576




This software is a tool for the design and geometry generation of axial turbines for rocket turbopumps using MATLAB.

The studio comprises of meanline turbine calculations and velocity triangles with various vortex laws to create rotor and vaned stators. Entire 3D CAD geometries can be developed through providing inputs which range from blisk diameters to total inlet temperatures. 

The studio can be used alongside CFD with optimisation schemes to further optimise precise rotor and stator geometries. The ability to alter the shape of these blades is thanks to multiple parameterised curves for each blade. The methodologies employed here are novel and were developed by S. Sudhof [1] who proposed improvements upon L. Pritchards eleven parameter blade shaping model [2].

Please access the "Rocket_Turbine_Studio.pdf" in the files for a description of the program regarding its inputs and outputs.


## Features
- Thermodynamic calculations for nozzle efficiency, mechanical efficiency, and total-to-total turbine efficiency.
- Velocity triangle calculations for inlet and outlet of stator and rotors.
- Custom vortex law implementations from free vortex to exponential vortex.
- Aids initial design using the more intuitive eleven parameter rotor target generation (user inputs smaller set of rotor inputs and the program generates the rest) from Pritchard on which Sudhofs parameters are fitted.
- Direct Sudhof rotor and stator section generation.
- Root, mean, and tip section generation.
- Cylindrical projection and spanwise quadratic lofting.
- 3D rotor and stator visualisation.
- Full- or partial-admission nozzle arrangement.
- CAD export of rotors and stators as .IGES file (can be converted to .STEP in almost all CAD programs. Used this method as it was compatible with MATLAB.)


## Method

The studio performs:

1. `solveMeanline` calculates flow angles.
2. `generateRotorGeometry` creates three Pritchard rotor target curves and fits the corresponding Sudhof rotor inputs.
3. `generateStatorGeometry` creates three Sudhof stator curves using spanwise outlet angles.
4. Both geometry generators will project the sections on cylindrical stream surfaces.
5. Root, mean, and tip curves are lofted into 3D blade mantles.
6. Once run, the program produces plots of blade angle spanwise distributions, spanwise reaction, stator and rotor cross sections (including those with stagger visualised from root to tip), and a 3D visual rotor blisk with stators positioned upstream.
7. Can set the CAD generation to true for .IGES file outputs.

Choosing initial values for Sudhof's rotor parameters is difficult. To counteract this, the user can provide the simpler eleven required parameters to generate a Pritchard rotor (user inputs smaller set of rotor inputs and the program generates the rest), to which the Sudhof rotor is fitted using a custom algorithm which matches geometric, tangent, and curvature conditions. The output can sometimes be off on the order of millimetres, but this is fine as this rotor is just the initial shape. Once the Sudhof rotor is fitted, its input parameters are shown in MATLAB's Command Window, which can then be used to to supply the rotor geometry to a CFD program for further optimisation of the curves on leading edge, trailing edge, suction, and pressure surfaces. After CFD optimisation, the new Sudhof inputs can be fed back into the studio to generate a CAD geometry of the final rotor (or stator) shape. The convention used is shown below in Fig. 1 and 2, and more in-depth in "Rocket_Turbine_Studio.pdf" (please refer to it before using the studio).


<img width="1297" height="902" alt="image" src="https://github.com/user-attachments/assets/f0b4c24e-6c56-4944-82a3-a98adf4ac6f7" />
Figure 1: Rotor which comprises of various curves.


<img width="1075" height="837" alt="image" src="https://github.com/user-attachments/assets/a2e16d8a-8954-4f6f-a82b-90fc828de291" />

Figure 2: Stator which comprises of various curves.

## Coordinates and units
Refer to "Rocket_Turbine_Studio.pdf".


## Requirements

MATLAB toolboxes required:

- NURBS toolbox (NURBS Toolbox by D.M. Spink) .
- IGES toolbox (igesout).
  
## Installation

1. Clone the repository
2. Open MATLAB
3. Provide inputs in runAliAerospaceRocketTurbineStudio.m
4. Run runAliAerospaceRocketTurbineStudio.m


## CAD export

On line 11 in runAliAerospaceRocketTurbineStudio.m, set generate_IGES = true; for it to compute and save the .IGES rotor and stator files in the "CAD output" folder.


## Main files
| File | Purpose |
|---|---|
| `runAliAerospaceRocketTurbineStudio.m` | Main script. Inputs entered and studio is run |
| `solveMeanline.m` | Calculates meanline thermodynamics, velocity triangles, and spanwise radial equilibrium distributions |
| `createPritchardRotor.m` | Generates 2D Pritchard rotor target profile |
| `fitPritchardToSudhof.m` | Fits the Pritchard rotor profile to the Sudhof rotor parameters |
| `createSudhofRotor.m` | Generates 2D Sudhof rotor profile |
| `createSudhofStator.m` | Generates 2D Sudhof stator profile |
| `calculateStaggerAngle.m` | Calculates stagger angle from blade metal angles |
| `generateRotorGeometry.m` | Generates, projects, and lofts 3D rotor |
| `generateStatorGeometry.m` | Generates, projects, and lofts 3D stator |
| `exportBladeRowsIGES.m` | Exports rotor and stator blade rows as IGES surfaces |
| `CAD output/` | Stores the generated IGES files |

## Important
Do not mix metres and millimetres. The workflow expects "cycle.d_m" in millimetres while thermodynamic quantities such as pressure, temperature, and velocity are SI.

## Rotor and stator convention
Refer to "Rocket_Turbine_Studio.pdf" for nomenclature and other important descriptions.

<img width="746" height="1079" alt="image" src="https://github.com/user-attachments/assets/2da56a28-6f81-4fea-934c-49173476f468" />




## References

1. S. Südhof, *Development Techniques for Supersonic Turbines*, doctoral dissertation, Institut für Raumfahrtsysteme, University of Stuttgart, 2020.
2. L. Pritchard, “An Eleven Parameter Axial Turbine Airfoil Geometry Model,” ASME 1985 International Gas Turbine Conference and Exhibit, Paper 85-GT-219, 1985.
3. R. H. Aungier, *Turbine Aerodynamics: Axial-Flow and Radial-Inflow Turbine Design and Analysis*, ASME Press, 2006.
4. M. Paluszny, F. Tovar, and R. R. Patterson, “G2 Composite Cubic Bézier Curves,” *Journal of Computational and Applied Mathematics*, vol. 102, no. 1, pp. 49–71, 1999.
5. L. Piegl and W. Tiller, *The NURBS Book*, 2nd ed., Springer, 1997.
