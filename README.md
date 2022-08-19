# CD-Nozzle
Shape optimization study of a rocket nozzle using MATLAB


Nozzles are aerospace components found in propulsion devices such as turbojet engines or rocket engines. Its purpose is to accelerate high-pressure, high-temperature gases as much as possible to get optimal thrust for vehicle propulsion, according to Newton's 3rd law.

In this project, we are going to simulate 3 different convergent-divergent nozzle profiles under different conditions to analyze thrust and specific impulse performances, solving compressible Euler's equations for calorically perfect gas, using MATLAB.

The analysis will be one-dimensional, meaning that all variables will be assumed uniform in each cross-section on the nozzle, depending only on axial position and on time. 

The discretization method will be Finite Differences with uniform space-step and variable time-step, and the scheme will be MacCormack's, which is an explicit 2-step scheme.

Two boundary conditions will be used for the supersonic boundary:

floating variables (extrapolation) which is valid for high chamber pressure and/or low external pressure. Only valid for full supersonic flow at the boundary. No shock capturing.
applied external pressure: more accurate and general (valid for subsonic and supersonic flow) but unstable for high chamber pressure and/or low external pressure. This boundary condition can capture shock formation after the throat, as we will show later.
The adimensianlized nozzle profiles that we will compare are:

