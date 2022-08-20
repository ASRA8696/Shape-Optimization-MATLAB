# CD-Nozzle
Shape optimization study of a rocket nozzle using MATLAB


Nozzles are aerospace components found in propulsion devices such as turbojet engines or rocket engines. Its purpose is to accelerate high-pressure, high-temperature gases as much as possible to get optimal thrust for vehicle propulsion, according to Newton's 3rd law.

In this project, we are going to simulate 3 different convergent-divergent nozzle profiles under different conditions to analyze thrust and specific impulse performances, solving compressible Euler's equations for calorically perfect gas, using MATLAB.

The analysis will be one-dimensional, meaning that all variables will be assumed uniform in each cross-section on the nozzle, depending only on axial position and on time. 

The discretization method will be Finite Differences with uniform space-step and variable time-step, and the scheme will be MacCormack's, which is an explicit 2-step scheme.


The non-dimensianlized nozzle profiles that we will compare are:


   ![shape-function](https://user-images.githubusercontent.com/79316741/185670660-4f5cfa6d-72b7-44cf-ac4f-07a98414e1c7.jpg)


STEADY STATE PROFILES FOR FULLY EXPANDED FLOWS WITHOUT SHOCK

https://user-images.githubusercontent.com/79316741/185671430-132acfc0-2eb7-4e4f-8402-ce1c7582be07.mp4


STEADY STATE PROFILES FOR FLOWS WITH A NORMAL SHOCK

https://user-images.githubusercontent.com/79316741/185671587-3e933def-cb2b-42c2-9d64-8697a22cac29.mp4

2D CONTOUR PROFILES OF NORMAL SHOCK

https://user-images.githubusercontent.com/79316741/185672420-423e9793-fec5-4823-a87d-82bb4d77724d.mp4

THRUST & SPECIFIC IMPULSE CALCULATION

There are two key indicators for performance of an engine nozzle: Thrust force and Specific impulse.

Thrust is the force that will accelerate the rocket. It strongly depends on the scale of the nozzle (the bigger, the more thrust)

Specific impulse is a measure of the thrust performance of a nozzle, regardless of its scale, so it only depends on the nozzle shape and operating conditions (reservoir pressure and temperature, and external pressure).


![thrust](https://user-images.githubusercontent.com/79316741/185731575-6b8e2f81-f045-4353-9263-e122a8df4396.jpg)


![impulse](https://user-images.githubusercontent.com/79316741/185731580-46eeea92-7598-4f0c-a30a-6120b655f65f.jpg)




