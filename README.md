# Ionizer 45
Tools used to optimize the design of the 45-degree extracting surface ionizer.

## Description
- bpm45.cpp
A script to show the beam profile at a given P(w)

- rms45.cpp
A script to search the beam focal point (horizontal/vertical) and draw the beam profile histo

- trajectory_display.cpp
A script to draw the trajectory as a X-Z graph and a X-Y graph

## Run
Make use of the MakeFile provided
(Works with ROOT and C++)

## How it works
- The SIMION coordinates differ from the Inventor coordinates.
- The ion source position and its height from the chamber center (h) can be used to define the "chamber center" (=Inventor origin) with respect to the SIMION origin: \vec{v_{SI}}
- The test plane P(w) is defined by z = -x + \sqrt{2}w, y = (any real number) in Inventor coordinates.
- Each step of the ion \vec{P_f} = (stepfrontX,stepfrontY,stepfrontZ)_{SIMION} = (stepfrontX-src_CPx, stepfrontY-src_CPy, stepfrontZ-src_CPz-h)_{Inventor} and its predecessor \vec{P_b} = (stepbackX,stepbackY,stepbackZ)_{SIMION} = (stepbackX-src_CPx, stepbackY-src_CPy, stepbackZ-src_CPz-h)_{Inventor} are picked out from the CSV data. Step forward until \vec{P_f} passes the plane. This "distance" is defined by (x_0+x_d+z_0+z_d)/sqrt2 - w. Here, x_0, y_0, z_0 are the x,y,z components of \vec{P_b} in Inventor coordinates and x_d,y_d,z_d are the x,y,z components of \vec{P_f} - \vec{P_b} in Inventor coordinates.
- The intercept of the trajectory and P(w) can be given by \vec{t_{trj}}(w) = t_{trj} (\vec{P_f} - \vec{P_b}) + \vec{P_b} with t_{trj} = (\sqrt{w} - x_0 - z_0)/(x_d+z_d).
- Shift to the testplane coordinates by moving \vec{t_trj}(w) to \vec{t_trj}(w) - (w/sqrt2, 0, w/sqrt2).
- +45 deg rotation wrt Y and -90 deg rotation wrt Z will give the X-Y coordinates of the trajectory points on the beam profile monitor (MCP). This yields (trjpty, (trjptz-trjptx)/sqrt2).
