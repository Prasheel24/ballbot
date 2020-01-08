# The Ballbot
## Disclaimer
```
This project is shared for informative purposes. The original records of this project are present on the University server. Please use this project for reference only.
```
## Author
Prasheel Renkuntla - [GitHub](https://github.com/Prasheel24)

I am pursuing my Master's in Robotics at the University of Maryland, College Park. My primary area of interest is in Vision integrated Robot Systems.

## Overview
The ballbot is an omnidirectional balancing mobile robot. This project discusses the author's approach to the design and development of the dynamic model and control architecture of the robot. Balancing Control and Station Keeping Control are two main controllers that can be simulated without a real robot. The results are verified with simulation in MATLAB and SIMULINK.

## Dependencies

* MATLAB - [Installation](https://www.mathworks.com/downloads/)

## Run
To run the main program and design the State Space, run the file - [mainBallBot](https://github.com/Prasheel24/ballbot/blob/master/code/mainBallBot.m)
To verify the simulation results of balancing control, run the file - [balancingControl](https://github.com/Prasheel24/ballbot/blob/master/code/balancingControl.slx)
To run the simulation of the stationkeeping control, use this file - [stationKeepingControl](https://github.com/Prasheel24/ballbot/blob/master/code/stationKeepingControl.slx)

## Demo
The output from the balancing control for the stabilisation of body angles (phi) can be seen in the image below-
<p align="center">
<h5>Balancing Control on Body Angle</h5>
<img src="/output/PhiPlot.PNG" width="70%">
</p>
Following image shows the stabilisation of body angle (phi) when subject to stationkeeping at a given constant point -
<p align="center">
<h5>Station Keeping Control output</h5>
<img src="/output/PhiStationKeeping.PNG" width="70%">
</p>
The trajectory plot is also shown in the output folder.

## References
* T. B. Lauwers, G. A. Kantor, and R. L. Hollis. A dynamically stable single-wheeled
mobile robot with inverse mouse-ball drive. In Proceedings 2006 IEEE International
Conference on Robotics and Automation, 2006. ICRA 2006., pages 2884-2889. IEEE, 2006.
* U. Nagarajan, G. Kantor, and R. Hollis. The ballbot: An omnidirectional balancing
mobile robot. The International Journal of Robotics Research, 33(6):917-930, 2014.
*  U. Nagarajan, G. Kantor, and R. L. Hollis. Trajectory planning and control of an
underactuated dynamically stable single spherical wheeled mobile robot. In 2009 IEEE
International Conference on Robotics and Automation, pages 3743-3748. IEEE, 2009.

