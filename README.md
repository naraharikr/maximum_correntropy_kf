# Kalman Filtering based on Maximum Correntropy criterion


Ultra-Wideband (UWB) communication is used to localize the robot.  As UWB radio technologies have significant interference within indoor environments, the idea behind this project was to filter out the shot and impulse noises from the location information.  This was part of the localization system, whose filtered data was used to combine with Odometry and IMU data for sensor fusion to precisely get the final location.  The computation of filtering is considered to take place in the UWB tag, which runs on fewer computation capabilities and are powered by battery power.  The aim of using this filtering technique is to obtain the accurate location of the tag (moving service robot within a predetermined map) which is not computational costly.

The implementaion follows the work proposed in the reserach paper.

R. Izanloo, S. A. Fakoorian, H. S. Yazdi and D. Simon, "Kalman filtering based on the maximum correntropy criterion in the presence of non-Gaussian noise," 2016 Annual Conference on Information Science and Systems (CISS), Princeton, NJ, 2016, pp. 500-505, [doi: 10.1109/CISS.2016.7460553] (https://ieeexplore.ieee.org/abstract/document/7460553).


# Instructions

Run mcc_ekf.m


### PS: Copyright 2015 Seyed Abolfazl Fakoorian. All rights reserved. This code can be freely used for noncommercial purposes. This code can be freely distributed provided that this file is included with the distribution. Please provide proper acknowledgment of all uses of this code.
