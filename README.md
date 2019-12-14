# ROBOTICS_2DSCAN_SLAM

### Demo
Click below to play the demo video:
<div align="center">
  <a href="https://www.youtube.com/watch?v=qH-ZQsRhMQU"><img src="https://github.com/xiaofeng419/ROBOTICS_2DSCAN_SLAM/tree/master/Image/SLAM_DEMO.png" alt="IMAGE ALT TEXT"></a>
</div>


### Dataset
The dataset is collected by Dirk Hähnel[1]. The robot platform is equipped with a 180deg FOV 2D lidar and a wheel odometry. There are 910 readings total and each reading contains robot’s x, y and theta(orientation) from odometry and 180 range-bearing reading spanning from -90deg to +90 deg. Both ground truth and raw data are provided with matching timestamp for each data point collected. Figure below shows the plot of both ground truth and raw data respectively. The big red dot is where robot started and each little dot indicates a location where a datapoint is collected. The color of the dots changes from red to purple as robot moving along its trajectory. The robot circled around four laps and ended where the big purple dot is.  

<p align="center">
  <img width="50%" src="https://github.com/xiaofeng419/ROBOTICS_2DSCAN_SLAM/tree/master/Image/GroundTruthLabel.png">
  <img width="48%" src="https://github.com/xiaofeng419/ROBOTICS_2DSCAN_SLAM/tree/master/Image/rawOdometryLabel.png"><br>
  <b>Fig.1 Ground Truth (Left) And Raw Odometry Data (Right) </b><br>
</p>

### Scan Matching 
As we can see from Fig. 1 that the odometry data from raw reading drifted and cannot be directly used. However, laser scan usually gives a much more accurate reading and its accuracy is within 5cm. Therefore, it is more accurate to obtain the robot’s pose via laser reading. 

We applied multi-resolution strategy[2], which is to first downsample the space into 10x more coarse grid and find the best match. That approximately speeds up by 100x. And then we do a more refined search in a much smaller space around the calculated optimal location in coarse scale.  

Fig. 2 illustrates the idea of multi-resolution search. The left image shows the optimal search result under the coarse grid, where yellow cells are scan results of existing map and the red dots are the new laser reading. The coarse search locates the scan result as best as it can under that resolution and then a more refined search as shown on the right side fine tunes its final location to precisely align the two readings. 
 
<p align="center">
  <img  src="https://github.com/xiaofeng419/ROBOTICS_2DSCAN_SLAM/tree/master/Image/ScanMatching.png"><br>
  <b>Fig.2 Scan Matching of Coarse Grid (Left) and Fine Grid (Right)</b><br>
</p>


### FastSLAM[3]:
We approximate the robot’s pose distribution by particles. In other words, at each timestep, the robot’s pose is not a single value but a belief and the belief is updated at each timestep by particle filter. As we know that the SLAM problem can be decomposed as:  
 
  <img src="https://latex.codecogs.com/svg.latex?\Large&space;p(x_{0:t},m\vert{z_{1:t}},u_{1:t})=p(x_{0:t}\vert{z}_{1:t},u_{1:t})p(m\vert{z}_{1:t},x_{1:t})\quad{(1)}"/>
 
where the first term is localization and once we have the trajectory from localization, the second term is a simple mapping problem. Therefore, the key of SLAM is to solve the robot’s trajectory x<sub>0:t</sub>given observation z<sub>1:t</sub>  and moving command u<sub>1:t</sub>. We can apply Bayes’ rule to further break down the first term as:

<img src="https://latex.codecogs.com/svg.latex?\Large&space;p(x_{0:t}\vert{z}_{1:t},u_{1:t})=p(z_{t}\vert{x}_{0:t},z_{1:t}-1,u_{1:t})p(x_{0:t}\vert{z}_{1:t-1},u_{1:t})"/>
 
<img src="https://latex.codecogs.com/svg.latex?\Large&space;=p(z_t\vert{x}_t,z_{1:t-1})p(x_t\vert{x}_{0:t-1},z_{1:t-1},u_{1:t})p(x_{0:t-1}\vert{z}_{1:t-1},u_{1:t-1})"/>

<img src="https://latex.codecogs.com/svg.latex?\Large&space;=p(z_t\vert{x}_t,z_{1:t-1})p(x_{t}\vert{x}_{t-1},u_{1:t})p(x_{0:t-1}\vert{z}_{1:t-1},u_{1:t-1})\quad{(2)"/> 

It is hard to directly sample from (2) and instead we can use importance sampling. Our proposal distribution is:\

<img src="https://latex.codecogs.com/svg.latex?\Large&space;\frac{p(z_t\vert{x}_t,z_{1:t-1})p(x_t\vert{u}_t,x_{t-1})}{\int{p}(z_t│x_t,z_{1:t-1})p(x_t\vert{u}_t,x_{t-1})dx_t}"/> 


which is a combination of both scan matching and motion model and the denominator is just a normalizer. As we use scan matching method in [2], the integration in the denominator is simply the sum of the results of the entire searching area. The weight of importance sample is therefore: 
 
<img src="https://latex.codecogs.com/svg.latex?\Large&space;w_t=\frac{p(z_{t}\vert{x}_t,z_{1:t-1})p(x_t\vert{u}_t,x_{t-1})p(x_{0:t-1}\vert{z}_{1:t-1},u_{1:t-1})}{p(z_t\vert{x}_t,z_{1:t-1})p(x_t\vert{u}_t,x_{t-1})}"/> 

<img src="https://latex.codecogs.com/svg.latex?\Large&space;=p(x_{0:t-1}\vert{z}_{1:t-1},u_{1:t-1})\int{p}(z_t\vert{x}_t,z_{1:t-1})p(x_t\vert{u}_t,x_{t-1})dx_t"/> 

<img src="https://latex.codecogs.com/svg.latex?\Large&space;=w_{t-1}\int{P}(z_t\vert{x}_t,z_{1:t-1})P(x_t\vert{u}_t,x_{t-1})dx_t"/>
 
where w<sub>t-1</sub> is the weight of particles of last time step. Given proposal distribution and the corresponding weight, we can perform particle filtering to update robot’s pose belief given new moving command and Lidar observation. Fig.9 is the result we obtained using 15 particles. We can clearly see that the algorithm successfully produces a globally consistent map and it gives a correct loop closure result when the robot revists the pre-visited area. The algorithm’s speed is 1 update per second. 


### Reference 
[1]. http://ais.informatik.uni-freiburg.de/slamevaluation/datasets.php \
[2]. Edwin B. Olson, Real-time correlative scan matching, 2009 IEEE International Conference on Robotics and Automation
 \
[3]. M. Montemerlo, S. Thrun, D. Koller, and B. Wegbreit, “FastSLAM: A Factored Solution to the Simultaneous Localization and Mapping Problem,” Proc. AAAI Nat’l Conf. Artificial Intelligence, 2002.
