Summary
This project is an attempt of simulating the Sun Earth James Webb system for the efforts of extracting lyapunov exponents. 
I wrote this after my 3dSim in an attempt to improve computational efficiency and the way I implemented the system. 
Unfortunately, this approach does not work fundamentally as the James Webb Telescope in the L2 lagrange point orbits in a plane perpendicular to the two other coordinates.
Thus, the system is fundamentally 3 dimensions. 

Directions for using the code:

The computational part of this code is all written in C++, and is included in the main script and src directory. All headers and implementations are included.
However, #include paths will need to be customized for your own computer.
Additionally, the data file output paths in the "dataToFile" function of "SingleProbeChaos.cpp" will need to be customized for your own hardware.

The main script allows customization of time-step, total steps taken (timestep * total steps = total time), and how many time steps you want to skip before you store
data. The Sun, Earth, and James Webb are already set up in their initial positions, a total number of perturbations (offset james webb satellites), and distance of offset
is chosen. The main function "SingleProbeChaos" returns a pointer to an array of two doubles, which are supposed to be the lyapunov exponents after the simulation has
completed the number of timesteps specified. 

Because the lyapunovs cannot be calculated accurately with this method, the real power in this program is it's ability to plot trajectories over time of different
planets and the Sun. Data is stored in a 2d vector2 array every timestep, and then written to file at the end of the simulation. 

For data analysis purposes, I have included a number of Python files I used to plot and animate different output data. 
The file paths will need to be customized for your own hardware, and for the creation of mp4 files, I used a free software called FfmPeg. 
If you desire to create animated plots, you will need to download this software and replace the path to the executable under "Desktop ffmpeg path".

 While the project doesn't correctly simulate the Sun Earth James Webb system, it does correctly simulate the Sun Earth system and I believe it improved my skills
 in designing what was originally a messy and unorganized approach. 
