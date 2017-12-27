# Kidnapped Vehicle | Term 2 P3

In this project, I used a *Particle Filter* to localize a vehicle moving on a simulated map amongst various landmarks. Some key parameters were that I only needed 5 particles or so for decent error values and successfuly localization. Bumping the number of particles above that, to say ~20, decreased the error but not terribly. 

To compile and run the C++ code, we can use the shell scripts provided (in Ubuntu):

1. ./clean.sh
2. ./build.sh
3. ./run.sh

Then, run the simulator and my code will localize the simluated vehicle within the map.