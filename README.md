# H_MECs_2D
2D Hydro-Mechanical Earthquake Cycle

Computational Earthquake Physics
ETH Zurich, 2022
Dal Zilio, L.,  Hegyi, B.,  Behr, W. M., Gerya, T. (2022)
Hydro-mechanical earthquake cycles in a
poro-visco-elasto-plastic fluid-bearing fault structure
DOI: https://doi.org/10.1016/j.tecto.2022.229516
=====================================================

This project uses the Eigen library.
Install it on Linux/WSL with:
```
sudo apt install libeigen3-dev
```
Install it on Mac OS X with:
```
brew install eigen
```
Additionally this project uses the intel compiler (icpc) to make use the Pardiso solver that uses MKL.
Be aware that you need a valid intel licence to use those.

=====================================================

The c++ version of this simulation works on the EULER cluster only.
To run it on a local device or another cluster, the .sh scripts would need to be changed accordingly.

libraries used:

 - Eigen
 - OpenMP
 - Libszip (similar libraries might work too)
 - HDF5
  
some important flags:

  -O3 (to activate optimizer -> greatly increases speed)

=====================================================

To run the code on EULER use the
```
env2lmod
```
and the .sh script provided in the corresponding folder of this project. (These scripts are written for the EULER cluster only!)

=====================================================

Please use the newest version.
Older versions can be found int the folder H_MEC_RSF/cpp/old_versions. Be careful using these as not all of them work nor guarantee correct results.

=====================================================

Most of the initial conditions (e.g., grid size, number of timesteps and savestep) can be found in the constants.hpp file and have to be changed there.

