# comet
C++ code for optimising coefficients of individual CVs for PCV MD simulations.

Requires the Eigen library.

The following command is sufficient to compile it on my WSL (Ubuntu):

```
g++ *.cpp -O3 -std=c++14 -lstdc++fs -I/path/to/Eigen -o COMet_2
```
