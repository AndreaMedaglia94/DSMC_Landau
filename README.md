# **DSMC for the space homogeneous Landau equation**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This MATLAB code solves the space-homogeneous Landau equation in the deterministic scenario with a Direct Simulation Monte Carlo (DSMC) approach. 
It is based on the paper [A. Medaglia, L. Pareschi, M. Zanella, _Particle simulation methods for the Landau-Fokker-Planck equation with uncertain data_,  J. Comput. Phys., **503** (2024), 112845](https://www.sciencedirect.com/science/article/pii/S0021999124000949?ref=pdf_download&fr=RR-2&rr=867703bb2d980e46).

The novelty of the code is that it avoids the use of iterative solvers to sample the approximate surrogate collisional kernel, as in the original paper [A. V. Bobylev, K. Nanbu, _Theory of collision algorithms for gases and plasmas based on the Boltzmann equation and the Landau-Fokker-Planck equation_, Phys. Rev. E, 61(4), 4576](https://journals.aps.org/pre/pdf/10.1103/PhysRevE.61.4576).
