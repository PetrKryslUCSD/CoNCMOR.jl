# CoNCMOR: Coherent Nodal Cluster Model Order Reduction

Library for constructing the "Ritz" transformation matrices for the Coherent Node Cluster model-order reduction.

Based on the paper accepted for publication as NME-Jul-19-0522.R1 - Rapid Free-vibration analysis with Model Reduction based on Coherent Nodal Clusters, in the International Journal for Numerical Methods in Engineering, 03/03/2020. Authors: Krysl, Petr, and Sivapuram, Raghavendra; both of University of California, San Diego , Jacobs School of Engineering Department of Structural Engineering; and
Abawi, Ahmad, of Heat, Light, and Sound Research, Inc.

This package works with [`FinEtools`](https://github.com/PetrKryslUCSD/FinEtools.jl). If you need the model reduction package to work with [`Elfel`](https://github.com/PetrKryslUCSD/Elfel.jl), please use [`CoNCMOR2`](https://github.com/PetrKryslUCSD/CoNCMOR2.jl).

## News

- 02/27/2024: Update for   FinEtools 8. 
- 06/21/2023: Update for   FinEtools 7.0. 
- 07/18/2022: Incompatible changes to allow clusters to have different numbers
  of basis functions (v0.3.0).
