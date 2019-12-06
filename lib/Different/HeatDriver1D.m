% Driver script for solving the 1D advection equations
% with variable coefficient
Globals1D;
% Polynomial order used for approximation
N = 4;
% Read in Mesh
[Nv, VX, K, EToV] = MeshGen1D(0,2*pi,100);
% Initialize solver and construct grid and metric
StartUp1D;
% Set initial conditions
u = sin(x);
% Solve Problem
FinalTime = 0.1;
[u,time] = Heat1D(u,FinalTime);