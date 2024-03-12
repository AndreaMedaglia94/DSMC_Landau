%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct Simulation Monte Carlo for the space homogeneous Landau equation %
% Andrea Medaglia (last update: 12/03/2024)                               %
% https://doi.org/10.1016/j.jcp.2024.112845                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

addpath(genpath('functions')); 

% Set the structs containing the parameters
[p_phys, p_sim, p_sch] = set_parameters ;

% Set initial data
[vx, vy, vz] = InitialData(p_sch, p_sim);

% Initialize observables
[obs, distr] = set_observables(vx, vy, vz, p_sim, p_phys);

% Solve the equation
[obs, distr] = Solve(vx, vy, vz, p_sch, p_sim, p_phys, obs, distr) ;

% Post processing 
postprocessingplots(obs,distr,p_sim,p_phys,p_sch) ;


fprintf('Direct Simulation Monte Carlo for the space homogeneous Landau equation \n\n');
fprintf('Code written by Andrea Medaglia (last update 12/03/2024) \n');
fprintf('https://doi.org/10.1016/j.jcp.2024.112845 \n \n');

