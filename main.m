close all
clear all
clc

addpath(genpath('functions')); 

% Set the structs containing the parameters
[p_phys, p_sim, p_sch] = set_parameters ;

% Set initial data
[vx, vy, vz] = InitialData(p_sch, p_sim, p_phys);

% Initialize observables
[obs, distr] = set_observables(vx, vy, vz, p_sim, p_phys);

% Solve the equation
[obs, distr] = Solve(vx, vy, vz, p_sch, p_sim, p_phys, obs, distr) ;

% plot results
M4_numerics = cell2array(obs,5) ;    
[~, ~, M4] = BKW_analytic(p_sim, p_sim.time_obs) ;

plot(p_sim.time_obs, M4, 'k-')
hold on
plot(p_sim.time_obs, M4_numerics, 'ro')
hold off




