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

% % % post processing % % %
% M4_numerics = cell2array(obs,5) ;    
% [~, ~, M4] = BKW_analytic(p_sim, p_sim.time_obs) ;

DT = Trubnikov_analytic(p_sim, p_phys, p_sch, p_sim.time_obs);
Txyz_numerics = cell2array(obs,3) ;

plot(p_sim.time_obs, DT, 'k-')
hold on
plot(p_sim.time_obs, (Txyz_numerics(1,:)-Txyz_numerics(3,:))./(Txyz_numerics(1,1)-Txyz_numerics(3,1)), 'ro')
hold off




