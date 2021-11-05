function [M_memb_tot, A_LU,V_membrane_displacement] = Hollow_Fiber(N_train, Cas_per_tank, Mod_per_cas,Module_SA)
% Input: 
    % Flow rate, Q_mgd [mgd]
    % Flux, J [L/m^2-hr]
    % Number of membrane trains, N_train
    % Backwashing time per day, T_bw [hr]

% Output: 
    % Mass of membrane material, M_memb_tot [kg] 
    % Number of small unit, N_SU (In this case, small unit refers to Zenon-ZeeWeed*500D module)
    % Number of large unit, N_LU_pt 
    % Surface area of each large unit, A_LU [ft^2] 
    
    
% Assume Zenon-ZeeWeed*500D Cassette for hollow fiber membrane
% (Membrane module specs available at:http://www.gewater.com/products/zeeweed-500-membrane.html)

% J_fpd = J * 0.07874; % [ft^3/ft^2-day] Membrane flux
                     % unit conversion: 1 L/m^2-hr = 0.07874 ft^3/ft^2-day
                     % (Checked by Brian - it's good)

% Q_cfd = Q_mgd * 133680.556; % [ft^3/day] 
% Q_bw_cfd = Q_cfd * T_bw / 24; % [ft^3/day] Backwashing flow rate
% A_rqd_tot = (Q_cfd + Q_bw_cfd) / J_fpd; % [ft^2] Total membrane area

A_LU = Mod_per_cas * Module_SA; % [ft^2] Surface area of each large unit

% Membrane Material 
OD_fiber = 1.9 * 10^-3; % [m] Outer diameter of each membrane fiber
ID_fiber = 0.8 * 10^-3; % [m] Inner diameter of each membrane fiber
L_fiber = 2.198; % [m] Length of each fiber 
A_fiber = L_fiber * pi * OD_fiber; % [m^2] Surface area of each fiber
V_fiber = L_fiber * pi/4 * (OD_fiber^2 - ID_fiber^2); % [m^3] Volume of fiber material
V_membrane_displacement = L_fiber * pi/4 * (OD_fiber^2); % Reactor volume lost due to membranes 
Number_fibers_per_module = Module_SA/A_fiber;

V_memb_SU = Number_fibers_per_module * V_fiber; % [m^3] Volume of membrane material of each small unit 
density_memb = 1.78 * 10^3; % [kg/m^3] Density of membrane material
M_memb_SU = density_memb * V_memb_SU; % [kg] Mass of membrane material per small unit
M_memb_tot = Mod_per_cas * Cas_per_tank * N_train * M_memb_SU; % [kg] Total Mass of membrane material

end