function [M_memb_tot, Q_R_mgd] = Multi_Tube(Q_mgd, J, T_bw, VEL_xflow, N_train,Cas_per_tank, Mod_per_cas,Module_SA)
% Input: 
    % Flow rate, Q_mgd [mgd]
    % Flux, J [L/m^2-hr]
    % Backwashing time per day, T_bw [hr]
    % Cross-flow velocity, VEL_xflow [m/s]
    
% Output: 
    % Mass of membrane material, M_memb_tot [kg] 
    % Number of small unit, N_SU (In this case, small unit refers to Pentair X-flow model AQFMBR 30 module)
    % Number of large unit, N_LU_pt 
    % Retentate flow rate, Q_R_mgd [mgd] 
    
    
% Assume Pentair X-flow model AQFMBR 30
% (Membrane information available at: % Info from http://onlinembr.info/Membrane%20process/Airlift.htm
% http://biomath.ugent.be/biomath/publications/download/ratkovichnicolas_phd.pdf)

% J_m3pm2d = (J / 10^3) * 24; % [m^3/m^2-day] Membrane flux
                            % unit conversion: [L/m^2-hr] to [m^3/m^2-day]

% Q_cmd = Q_mgd * 3785.41178; % [m^3/day]
% Q_bw_cmd = Q_cmd * T_bw / 24; % [m^3/day] Backwashing flow rate
% A_rqd_tot = (Q_cmd + Q_bw_cmd) / J_m3pm2d; % [ft^2] Total membrane area

A_LU = Mod_per_cas * Module_SA;

% N_SU_pLU = 30; % Number of small units per large unit (assume 30 Pentair X-flow model AQFMBR 30 modules)
% A_LU = N_SU_pLU * A_module; % [ft^2] % Surface area of each large unit
% N_LU = A_rqd_tot / A_LU; % Total number of large units
% N_SU = N_LU * N_SU_pLU; % % Total number of small units
               
Q_xflow = 53.5 * VEL_xflow; % [m^3/hr] cross-flow flow rate per module, based on manuf. specs. for compact 33
Q_R_cmh = Mod_per_cas*Cas_per_tank*N_train * Q_xflow; % [m^3/hr] Total retentate flow rate 
Q_R_mgd = Q_R_cmh * 0.00634; % [mgd] unit conversion

% Membrane Material (correct)
OD_tube = 6 * 10^-3; % [m] Outer diameter of each membrane tube
ID_tube = 5.2 * 10^-3; % [m] Inner diameter of each membrane tube
L_tube = 3; % [m] Length of each membrane tube
N_tube = 700; % Number of tubes in each small unit
V_memb_tube = L_tube * pi/4 * (OD_tube^2 - ID_tube^2); % [m^3] Volume of each membrane tube
V_memb_SU = N_tube * V_memb_tube; % [m^3] Volume of membrane material per small unit
density_memb = 1.78 * 10^3; % [kg/m^3] Density of membrane material
M_memb_SU = density_memb * V_memb_SU; % [kg] [kg] Mass of membrane material per small unit
M_memb_tot = Mod_per_cas*Cas_per_tank*N_train * M_memb_SU; % [kg] Total Mass of membrane material

end
