function [M_memb_tot, Q_R_mgd] = Flat_Sheet_Xflow(Q_mgd, J, T_bw, VEL_xflow, N_train, Cas_per_tank, Mod_per_cas,Module_SA)
% Input: 
    % Flow rate, Q_mgd [mgd]
    % Flux, J [L/m^2-hr]
    % Backwashing time per day, T_bw [hr]
    % Cross-flow velocity, VEL_xflow [m/s]
    
% Output: 
    % Mass of membrane material, M_memb_tot [kg] 
    % Number of small unit, N_SU (In this case, small unit refers to flat-sheet panel)
    % Number of large unit, N_LU_pt 
    % Retentate flow rate, Q_R_mgd [mgd] 
   
    
% Assume Kubota-RM515 for flat sheet membrane 

J_fpd = J * 0.07874; % [ft^3/ft^2-day] Membrane flux
                     % unit conversion: 1 L/m^2-hr = 0.07874 ft^3/ft^2-day
% 
% Q_cfd = Q_mgd * 133680.556; % [ft^3/day] 
% Q_bw_cfd = Q_cfd * T_bw / 24; % [ft^3/day] Backwashing flow rate
% A_rqd_tot = (Q_cfd + Q_bw_cfd) / J_fpd; % [ft^2] Total membrane area

A_LU = Mod_per_cas * Module_SA;

% A_SU = 1.45 * 10.7639; % [ft^2] Membrane surface area of each small unit 
% N_SU_pLU = 150; % Number of small units per large unit (assume 150-200 panels, Kubota-RM515)
% A_LU = N_SU_pLU * A_SU; % [ft^2] Surface area of each large unit
% N_LU = ceil(A_rqd_tot / A_LU); % Total number of large units
% N_SU = N_LU * N_SU_pLU; % Total number of small units

Q_xflow = 53.5 * VEL_xflow; % [m^3/hr] cross-flow flow rate per module, based on manuf. specs. for compact 33 Mult-tube membrane (Assume it's the same efficiency for flat sheet)
Q_R_cmh = Mod_per_cas*Cas_per_tank*N_train * Q_xflow; % [m^3/hr] Total retentate flow rate 
Q_R_mgd = Q_R_cmh * 0.00634; % [mgd] unit conversion

% Membrane Material 
% Kubota-RM515 panel dimensions (L-W-Tickness): 1560-575-6 (mm)
L_SU = 1560 * 10^-3; % [m] Length of small unit
W_SU = 575 * 10^-3; % [m] Width of small unit
t_SU = 4 * 10^-3; % [m] Thickness of small unit
V_memb_SU = L_SU * W_SU * t_SU; % [m^3] Volume of membrane material per small unit
density_memb = 1.78 * 10^3; % [kg/m^3] Density of membrane material
M_memb_SU = density_memb * V_memb_SU; % [kg] Mass of membrane material per small unit
M_memb_tot = Mod_per_cas*Cas_per_tank*N_train * M_memb_SU; % [kg] Total Mass of membrane material

end
