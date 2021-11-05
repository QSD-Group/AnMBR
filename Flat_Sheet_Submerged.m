function [M_memb_tot, A_LU,V_membrane_displacement] = Flat_Sheet_Submerged(N_train, Cas_per_tank, Mod_per_cas,Module_SA)
% Input: 
    % Flow rate, Q_mgd [mgd]
    % Flux, J [L/m^2-hr]
    % Number of membrane trains, N_train
    % Backwashing time per day, T_bw [hr]

% Output: 
    % Mass of membrane material, M_memb_tot [kg] 
    % Number of small unit, N_SU (In this case, small unit refers to flat-sheet panel)
    % Number of large unit, N_LU_pt 
    % Surface area of each large unit, A_LU [ft^2] 
    
    
% Assume Kubota-RM515 for flat sheet membrane 

% J_fpd = J * 0.07874; % [ft^3/ft^2-day] Membrane flux
                     % unit conversion: 1 L/m^2-hr = 0.07874 ft^3/ft^2-day

% Q_cfd = Q_mgd * 133680.556; % [ft^3/day] 
% Q_bw_cfd = Q_cfd * T_bw / 24; % [ft^3/day] Backwashing flow rate
% A_rqd_tot = (Q_cfd + Q_bw_cfd) / J_fpd; % [ft^2] Total membrane area

A_LU = Mod_per_cas * Module_SA; % [ft^2] Surface area of each large unit

% Membrane Material 
% Kubota-RM515 panel dimensions (L-W-Tickness): 1560-575-6 (mm)
L_SU = 1560 * 10^-3; % [m] Length of small unit
W_SU = 575 * 10^-3; % [m] Width of small unit
t_SU = 4 * 10^-3; % [m] Thickness of small unit (2mm per side of sheet)
t_Permeate_Carrier = 0.002; %[m] 
V_memb_SU = L_SU * W_SU * t_SU; % [m^3] Volume of membrane material per small unit
V_membrane_displacement = (V_memb_SU+(L_SU*W_SU*t_Permeate_Carrier))*Mod_per_cas * Cas_per_tank * N_train;
density_memb = 1.78 * 10^3; % [kg/m^3] Density of membrane material
M_memb_SU = density_memb * V_memb_SU; % [kg] Mass of membrane material per small unit
M_memb_tot = Mod_per_cas * Cas_per_tank * N_train * M_memb_SU; % [kg] Total Mass of membrane material

end

