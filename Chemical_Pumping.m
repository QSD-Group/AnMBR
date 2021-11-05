function [Q_CHEM_mgd, NP_CHEM, P_input_CHEM, M_SS_CHEM, M_HDPE_CHEM_kg] = Chemical_Pumping(Q_CHEM_weekly)
% Input:
%  Weekly flow rate of chemical solution, Q_CHEM_weekly [gal/week]


% Chemical container (assume cubic in shape, HDPE in material)
V_CHEM = 2 * Q_CHEM_weekly * 0.00378541; % [m^3] Volume of container holding 2 weeks of chemicals
% 1 gal = 0.00378541 m^3
t_container = 0.003; % [m] Thickness of container
V_HDPE = t_container * (V_CHEM^(1/3))^2 * 6; % [m^3] Volume of container material
Ro_HDPE = 950; % [kg/m^3]
M_HDPE_CHEM_kg = Ro_HDPE * V_HDPE ; % [kg]

% Chemical pumping
Q_CHEM_mgd = Q_CHEM_weekly / 10^6 / 7; % [mgd] Total permeate pumping flow rate
NP_CHEM = 1; % # of permeate pumps in duty

% --> Static head
H_ss_CHEM = V_CHEM^(1/3) * 3.28084; % [ft] Suction Static Head
% 3.28084 ft = 1 m
H_ds_CHEM = 9 + 7/12 - 18/12; % [ft] Suction Static Head
% 9'-7" is the water level in membrane trains
% 18" is the distance from C/L of the pump to the ground
H_ts_CHEM = H_ds_CHEM - H_ss_CHEM; % [ft] Total Static Head

% --> Pressure head
H_p_CHEM = 0; % [ft]

% --> Suction side
% Suction pipe (permeate headers)
L_s_CHEM = 0; % [ft] Length of permeate header per train

% --> Discharge side
% Discharge pipe (Permeate collector)
L_d_CHEM = 30; % [ft] length of permeate collector

[P_input_CHEM, M_SS_CHEM] = Pumping(Q_CHEM_mgd, NP_CHEM, H_ts_CHEM, L_s_CHEM, L_d_CHEM, H_p_CHEM);







end