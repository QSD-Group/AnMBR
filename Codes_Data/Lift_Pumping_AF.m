function [Q_LIFT_mgd, NP_LIFT, P_input_LIFT, M_SS_LIFT] = Lift_Pumping_AF(Q_mgd, N_ANA, D_ANA)
% Input:
% Influent flow rate, Q_mgd [mgd]
% Number of anaerobic filters, N_ANA
% Depth of each anaerobic filter, D_ANA [ft]

Q_LIFT_mgd = Q_mgd; % [mgd] total lift pumping flow rate
NP_LIFT = N_ANA; % # of lift pumps in duty (assume 1 pump per filter)

% --> Static head
H_ss_LIFT = 0; % [ft] Suction Static Head
H_ds_LIFT = D_ANA; % [ft] Discharge Static Head
H_ts_LIFT = H_ds_LIFT - H_ss_LIFT; % [ft] Total Static Head

% --> Pressure head
H_p_LIFT = 0; % [ft] 

% --> Suction side 
% Suction pipe 
L_s_LIFT = 150; % [ft]Length of suction pipe per filter

% --> Discharge side 
% Discharge pipe
L_d_LIFT = 30; % [ft] pipe length/filter (!assumption)
                                   
% % Minor losses, H_m_LIFT [ft] 
% K1 = 0.30; % 90 degree bend
% K2 = 0.30; % tees    
% K3 = 1.0; % exit loss
% K4 = 0.5; % entrance loss
% g = 32.17405; % [ft/s2] Acceleration of gravity 
% 
% H_m_LIFT_3 = 1 * (K3 * VEL_LIFT ^2 / 2 / g); % exit loss
% H_m_LIFT = H_m_LIFT_3; % [ft] total minor losses

[P_input_LIFT, M_SS_LIFT] = Pumping(Q_LIFT_mgd, NP_LIFT, H_ts_LIFT, L_s_LIFT, L_d_LIFT, H_p_LIFT);
end