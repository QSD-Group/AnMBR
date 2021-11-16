function [NP_R, P_input_R, M_SS_R] = Retentate_Pumping_CSTR(Q_R_mgd, N_unit, D_train)
% Input:
% Retentate flow rate, Q_R_mgd [mgd]
% Number of cross-flow membrane units, N_unit
% Depth of each anaerobic filter, D_ANA [ft]

% Q_R_mgd was calculated above
NP_R = N_unit; % # of retentate pumps on duty (assume 1 pump per membrane unit)

% --> Static head
H_ss_R = D_train; % [ft] Suction Static Head
H_ds_R = D_train; % [ft] Discharge Static Head
H_ts_R = H_ds_R - H_ss_R; % [ft] Total Static Head 

% --> Pressure head
H_p_R = 0; % [ft] 

% --> Suction side 
% Suction pipe (permeate headers)
L_s_R = 100; % [ft] pipe length/module (!assumption)

% --> Discharge side 
% Discharge pipe (Permeate collector)
L_d_R = 30; % [ft] pipe length/filter (same as discharge side of lift pumping)

[P_input_R, M_SS_R] = Pumping(Q_R_mgd, NP_R, H_ts_R, L_s_R, L_d_R, H_p_R);
end