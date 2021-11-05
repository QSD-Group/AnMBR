function [Q_PERM_mgd, NP_PERM, P_input_PERM, M_SS_PERM] = Permeate_Pumping_Submerged(Q_mgd, N_train, D_train, L_train, W_N_trains, TMP,ia1)
% Input:
    % Influent flow rate, Q_mgd [mgd]
    % Number of trains, N_train
    % Depth of each train, D_train [ft]
    % Length of each train, L_train [ft]
    % Width of N trains, W_N_trains [ft]
    % Transmembrane pressure, TMP [psi]

Q_PERM_mgd = Q_mgd; % [mgd] Total permeate pumping flow rate
NP_PERM = N_train; % # of permeate pumps in duty (assume 1 pump per train)

% --> Static head
if ia1==3
    H_ss_PERM = 0;
else
    H_ss_PERM = D_train - 18/12; % [ft] Suction Static Head     
               % 9'-7" is the water level in membrane trains
               % 18" is the distance from C/L of the pump to the ground
end
H_ds_PERM = D_train - 18/12; % [ft] 
H_ts_PERM = H_ds_PERM - H_ss_PERM; % [ft] Total Static Head
                                 
% --> Pressure head
H_p_PERM = TMP * 2.3106; % [ft] TMP in water head

% --> Suction side 
% Suction pipe (permeate headers)
L_s_PERM = L_train + 4.5/17*(75 - 22); % [ft] Length of permeate header per train

% --> Discharge side 
% Discharge pipe (Permeate collector)
L_d_PERM = W_N_trains; % [ft] length of permeate collector 

[P_input_PERM, M_SS_PERM] = Pumping(Q_PERM_mgd, NP_PERM, H_ts_PERM, L_s_PERM, L_d_PERM, H_p_PERM);
end