function [Q_PERM_mgd, NP_PERM, P_input_PERM, M_SS_PERM] = Permeate_Pumping_Cross_Flow(Q_mgd, N_LU, TMP,D_train,ia1)
% Input:
    % Influent flow rate, Q_mgd [mgd]
    % Number of large membrane units, N_LU
    % Transmembrane pressure, TMP [psi]

% Permeate Pumping 
Q_PERM_mgd = Q_mgd; % [mgd] total permeate flow rate
NP_PERM = N_LU; % # of permeate pumps in duty (assume 1 pump per membrane unit)

% --> Static head
if ia1==3
    H_ss_PERM = 0;
else
    H_ss_PERM = D_train; % [ft] Suction Static Head     
end
H_ds_PERM = D_train; % [ft] Discharge Static Head (based on a unit height of 4.6 m)
H_ts_PERM = H_ds_PERM - H_ss_PERM; % [ft] Total Static Head
                                 
% --> Pressure head
H_p_PERM = TMP * 2.3106; % [ft] TMP in water head

% --> Suction side 
% Suction pipe (permeate headers)
L_s_PERM = 20; % [ft] length of permeate header of each module (based on a 30-module unit length 6 m)

% --> Discharge side 
% Discharge pipe (Permeate collector)
L_d_PERM = 10 * N_LU; % [ft] length of permeate collector (based on a 30-module unit width 1.6 m and space between modules)

[P_input_PERM, M_SS_PERM] = Pumping(Q_PERM_mgd, NP_PERM, H_ts_PERM, L_s_PERM, L_d_PERM, H_p_PERM);
end
