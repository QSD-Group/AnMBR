function [Q_IR_mgd, NP_IR, P_input_IR, M_SS_IR] = Recirculation_Pumping_Cross_Flow(Q_mgd, R_ANA, N_ANA, D_ANA, Dia_ANA)
% Input:
    % Influent flow rate, Q_mgd [mgd]
    % Number of anaerobic filters, N_ANA
    % Recirculation ratio, R_ANA
    % Depth of each anaerobic filter, D_ANA [ft]
    % Diameter of each anaerobic filter, Dia_ANA [ft]

Q_IR_mgd = Q_mgd * R_ANA; % [mgd] flow rate of each pump
NP_IR = N_ANA; % # of IR pumps in duty (assume 1 pump per filter)

% --> Static head
H_ss_IR = D_ANA; % [ft] Suction Static Head
H_ds_IR = D_ANA; % [ft] Discharge Static Head
H_ts_IR = H_ds_IR - H_ss_IR; % [ft] Total Static Head 

% --> Pressure head
H_p_IR = 0; % [ft] 

% --> Suction side 
% Suction pipe 
L_s_IR = D_ANA + Dia_ANA; % [ft] pipe length/filter (!assumption)

% --> Discharge side 
% Discharge pipe
L_d_IR = 30; % [ft] pipe length/filter (same as discharge side of lift pumping)

% % Minor losses, H_m_IR [ft] 
% H_m_IR_1 = 2 * (K1 * VEL_IR ^2 / 2 / g); % 90 degree bend
% H_m_IR_2 = 2 * (K2 * VEL_IR ^2 / 2 / g); % tees  
% H_m_IR = H_m_IR_1 + H_m_IR_2; % [ft] total minor losses

[P_input_IR, M_SS_IR] = Pumping(Q_IR_mgd, NP_IR, H_ts_IR, L_s_IR, L_d_IR, H_p_IR);
end
