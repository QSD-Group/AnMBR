function [Q_IR_AF_mgd, NP_IR, P_input_IR, M_SS_IR] = Recirculation_Pumping_AF(Q_mgd, R_ANA, N_ANA, D_ANA, Dia_ANA,Q_R_mgd,Upflow_vel_for_GAC,L_membrane_tank,W_membrane_tank,N_train,ia1)
% Input:
    % Influent flow rate, Q_mgd [mgd]
    % Number of anaerobic filters, N_ANA
    % Recirculation ratio, R_ANA
    % Depth of each anaerobic filter, D_ANA [ft]
    % Diameter of each anaerobic filter, Dia_ANA [ft]

Q_IR_AF_initial_mgd = (Q_mgd * R_ANA)-Q_R_mgd; % [mgd] flow rate of each pump
if Q_IR_AF_initial_mgd < 0
    Q_IR_AF_initial_mgd = 0;
end

%Additional flow to achieve adequate upflow velocity for GAC
if ia1==2
    Q_upflow_req = ((Upflow_vel_for_GAC*3.28084)*L_membrane_tank*W_membrane_tank*N_train*24)/133681; %[mgd]
    Q_additional_IR = Q_upflow_req-(Q_mgd+Q_IR_AF_initial_mgd);
    if Q_additional_IR <0
        Q_additional_IR = 0;
    end
    Q_IR_AF_mgd = Q_IR_AF_initial_mgd+Q_additional_IR;
else
    Q_IR_AF_mgd = Q_IR_AF_initial_mgd;
end

NP_IR = N_ANA; % # of IR pumps in duty (assume 1 pump per filter)

% --> Static head
H_ss_IR = D_ANA; % [ft] Suction Static Head
H_ds_IR = D_ANA; % [ft] Discharge Static Head
H_ts_IR = H_ds_IR - H_ss_IR; % [ft] Total Static Head 

% --> Pressure head
H_p_IR = 0; % [ft] 

% --> Suction side 
% Suction pipe 
L_s_IR = D_ANA + Dia_ANA; % [ft] pipe length/filter

% --> Discharge side 
% Discharge pipe
L_d_IR = 30; % [ft] pipe length/filter (same as discharge side of lift pumping)

% % Minor losses, H_m_IR [ft] 
% H_m_IR_1 = 2 * (K1 * VEL_IR ^2 / 2 / g); % 90 degree bend
% H_m_IR_2 = 2 * (K2 * VEL_IR ^2 / 2 / g); % tees  
% H_m_IR = H_m_IR_1 + H_m_IR_2; % [ft] total minor losses

[P_input_IR, M_SS_IR] = Pumping(Q_IR_AF_mgd, NP_IR, H_ts_IR, L_s_IR, L_d_IR, H_p_IR);
end
