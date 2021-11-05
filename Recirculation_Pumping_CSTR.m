function [Q_IR_initial_mgd, NP_IR, P_input_IR, M_SS_IR] = Recirculation_Pumping_CSTR(Q_mgd, IRR, L_train,Q_R_mgd,Upflow_vel_for_GAC,L_membrane_tank,W_membrane_tank,N_train,ia1)
% Input:
% Influent flow rate, Q_mgd [mgd]
% Internal recirculation ratio, IRR
% Length of each train, L_train [ft]

Q_IR_initial_mgd = (Q_mgd * IRR)-Q_R_mgd; % [mgd] Total IR pumping flow rate
if Q_IR_initial_mgd < 0
    Q_IR_initial_mgd = 0;
end
NP_IR = 1; % # of IR pumps in duty

%Additional flow to achieve adequate upflow velocity for GAC
if ia1==2
    Q_upflow_req = ((Upflow_vel_for_GAC*3.28084)*L_membrane_tank*W_membrane_tank*N_train*24)/133681; %[mgd]
    Q_additional_IR = Q_upflow_req-(Q_mgd+Q_IR_initial_mgd);
    if Q_additional_IR <0
        Q_additional_IR = 0;
    end
    Q_IR_mgd = Q_IR_initial_mgd+Q_additional_IR;
else
    Q_IR_mgd = Q_IR_initial_mgd;
end

% --> Static head 
H_ss_IR = 0; % [ft] Suction Static Head
H_ds_IR = 5; % [ft] Discharge Static Head
H_ts_IR = H_ds_IR - H_ss_IR; % [ft] Total Static Head

% --> Pressure head
H_p_IR = 0; % [ft] 

% --> Suction side
% Suction Pipe
L_s_IR = 0; % [ft] (ignore suction side)

% --> Discharge side 
% Discharge Pipe
L_d_IR = L_train; % [ft] pipe length/filter (!assumption)

[P_input_IR, M_SS_IR] = Pumping(Q_IR_mgd, NP_IR, H_ts_IR, L_s_IR, L_d_IR, H_p_IR);
end