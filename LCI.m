function [INV_CON, INV_OP, DEmis, Power_pct, E_input_kWh, E_offset_kWh, V_treated, Output_LCC] = LCI(ia,ia1,ib,ic,id,ih, DVariable, UParameter,interest)
% Input:
% step_A (Reactor type) input 'CSTR' or 'AF'
% step_A1: If there''s aerobic polishing filters added, input '+AER'
%          If there''s GAC added, input '+GAC'
% step_B (Configuration): input 'Submerged' or 'Cross-flow'
% step_C (Membrane type) input 'MT' for "multi-tube", or 'HF' for "hollow fiber", or 'FS' for "flat sheet"
% step_D (Membrane material) input 'PET' or 'PTFE'
% step_F (Physical cleaning)
% step_H (Soluble methane management)
% step_I (Methane processing)
% DVariable: A vector containing Decision variables
% UParameter: A vector containing Uncertainty parameters

% Output:
% Construction inventory matrix, INV_CON
% Operational inventory matrix, INV_OP
% Direct emission matrix, DEmis
% Energy source percentage matrix, Power_pct
% Total energy input over N years, E_input_kWh [kWh]
% Total energy offset over N years, E_offset_kWh [kWh]
% Volume of treated wastewater overy N years, V_treated [m^3]
% Output matrix for LCC, Output_LCC


Year = 30; % Assumed years of operation

%% Influent
Q_mgd = UParameter(1); % Influent volumetric flow rate [mgd]
V_treated = Q_mgd * 3785.41178 * Year * 365; % [m^3] Volume of treated wastewater overy N years
S_SO = UParameter(2); % [mg-COD/L] or [g/m^3]
X_SO = UParameter(3); % [mg-COD/L] or [g/m^3]
VEL_xflow = UParameter(4); % [m/s] Cross-flow velocity
Upflow_vel_for_GAC = UParameter(9);
COD_inf = S_SO + X_SO; % [mg-COD/L] or [g/m^3]

%% Effluent
COD_eff = 30; % [mg-COD/L] or [g/m^3]
NH3_eff = 0; % [mg-N/L]
NH4_eff = 25; % [mg-N/L] % Hospido,2012
Org_N_eff = 0; % [mg-N/L]
P_eff = 0; % [mg-P/L]
PO4_eff = 3; % [mg-P/L] % Hospido,2012
CH4_eff = 0; % [mg-CH4/L] % Hospido,2012

% Assign values to decision variables
if ia==1 && ia1==1 %CSTR_none
    J = DVariable(1);
    TMP = DVariable(2);
    IRR = DVariable(3);
    HRT = DVariable(4);
    HRT_for_GAC = DVariable(5);
elseif ia==1 && ia1==2 %CSTR_GAC
    J = DVariable(1);
    TMP = DVariable(2);
    IRR = DVariable(3);
    HRT = DVariable(4);
    HRT_for_GAC = DVariable(5);
elseif ia==2 && ia1==1 %AF_none
    J = DVariable(1);
    TMP = DVariable(2);
    OL_AF = DVariable(3);
    HL_AF = DVariable(4);
    R_AF = DVariable(5);
elseif ia==2 && ia1==2
    J = DVariable(1);
    TMP = DVariable(2);
    OL_AF = DVariable(3);
    HL_AF = DVariable(4);
    R_AF = DVariable(5);
    HRT_for_GAC = DVariable(6);
elseif ia==2 && ia1==3 %AF_AER
    J = DVariable(1);
    TMP = DVariable(2)/2; %Assume AeF allows for lower TMP
    OL_AF = DVariable(3);
    HL_AF = DVariable(4);
    R_AF = DVariable(5);
    OL_AER = DVariable(6);
    HL_AER = DVariable(7);
end

%% H&S Membrane tank design criteria
N_train=2;
if ic==1 %Hollow fiber
    Module_SA = 370; %Membrane surface area per Zenon Membrane Module [ft2]
    Mod_per_cas = 30; %Zenon 500D Modules per Cassette (Max = 48)
    Contact_SA = Module_SA*Mod_per_cas; %Contact Surface Area per Zenon Membrane Cassette [ft2]
    Cas_per_tank = 16; %Cassettes per membrane tank
    Spare_cassettes = 2; %Spare cassettes per membrane tank
    %     Flux_all = (Q_mgd*1000000)/(N_train*Cas_per_tank*Contact_SA);
    Flux_oneoff = (Q_mgd*1000000)/((N_train-1)*Cas_per_tank*Contact_SA);
    while Flux_oneoff>J*0.5754 %J converted to gal/ft^2-d by multiplying by 24 hr/d and dividing by (10.7639 ft2/m2 and 3.785 L/gal)
        Mod_per_cas=Mod_per_cas+1;
        if Mod_per_cas==49
            if Cas_per_tank==23
                N_train=N_train+1;
                Cas_per_tank=16;
                Mod_per_cas = 30;
            else
                Cas_per_tank = Cas_per_tank+1;
                Mod_per_cas = 30;
            end
        end
        Contact_SA = Module_SA*Mod_per_cas;
        Flux_oneoff = (Q_mgd*1000000)/((N_train-1)*Cas_per_tank*Contact_SA);
    end
elseif ic==2 %Flat sheet
    Module_SA = 1.45 * 10.7639; %Membrane surface area per sheet (assume 150-200 panels, Kubota-RM515) [ft2]
    Mod_per_cas = 150; %(Max = 200)
    Contact_SA = Module_SA*Mod_per_cas; %Contact Surface Area per Kubota sheet [ft2]
    Cas_per_tank = 16; %Cassettes per membrane tank
    Spare_cassettes = 2; %Spare cassettes per membrane tank
    %         Flux_all = (Q_mgd*1000000)/(N_train*Cas_per_tank*Contact_SA);
    Flux_oneoff = (Q_mgd*1000000)/((N_train-1)*Cas_per_tank*Contact_SA);
    while Flux_oneoff>J*0.5754 %J converted to gal/ft^2-d by multiplying by 24 hr/d and dividing by (10.7639 ft2/m2 and 3.785 L/gal)
        Mod_per_cas=Mod_per_cas+1;
        if Mod_per_cas==201
            if Cas_per_tank==23
                N_train=N_train+1;
                Cas_per_tank=16;
                Mod_per_cas = 150;
            else
                Cas_per_tank = Cas_per_tank+1;
                Mod_per_cas = 150;
            end
        end
        Contact_SA = Module_SA*Mod_per_cas;
        Flux_oneoff = (Q_mgd*1000000)/((N_train-1)*Cas_per_tank*Contact_SA);
    end
elseif ic==3 %MT
    Module_SA = 32*10.7639; %Membrane surface area per tube; assume Pentair X-flow model AQFMBR 30 [ft2]
    Mod_per_cas = 44; %Max = 48
    Contact_SA = Module_SA*Mod_per_cas; %Contact Surface Area per Pentair X-flow [ft2]
    Cas_per_tank = 16; %Cassettes per membrane tank
    Spare_cassettes = 2; %Spare cassettes per membrane tank
    %     Flux_all = (Q_mgd*1000000)/(N_train*Cas_per_tank*Contact_SA);
    Flux_oneoff = (Q_mgd*1000000)/((N_train-1)*Cas_per_tank*Contact_SA);
    while Flux_oneoff>J*0.5754 %J converted to gal/ft^2-d by multiplying by 24 hr/d and dividing by (10.7639 ft2/m2 and 3.785 L/gal)
        Mod_per_cas=Mod_per_cas+1;
        if Mod_per_cas==49
            if Cas_per_tank==23
                N_train=N_train+1;
                Cas_per_tank=16;
                Mod_per_cas = 44;
            else
                Cas_per_tank = Cas_per_tank+1;
                Mod_per_cas = 44;
            end
        end
        Contact_SA = Module_SA*Mod_per_cas;
        Flux_oneoff = (Q_mgd*1000000)/((N_train-1)*Cas_per_tank*Contact_SA);
    end
end

N_cassettes = Cas_per_tank*N_train;
L_membrane_tank = ceil(Cas_per_tank*3.4+Spare_cassettes*3.4); %[ft]
W_membrane_tank = 21; %[ft]
D_membrane_tank = 12; %[ft]

%Design_criteria = input('Enter design criterion [gal/ft2/d]: ');
% Tank_footprint = N_train*L_membrane_unit*W_membrane_unit;

%Membrane Sparging (H&S)
SGD = UParameter(5); % Specific gas demand [Nm^3 gas/m^2 membrane area-h] (Robles, 2012=0.23)
SGD_cfm_per_module = (SGD*(Module_SA/(3.28084^2))*(3.28084^3))/60; % Converting SGD from m3/m2-hr to cfm/module
Gas_demand = SGD_cfm_per_module; %[cfm/module]
% Air_demand_10_30 = 3; %[cfm/module]
Air_demand_instantaneous = ceil(N_train*Gas_demand*(Mod_per_cas*Cas_per_tank)); %[cfm/module]
% Air_demand_10_30_instantaneous = ceil(N_train*Air_demand_10_30*(Mod_per_cas*Cas_per_tank)); %[cfm/module]

if ia1~=2 && ib==1
    TCFM=Air_demand_instantaneous;
    if (TCFM<=30000)
        N=1;
        CFMB = TCFM/N;
        while (CFMB>7500)
            N=N+1;
            CFMB = TCFM/N;
        end
        NB = N +1;
    elseif (TCFM>30000)&&(TCFM<=72000)
        N=1;
        CFMB = TCFM/N;
        while (CFMB>18000)
            N=N+1;
            CFMB = TCFM/N;
        end
        NB = N +1;
    elseif (TCFM>72000)
        N=1;
        CFMB = TCFM/N;
        while (CFMB>100000)
            N=N+1;
            CFMB = TCFM/N;
        end
        NB = N +1;
    end
end



%% Step A - Reactor Type
if ia==1
    if ib==1
        % Reactor size depends on # of membrane modules submerged, so membrane functions are called here.
        if ic==1
            [M_Membrane_kg, A_LU,V_membrane_displacement] = Hollow_Fiber(N_train, Cas_per_tank, Mod_per_cas,Module_SA); % Call membrane function (hollow fiber)
        elseif ic==2
            [M_Membrane_kg, A_LU,V_membrane_displacement] = Flat_Sheet_Submerged(N_train, Cas_per_tank, Mod_per_cas,Module_SA); % Call membrane function (flat sheet)
        end
        V_membrane_tank = N_train * (D_membrane_tank * L_membrane_tank * W_membrane_tank-V_membrane_displacement); % [ft^3] Volume of reactor
        HRT_membrane_tank = V_membrane_tank/(Q_mgd * 133681 / 24); % [hr];
        [D_train, L_CSTR, W_train, W_N_trains, W_PB, L_BB, VWC_CSTR, VSC_CSTR, VEX_CSTR] = CSTR_Submerged(Q_mgd, HRT, N_train, Cas_per_tank,HRT_membrane_tank,L_membrane_tank,ia1,HRT_for_GAC); % Call CSTR function (submerged membrane)

        [~, ~, P_input_PERM, M_SS_PERM] = Permeate_Pumping_Submerged(Q_mgd, N_train, D_train, L_membrane_tank, W_N_trains, TMP,ia1); % Call permeate pumping function, ~ are Q_PERM_mgd, NP_PERM

        if ia1==1
            freq = UParameter(6); % freq: Gas sparging frequency [sparging time/total time]
            [~, P_input_blower, M_SS_gas, ~, ~] = Gas_Sparging_Submerged(N_train, Cas_per_tank, A_LU, SGD, L_membrane_tank, W_PB, L_BB, freq); % Call gas sparging function, ~ are Q_gas_cfm, OD_gh, OD_gsm
        end
        Q_R_mgd = 0;

    elseif ib==2
        [D_train, L_CSTR, W_train, ~, ~, ~, VWC_CSTR, VSC_CSTR, VEX_CSTR] = CSTR_Cross_Flow(Q_mgd, HRT, N_train,L_membrane_tank); % Call CSTR function (with cross-flow membrane), ~ are W_N_trains, W_PB, L_BB
        if ic==3
            [M_Membrane_kg, Q_R_mgd] = Multi_Tube(Q_mgd, J, 0.4, VEL_xflow,N_train,Cas_per_tank, Mod_per_cas,Module_SA);
            [~, ~, P_input_PERM, M_SS_PERM] = Permeate_Pumping_Cross_Flow(Q_mgd, Cas_per_tank, TMP,D_train,ia1); % ~ are Q_PERM_mgd, NP_PERM
        elseif ic==2
            [M_Membrane_kg, Q_R_mgd] = Flat_Sheet_Xflow(Q_mgd, J, 0.4, VEL_xflow,N_train,Cas_per_tank, Mod_per_cas,Module_SA);
            [~, ~, P_input_PERM, M_SS_PERM] = Permeate_Pumping_Cross_Flow(Q_mgd, Cas_per_tank, TMP,D_train,ia1); % ~ are Q_PERM_mgd, NP_PERM
        end
        [~, P_input_R, M_SS_R] = Retentate_Pumping_CSTR(Q_R_mgd, Cas_per_tank, D_train); % ~ is NP_R
    end

    V_rctr = N_train * D_train * L_CSTR * W_train; % [ft^3] Volume of reactor
    [~, ~, P_input_IR, M_SS_IR] = Recirculation_Pumping_CSTR(Q_mgd, IRR, L_membrane_tank,Q_R_mgd,Upflow_vel_for_GAC,L_membrane_tank,W_membrane_tank,N_train,ia1); % ~ are Q_IR_mgd,NP_IR

elseif ia==2
    if ib==1
        [V_m_AF, D_AF, Dia_AF, N_AF, VWC_AF, VSC_AF, VEX_AF, W_PB, L_BB,W_N_trains] = AF_Submerged(Q_mgd, S_SO, X_SO, OL_AF, HL_AF, R_AF,Cas_per_tank,N_train,L_membrane_tank); % Call anaerobic filter function
        V_rctr = N_AF * V_m_AF / 0.5; % [ft^3] Volume of reactor (assume 50% volume is occupied by packing media)
        [M_LDPE_AF_kg, M_HDPE_AF_kg] = Packing_Media (V_m_AF); % Call packing media function
        [~, ~, P_input_LIFT, M_SS_LIFT] = Lift_Pumping_AF(Q_mgd, N_AF, D_AF); % Call lift pumping function, ~ are Q_LIFT_mgd, NP_LIFT
        if ic==1
            [M_Membrane_kg, A_LU,V_membrane_displacement] = Hollow_Fiber(N_train, Cas_per_tank, Mod_per_cas,Module_SA); % Call membrane function (hollow fiber)
        elseif ic==2
            [M_Membrane_kg, A_LU,V_membrane_displacement] = Flat_Sheet_Submerged(N_train, Cas_per_tank, Mod_per_cas,Module_SA); % Call membrane function (flat sheet)
        end
        V_membrane_tank = N_train * (D_membrane_tank * L_membrane_tank * W_membrane_tank-V_membrane_displacement); % [ft^3] Volume of reactor
        %         HRT_membrane_tank = V_membrane_tank/(Q_mgd * 133681 / 24); % [hr];

        [~, ~, P_input_PERM, M_SS_PERM] = Permeate_Pumping_Submerged(Q_mgd, N_train, D_AF, L_membrane_tank, W_N_trains, TMP,ia1); % Call permeate pumping function, ~ are Q_PERM_mgd, NP_PERM

        if ia1~=2
            freq = UParameter(6); % freq: Gas sparging frequency [sparging time/total time]
            [~, P_input_blower, M_SS_gas, ~, ~] = Gas_Sparging_Submerged(N_train, Cas_per_tank, A_LU, SGD, L_membrane_tank, W_PB, L_BB, freq); % Call gas sparging function, ~ are Q_gas_cfm, OD_gh, OD_gsm
        end
        Q_R_mgd = 0;

    elseif ib==2
        [V_m_AF, D_AF, Dia_AF, N_AF, VWC_AF, VSC_AF, VEX_AF] = AF_Crossflow(Q_mgd, S_SO, X_SO, OL_AF, HL_AF, R_AF); % Call anaerobic filter function
        V_rctr = N_AF * V_m_AF / 0.5; % [ft^3] Volume of reactor (assume 50% volume is occupied by packing media)
        [M_LDPE_AF_kg, M_HDPE_AF_kg] = Packing_Media (V_m_AF); % Call packing media function
        [~, ~, P_input_LIFT, M_SS_LIFT] = Lift_Pumping_AF(Q_mgd, N_AF, D_AF); % Call lift pumping function, ~ are Q_LIFT_mgd, NP_LIFT
        if ic==3
            [M_Membrane_kg, Q_R_mgd] = Multi_Tube(Q_mgd, J, 0.4, VEL_xflow,N_train,Cas_per_tank, Mod_per_cas,Module_SA);
            [~, ~, P_input_PERM, M_SS_PERM] = Permeate_Pumping_Cross_Flow(Q_mgd, Cas_per_tank, TMP,0,ia1); %D_train is 0 for AF, this value is used in calculation of static head loss through X-flow membranes, ~ are Q_PERM_mgd, NP_PERM
        elseif ic==2
            [M_Membrane_kg, Q_R_mgd] = Flat_Sheet_Xflow(Q_mgd, J, 0.4, VEL_xflow,N_train,Cas_per_tank, Mod_per_cas,Module_SA);
            [~, ~, P_input_PERM, M_SS_PERM] = Permeate_Pumping_Cross_Flow(Q_mgd, Cas_per_tank, TMP,0,ia1); %D_train is 0 for AF, this value is used in calculation of static head loss through X-flow membranes, ~ are Q_PERM_mgd, NP_PERM
        end
        % Call retentate pumping function
        [~, P_input_R, M_SS_R] = Retentate_Pumping_AF(Q_R_mgd, Cas_per_tank, D_AF); % ~ is NP_R
    end
    [~,~, P_input_IR, M_SS_IR] = Recirculation_Pumping_AF(Q_mgd, R_AF, N_AF, D_AF, Dia_AF,Q_R_mgd,Upflow_vel_for_GAC,L_membrane_tank,W_membrane_tank,N_train,ia1); % Call recirculation pumping function, ~ are Q_IR_mgd,NP_IR
end

%% Step A1 - Aerobic Polishing Filter/GAC
if ia1==3
    [V_m_AER, ~, ~, ~, VWC_AER, VSC_AER, VEX_AER] = AER_Filter(Q_mgd, S_SO, X_SO, OL_AER, HL_AER); % Call aerobic filter function, ~ are D_AER, Dia_AER, N_AER
    [M_LDPE_AER_kg, M_HDPE_AER_kg] = Packing_Media (V_m_AER); % Call packing media function
elseif ia1==2 && ib==1 %GAC in submerged membrane
    C_GAC = DVariable(length(DVariable)); % [g/L] or [kg/m^3] Concentration of GAC (Yoo, 2012)
    M_GAC_kg = C_GAC * V_membrane_tank * 0.0283168; % [kg/m^3] Mass of GAC
end


%% Step B, C, D - Membrane Configuration & Type & Material
% if ib==1
% %     if ic==1
% %         [M_Membrane_kg, A_LU,V_membrane_displacement] = Hollow_Fiber(Q_mgd, J, N_train, 0.4, Cas_per_tank, Mod_per_cas,Module_SA);
% %     elseif ic==2
% %         [M_Membrane_kg, A_LU,V_membrane_displacement] = Flat_Sheet_Submerged(Q_mgd, J, N_train, 0.4, Cas_per_tank, Mod_per_cas,Module_SA);
% %     end
%     [Q_PERM_mgd, NP_PERM, P_input_PERM, M_SS_PERM] = Permeate_Pumping_Submerged(Q_mgd, N_train, D_train, L_membrane_tank, W_N_trains, TMP); % Call permeate pumping function
%
%     if ia1==1
%         freq = UParameter(6); % freq: Gas sparging frequency [sparging time/total time]
%         [Q_gas_cfm, NB1, P_input_blower, M_SS_gas, OD_gh, OD_gsm] = Gas_Sparging_Submerged(N_train, Cas_per_tank, A_LU, SGD, L_membrane_tank, W_PB, L_BB, freq); % Call gas sparging function
%     end
%
% % elseif ib==2
%     if ic==3
%         [M_Membrane_kg, N_SU, N_LU, Q_R_mgd] = Multi_Tube(Q_mgd, J, 0.4, VEL_xflow);
%         [Q_PERM_mgd, NP_PERM, P_input_PERM, M_SS_PERM] = Permeate_Pumping_Cross_Flow(Q_mgd, N_LU, TMP);
%     elseif ic==2
%         [M_Membrane_kg, N_SU, N_LU, Q_R_mgd] = Flat_Sheet_Xflow(Q_mgd, J, 0.4, VEL_xflow);
%         [Q_PERM_mgd, NP_PERM, P_input_PERM, M_SS_PERM] = Permeate_Pumping_Cross_Flow(Q_mgd, N_LU, TMP);
%     end
%
%     % Call retentate pumping function
%     if ia==2 && ic==3 || ia==2 && ic==2
%         [NP_R, P_input_R, M_SS_R] = Retentate_Pumping_AF(Q_R_mgd, N_LU, D_AF);
%     elseif ia==1 && ic==3 || ia==1 && ic==2
%         [NP_R, P_input_R, M_SS_R] = Retentate_Pumping_CSTR(Q_R_mgd, N_LU, D_train);
%     end
% end


%% Step G - Chemical Cleaning
[M_NaOCl_kg, Q_NaOCl_weekly, M_CA_kg, Q_CA_weekly] = Chemical_Cleaning (Q_mgd, Year); % Call chemical cleaning function

% NaOCl Pumping
[~, ~, P_input_NaOCl, M_SS_NaOCl, M_HDPE_NaOCl_kg] = Chemical_Pumping(Q_NaOCl_weekly); % Call chemical pumping function, ~ are Q_NaOCl_mgd, NP_NaOCl

% Citric acid Pumping
[~, ~, P_input_CA, M_SS_CA, M_HDPE_CA_kg] = Chemical_Pumping(Q_CA_weekly); % Call chemical pumping function, ~ are Q_CA_mgd, NP_CA

M_HDPE_CHEM_kg = M_HDPE_NaOCl_kg + M_HDPE_CA_kg;
M_SS_CHEM = M_SS_NaOCl + M_SS_CA;


%% Step H&I - Methane Processing
Q_CH4_per_kg_COD = 0.2; % [m^3-CH4/kg-COD removed] Unit CH4 production rate
Q_CO2_per_kg_COD = 0.108; % [m^3-CO2/kg-COD removed] Unit CO2 production rate, based on assumption that methane is 65% of gas produced and CO2 is 35%

Q_cmd = Q_mgd * 3785.41178; % [m^3/day] Unit conversion
M_COD_removed = (COD_inf - COD_eff)/10^3 * Q_cmd; % [kg-COD/day] Daily COD removal

if ih==1
    Q_CH4 = Q_CH4_per_kg_COD * M_COD_removed; % [m^3-CH4/day], assuming DM is 100% efficient
    n_DM = Q_cmd/(24*30); % [unitless] Calculates number of DMs needed assuming each can handle 30 m3/hr (from spec sheet)
    E_DM = 3*n_DM; % [kW] Calculates total power of all DMs  (assumes 3 kW/DM) (NOT: 42 W/DM, Bandara et al., 2011)
    Q_CH4_1 = Q_CH4_per_kg_COD * M_COD_removed;
else
    Q_CH4 = Q_CH4_per_kg_COD * M_COD_removed*UParameter(8); % [m^3-CH4/day]
    Q_CH4_1 = Q_CH4_per_kg_COD * M_COD_removed;
    CH4_eff = (Q_CH4_1-Q_CH4)*0.668*1000/(Q_mgd*3785.41178); %[g/m3]; Methane released from effluent; Density of methane = 0.668 kg/m3
end

[P_offset] = CHP(Q_CH4); % Call combined heat and power function
E_offset_kWh = P_offset * Year * 365 * 24; % [kWh] Total electricity offset over N years
Q_CO2 = Q_CO2_per_kg_COD*M_COD_removed+Q_CH4*1.49; %[m^3-CO2/day], accounting for CO2 emitted by burning methane (1.49 m3 CO2 per m3 CH4 burned)
CO2_eff = (Q_CO2*1.842*1000/(Q_mgd*3785.41178))*0.25; %[g/m3]; Density of CO2 = 1.842 kg/m3, assume 25% of carbon in WW is non-biogenic according to Griffith, D.R. et al. 2009, ES&T, vol 43, issue 15, 5647-5651

%% Sludge Handling
C_rctr = 8500; % [mg/L] MLSS concentration in the reactor (number suggested by H&S)
C_WAS = 10500; % [mg/L] MLSS concentration in the WAS flow (number suggested by H&S)
yield = UParameter(7); % [day] H&S suggested 20 - 50 day SRT
M_biomass_wastage = (COD_inf*yield*Q_cmd)/1000; %[kg/d], 0.03 is the assumed yield of kg biomass/kg COD degraded
Q_WAS = M_biomass_wastage/((C_WAS*28.3168)/10^6); % [ft^3/day] 28.3168 L/ft3
Q_WAS_mgd = Q_WAS * 7.48052 / 10^6; % [mgd] unit conversion
[P_input_sludge, M_SS_GBT, M_COD_WAS] = Sludge_Handling (Q_CH4_1, Q_WAS_mgd, M_COD_removed);


%% Life Cycle Inventory Summary
%% Construction Phase
% Concrete
if ia==1
    VSC = VSC_CSTR;
    VWC = VWC_CSTR;
    VEX = VEX_CSTR;
elseif ia==2 && ia1~=3
    VSC = VSC_AF;
    VWC = VWC_AF;
    VEX = VEX_AF;
elseif ia==2 && ia1==3
    VSC = VSC_AF + VSC_AER;
    VWC = VWC_AF + VWC_AER;
    VEX = VEX_AF + VEX_AER;
end

VC_m3 = (VSC + VWC) * 0.0283168; % [m^3] Unit conversion from [ft^3] to [m^3]

% Excavation
VEX_m3 = VEX * 0.0283168; % [m^3] Unit conversion from [ft^3] to [m^3]

% GAC
if ia1~=2
    M_GAC_kg = 0;
end

% Stainless Steel and  LDPE & HDPE
if ia==1
    M_LDPE_kg = 0;
    M_HDPE_kg = M_HDPE_CHEM_kg;
    if ib==1 && ia1==1
        M_Steel_kg = (M_SS_PERM + M_SS_IR + M_SS_gas + M_SS_CHEM + M_SS_GBT);
    elseif ib==1 && ia1==2
        M_Steel_kg = (M_SS_PERM + M_SS_IR + M_SS_CHEM + M_SS_GBT);
    elseif ib==2
        M_Steel_kg = (M_SS_PERM + M_SS_IR + M_SS_R + M_SS_CHEM + M_SS_GBT);
    end

elseif ia==2
    if ib==1 && ia1~=2
        M_Steel_kg = (M_SS_LIFT + M_SS_PERM + M_SS_IR + M_SS_gas + M_SS_CHEM + M_SS_GBT);
    elseif ib==1 && ia1==2
        M_Steel_kg = (M_SS_LIFT + M_SS_PERM + M_SS_IR + M_SS_CHEM + M_SS_GBT);
    elseif ib==2
        M_Steel_kg = (M_SS_LIFT + M_SS_PERM + M_SS_IR + M_SS_R + M_SS_CHEM + M_SS_GBT);
    end
    if ia1==3
        M_LDPE_kg = M_LDPE_AF_kg + M_LDPE_AER_kg;
        M_HDPE_kg = M_HDPE_AF_kg + M_HDPE_AER_kg + M_HDPE_CHEM_kg;
    else
        M_LDPE_kg = M_LDPE_AF_kg;
        M_HDPE_kg = M_HDPE_AF_kg + M_HDPE_CHEM_kg;
    end

end

% Call the "LCI_CON" function to generate construction inventory matrix
INV_CON = LCI_CON (VC_m3, VEX_m3, M_Steel_kg, M_Membrane_kg, M_LDPE_kg, M_HDPE_kg, M_GAC_kg, id);


%% Operational Phase
% Total Energy Consumption
if ih==1
    if ia==1
        if ia1==1 && ib==1
            P_input_tot = P_input_PERM + P_input_IR + P_input_blower + P_input_NaOCl + P_input_CA + P_input_sludge+E_DM; % [Kw] Total power input
            Power_pct = [P_input_PERM/P_input_tot, P_input_IR/P_input_tot, ...
                (P_input_NaOCl + P_input_CA)/P_input_tot, P_input_sludge/P_input_tot,P_input_blower/P_input_tot,E_DM/P_input_tot];
        elseif ib==2
            P_input_tot = P_input_PERM + P_input_IR + P_input_R + P_input_NaOCl + P_input_CA + P_input_sludge+E_DM; % [Kw] Total power input
            Power_pct = [P_input_PERM/P_input_tot, P_input_IR/P_input_tot, ...
                (P_input_NaOCl + P_input_CA)/P_input_tot, P_input_sludge/P_input_tot,P_input_R/P_input_tot,E_DM/P_input_tot];
        elseif ia1==2 && ib==1
            P_input_tot = P_input_PERM + P_input_IR + P_input_NaOCl + P_input_CA + P_input_sludge+E_DM; % [Kw] Total power input
            Power_pct = [P_input_PERM/P_input_tot, P_input_IR/P_input_tot, ...
                (P_input_NaOCl + P_input_CA)/P_input_tot, P_input_sludge/P_input_tot,E_DM/P_input_tot];
        end
    elseif ia==2
        if ia1~=2 && ib==1
            P_input_tot = P_input_PERM + P_input_IR + P_input_blower + P_input_NaOCl + P_input_CA + P_input_sludge+P_input_LIFT+E_DM; % [Kw] Total power input
            Power_pct = [P_input_PERM/P_input_tot, P_input_IR/P_input_tot, ...
                (P_input_NaOCl + P_input_CA)/P_input_tot, P_input_sludge/P_input_tot,P_input_blower/P_input_tot,P_input_LIFT/P_input_tot,E_DM/P_input_tot];
        elseif ib==2
            P_input_tot = P_input_PERM + P_input_LIFT + P_input_IR + P_input_R + P_input_NaOCl + P_input_CA + P_input_sludge+E_DM; % [Kw] Total power input
            Power_pct = [P_input_PERM/P_input_tot, P_input_IR/P_input_tot,...
                (P_input_NaOCl + P_input_CA)/P_input_tot, P_input_sludge/P_input_tot, P_input_LIFT/P_input_tot,P_input_R/ P_input_tot,E_DM/P_input_tot];
        elseif ia1==2 && ib==1
            P_input_tot = P_input_PERM + P_input_LIFT + P_input_IR + P_input_NaOCl + P_input_CA + P_input_sludge+E_DM; % [Kw] Total power input
            Power_pct = [P_input_PERM/P_input_tot, P_input_IR/P_input_tot, ...
                (P_input_NaOCl + P_input_CA)/P_input_tot, P_input_sludge/P_input_tot, P_input_LIFT/P_input_tot,E_DM/P_input_tot];
        end
    end
else
    if ia==1
        if ia1==1 && ib==1
            P_input_tot = P_input_PERM + P_input_IR + P_input_blower + P_input_NaOCl + P_input_CA + P_input_sludge; % [Kw] Total power input
            Power_pct = [P_input_PERM/P_input_tot, P_input_IR/P_input_tot, ...
                (P_input_NaOCl + P_input_CA)/P_input_tot, P_input_sludge/P_input_tot,P_input_blower/P_input_tot];
        elseif ib==2
            P_input_tot = P_input_PERM + P_input_IR + P_input_R + P_input_NaOCl + P_input_CA + P_input_sludge; % [Kw] Total power input
            Power_pct = [P_input_PERM/P_input_tot, P_input_IR/P_input_tot, ...
                (P_input_NaOCl + P_input_CA)/P_input_tot, P_input_sludge/P_input_tot,P_input_R/P_input_tot];
        elseif ia1==2 && ib==1
            P_input_tot = P_input_PERM + P_input_IR + P_input_NaOCl + P_input_CA + P_input_sludge; % [Kw] Total power input
            Power_pct = [P_input_PERM/P_input_tot, P_input_IR/P_input_tot, ...
                (P_input_NaOCl + P_input_CA)/P_input_tot, P_input_sludge/P_input_tot];
        end
    elseif ia==2
        if ia1~=2 && ib==1
            P_input_tot = P_input_PERM + P_input_IR + P_input_blower + P_input_NaOCl + P_input_CA + P_input_sludge+P_input_LIFT; % [Kw] Total power input
            Power_pct = [P_input_PERM/P_input_tot, P_input_IR/P_input_tot, ...
                (P_input_NaOCl + P_input_CA)/P_input_tot, P_input_sludge/P_input_tot,P_input_blower/P_input_tot,P_input_LIFT/P_input_tot];
        elseif ib==2
            P_input_tot = P_input_PERM + P_input_LIFT + P_input_IR + P_input_R + P_input_NaOCl + P_input_CA + P_input_sludge; % [Kw] Total power input
            Power_pct = [P_input_PERM/P_input_tot, P_input_IR/P_input_tot,...
                (P_input_NaOCl + P_input_CA)/P_input_tot, P_input_sludge/P_input_tot, P_input_LIFT/P_input_tot,P_input_R/ P_input_tot];
        elseif ia1==2 && ib==1
            P_input_tot = P_input_PERM + P_input_LIFT + P_input_IR + P_input_NaOCl + P_input_CA + P_input_sludge; % [Kw] Total power input
            Power_pct = [P_input_PERM/P_input_tot, P_input_IR/P_input_tot, ...
                (P_input_NaOCl + P_input_CA)/P_input_tot, P_input_sludge/P_input_tot, P_input_LIFT/P_input_tot];
        end
    end
end

E_input_kWh = P_input_tot * Year * 365 * 24; % [kWh] Total electricity consumption over N years

INV_OP = [0; 0; 0; M_NaOCl_kg; E_input_kWh; 0; 0; M_COD_WAS; M_CA_kg]; % Operational inventory matrix


%% Direct Emissions
% Total emissions over N years
COD_water = COD_eff /10^3 * V_treated; % [kg-COD]
NH3_water = NH3_eff /10^3 * V_treated; % [kg-N]
NH4_water = NH4_eff /10^3 * V_treated; % [kg-N]
Org_N_water = Org_N_eff /10^3 * V_treated; % [kg-N]
P_water = P_eff /10^3 * V_treated; % [kg-N]
PO4_water = PO4_eff /10^3 * V_treated; % [kg-N]
CH4_air = CH4_eff /10^3 * V_treated; % [kg-CH4]
CO2_air = CO2_eff /10^3 * V_treated; % [kg-CO2]

DEmis = [COD_water; NH3_water; NH4_water; Org_N_water; P_water; PO4_water; CH4_air; CO2_air];


%% Output for Cost Estimation
%NOTE: Equations are either from H&S (indicated in parentheses) or from
%Capdet (all others). Some equations from Capdet have been modified to be
%smooth.
Output_LCC=struct();
total_OM=0;
total_const=0;
total=0;

Output_LCC.Air_Piping=0;
Output_LCC.Blowers=0;

if ia1~=2 && ib==1
    %Air piping (CapdetWorks, unchanged and H&S)
    CFMD = Air_demand_instantaneous; % Design Capacity of Blowers, scfm
    if (CFMD<=1000)
        COSTAP = 617.2 * 3.33*CFMD^0.2553; % 3.33: Air Flow Fraction, calculated as STE/6, and STE stands for Standard Oxygen Transfer Efficiency _%
    elseif (CFMD>1000)&&(CFMD<=10000)       % The default STE value is 20. If users have their own STE value, then plug the value into STE/6.
        COSTAP = 1.43 * 3.33* CFMD^1.1337;  % If STE/6 > 1, take this value and replace 3.33, if the value is smaller than 1, then use 1 to replace 3.33.
    elseif (CFMD>10000)
        COSTAP = 28.59 * 3.33*CFMD^0.8085; % In calculating COSTAP, the CEPCIP (Current CE Plant Cost Index for Pipe, Valves, etc) value is set to be 241
    end
    Output_LCC.Air_Piping=COSTAP;
    total_const=total_const+COSTAP;

    %Blowers (CapdetWorks, changed and H&S)
    %     NB=2;
    %     BLW_size =ceil(Air_demand_10_10_instantaneous/(NB-1)); %Blower size to handle 10:10
    %     Tot_cfm=BLW_size*(NB-1);
    %     while Tot_cfm<Air_demand_10_10_instantaneous
    %         NB=NB+1;
    %         Tot_cfm=BLW_size*(NB-1);
    %     end
    if (TCFM<=30000)
        COSTRO = 0.7*CFMB^0.6169; %These equations are changed to smooth out the curve
    elseif (TCFM>30000)&&(TCFM<=72000)
        COSTRO = 0.377* CFMB^0.5928;
    elseif (TCFM>72000)
        COSTRO = 0.964* CFMB^0.4286;
    end

    if (TCFM<=30000)
        COSTBE  = 58000*(COSTRO) / 100; % Purchase Cost of Blower of CFMB Capacity, $, 50000:Purchase cost of 3,000 scfm at 8 psig
    elseif (TCFM>30000)&&(TCFM<=72000)
        COSTBE  = 218000*(COSTRO) / 100; % Purchase Cost of 12,000 scfm at 8 psig
    elseif (TCFM>72000)
        COSTBE  = 480000*(COSTRO) / 100; % Purchase Cost of 50,000 scfm at 8 psig
    end

    IBC = 2 * COSTBE; % IBC: Installed Blower Costs, $
    BBA = 128* TCFM^0.256; % Blower Building Area, ft^2
    BBC = BBA * 90; % 90: Unit Price Input for Building Costs, $/ft^2
    TBCC_BLW = (IBC * NB*2 + BBC) * 1.11; % 1.11: Correction Factor for Other Minor Construction Costs includes piping, concrete, steel, electrical, paint and installation labor. Multiply NB by 2 to represent replacement after 15 years (Daigger, 2011).
    Output_LCC.Blowers=TBCC_BLW;
    total_const=total_const+TBCC_BLW;
end

%Chemical Feed System (Hazen and Sawyer)
Usage_rate_NaOCl = 2200; %[gal/yr/MGD]
Chem_cost_NaOCl = 0.54; %[$/gal], 12.5% solution
Annual_cost_NaOCl = Chem_cost_NaOCl*Q_mgd*Usage_rate_NaOCl*(12.5/15); %[$]
Usage_rate_citric = 600; %[gal/yr/MGD]
Chem_cost_citric = 0.85; %[$/gal], 100% solution, 13.8 lb/gal
Annual_cost_citric = Chem_cost_citric*Q_mgd*Usage_rate_citric*(13.8); %[$]
Usage_rate_Bisulfite = 350; %[gal/yr/MGD]
Chem_cost_Bisulfite = 0.3; %[$/gal], 38% solution, 3.5 lb/gal
Annual_cost_Bisulfite = Chem_cost_Bisulfite*Q_mgd*Usage_rate_Bisulfite*3.5*0.38; %[$]
Cost_cleaning = Annual_cost_NaOCl+Annual_cost_citric+Annual_cost_Bisulfite;
Output_LCC.Chemical_Cleaning=Cost_cleaning;
total_OM=total_OM+Cost_cleaning;

%Pumping (CapdetWorks, changed)
RSR = 4;%input('Enter the Return Sludge Ratio to Average Wastewater Flow: ');
Qavg=Q_mgd;
GPMI = (2*Qavg*10^6) / 1440; % GPMI: Design Capacity of Intermediate Pumps in gpm. 2 = excess capacity factor to handle peak flows
GPMR = (Qavg*RSR*10^6) / 1440; % GPMR: Design Capacity of Return Sludge Pumps in gpm
FPCI = 1440*GPMI / (10^6); % FPC: Firm Pumping Capacity in mgd
FPCR = 1440*GPMR / (10^6);
if (FPCI == 0)
    OLCI = 0;
elseif (FPCI > 0) && (FPCI <= 7)
    OLCI = 440*25* FPCI^0.1285;
elseif (FPCI > 7) && (FPCI <= 41) %Changed the upper limit to 41 (from 30)
    OLCI = 294.4*25* FPCI^0.3335;
elseif (FPCI > 41) && (FPCI <= 80)
    OLCI = 40.5 *25*FPCI^0.8661;
elseif (FPCI > 80)
    OLCI = 21.3*25*FPCI^1.012;
end
if (FPCI == 0)
    MLCI = 0;
elseif (FPCI > 0) && (FPCI <= 7)
    MLCI = 360 *25*FPCI^0.1478;
elseif (FPCI > 7) && (FPCI <= 41) %Changed the upper limit to 41 (from 30)
    MLCI = 255.2* 25*FPCI^0.3247;
elseif (FPCI > 41) && (FPCI <= 80)
    MLCI = 85.7* 25*FPCI^0.6456;
elseif (FPCI > 80)
    MLCI = 30.6*25* FPCI^0.8806;
end
if (FPCR == 0)
    OLCR = 0;
elseif (FPCR > 0) && (FPCR <= 7)
    OLCR = 440* 25*FPCR^0.1285;
elseif (FPCR > 7) && (FPCR <= 41) %Changed the upper limit to 41 (from 30)
    OLCR = 294.4* 25*FPCR^0.3335;
elseif (FPCR > 41) && (FPCR <= 80)
    OLCR = 40.5 *25*FPCR^0.8661;
elseif (FPCR > 80)
    OLCR = 21.3*25*FPCR^1.012;
end
if (FPCR == 0)
    MLCR = 0;
elseif (FPCR > 0) && (FPCR <= 7)
    MLCR = 360 *25*FPCR^0.1478;
elseif (FPCR > 7) && (FPCR <= 41) %Changed the upper limit to 41 (from 30)
    MLCR = 255.2*25* FPCR^0.3247;
elseif (FPCR > 41) && (FPCR <= 80)
    MLCR = 85.7*25* FPCR^0.6456;
elseif (FPCR > 80)
    MLCR = 30.6*25* FPCR^0.8806;
end
TOLC = OLCI + OLCR;
TMLC = MLCI + MLCR;
a = 1;
GPMBI = GPMI/a;
while (GPMBI > 80000)
    a = a+1;
    GPMBI = GPMI/a;
end
if (GPMR == 0)
    b = 0;
elseif (GPMR>0)
    b = 1;
    GPMBR = GPMR/b;
    while (GPMBR > 80000)
        b = b+1;
        GPMBR = GPMR/b;
    end
end
c=2;
GPMPI = GPMBI / c;
while (GPMPI >20000)
    c=c+1;
    GPMPI = GPMBI / c;
end
c=c+1;
if (GPMR == 0)
    d = 0;
elseif (GPMR>0)
    d = 2;
    GPMPR = GPMBR/d;
    while (GPMPR > 20000)
        d = d+1;
        GPMPR = GPMR/d;
    end
    d = d+1;
end
PBAI = (0.0284*GPMBI+640)*a;
PBAR = (0.0284*GPMBR+640)*b;
COSTPB=(PBAI+PBAR)*90; % 90: Unit price for pump building cost. Information from CapdetWorks 2007.

if (GPMPI>0)&&(GPMPI<5000)
    COSTROI=2.93*GPMPI^0.4404;
elseif (GPMPI>5000)
    COSTROI=0.0064*GPMPI^1.16;
end
if (GPMPR>0)&&(GPMPR<5000)
    COSTROR=2.93*GPMPR^0.4404;
elseif (GPMPR>5000)
    COSTROR=0.0064*GPMPR^1.16;
end
IPC = 2.065*10^5+7.721*10^4*Qavg; % Fit a curve to the data instead of using the equation below
%IPC = 2*(COSTROI/100) *48549.837*a*c + 2*(COSTROR/100) *48549.837*b*d; % IPC: Installed equipment cost. 48549.837:Standard Cost of a 3000 gpm pump and driver.
COSTE = VEX/27*0.3; %Cost of earthwork
TBCC_PUP = (COSTE + COSTPB + IPC*2)*1.18; %Multiply cost of pumps by 2 to represent replacement cost after 15-year useful life is up (Diagger, 2011).
MSC = TBCC_PUP*(0.007/100);
TBCC_PUP=TBCC_PUP+MSC;
Output_LCC.Pumping_Op_Labor=TOLC;
Output_LCC.Pumping_Main_Labor=TMLC;
Output_LCC.Pumps=TBCC_PUP;
total_const=total_const+TBCC_PUP+MSC;
total_OM=total_OM+TOLC+TMLC;

%Concrete/earthwork (CapdetWorks, unchanged)
COSTCW = VWC/27 * 650;  % Cost of wall concrete. 18.52: Unit price for wall concrete in dollar, information obtained from the 2007 vendor data in CapdetWorks
% 27: Conversion Factor from ft^3 to yd^3
COSTCS = VSC/27*350;    % Cost of slab concrete. 12.96: Unit price for slab concrete in dollar, information obtained from the 2007 vendor data in CapdetWorks
% 27: Conversion Factor from ft^3 to yd^3
COSTRC = COSTCW + COSTCS; % Total concrete cost
COSTE = VEX/27*8.00; %Cost of earthwork
Output_LCC.Concrete=COSTRC;
Output_LCC.Earthwork=COSTE;
total_const=total_const+COSTRC+COSTE;

%GAC
Cost_GAC = M_GAC_kg * 6.20*2.20462; % $6.20/lb GAC, 2.20462 lb/kg (CapdetWorks 4)
Output_LCC.GAC=Cost_GAC;
total_const=total_const+Cost_GAC;

%Packing Media
if ia==2
    if ia1==3
        Cost_Packing_Media = 5.5*(V_m_AF+V_m_AER); %$5.50/ft3 packing media
        Output_LCC.Packing_Media=Cost_Packing_Media;
        total_const=total_const+Cost_Packing_Media;
    else
        Cost_Packing_Media = 5.5*(V_m_AF);
        Output_LCC.Packing_Media=Cost_Packing_Media;
        total_const=total_const+Cost_Packing_Media;
    end
else
    Output_LCC.Packing_Media=0;
end

%Handrail (CapdetWorks, unchanged)
%LHR = input('Enter a value for Length of Handrails in place in ft: ');
%COSTHR = LHR * 75; % 75: Unit price for handrail in place in dollar, information obtained from the 2007 vendor data in CapdetWorks

%Electricity (Hazen and Sawyer)
% Fine_screens = 1*Q_mgd; %at 1 hp/MGD
% Est_TDH_perm = 40; %[ft]
% Est_eff_perm = 0.8;
% Avg_bhp_perm_pumps = (Q_mgd*1000000/1440)*Est_TDH_perm/(Est_eff_perm*3960);
% HP_subtotal = Fine_screens+Avg_bhp_perm_pumps;
%
% if strcmp(step_A1, 'none') && strcmp(step_B, 'Submerged')
%     Est_TDH_blower = 6; %[psig]
%     Est_eff_blower = 0.7;
%     Avg_bhp_mem_blowers = (CFMD*0.23*(((14.7+Est_TDH_blower)/14.7)^0.283-1))/Est_eff_blower; %Note, CFMD is actually something different in H&S''s code!
%     HP_subtotal = HP_subtotal+Avg_bhp_mem_blowers;
% end
%
% Misc_power = HP_subtotal*0.015;
% HP_total = HP_subtotal+Misc_power;
% Total_avg_power = HP_total*0.746*24;
Power_cost = 0.10; %$/kWh
% Daily_power_cost = Total_avg_power*Power_cost;
Cost_power = P_input_tot * 365 * 24*Power_cost; %[$]
Output_LCC.Electricity(1,1)=Cost_power;
total_OM=total_OM+Cost_power;

%Screenings disposal (additional screenings for MBR, i.e., in addition to 6
%mm screens) (Hazen and Sawyer)
Screen_rate_max = 4; %[CF/hr/MGD]
Screen_rate = 0.5;
Screen_per_day = Screen_rate_max*Q_mgd*24*Screen_rate; %[CF/day]
Compaction = 0.75;
Daily_screen = Screen_per_day*(1-Compaction);
Daily_screen_CY = Daily_screen/27;
Cost_disposal = 225; %$/20 CY container
Cost_screening_and_sludge = Daily_screen_CY*365*Cost_disposal/20+(M_COD_WAS/907.18474)*125*365; %$125/ton dewatered sludge disposal, 907.18474 conversion of kg to ton (US)
Output_LCC.Screening=Cost_screening_and_sludge;
total_OM=total_OM+Cost_screening_and_sludge;

%Membrane replacement cost (Hazen and Sawyer)
Cost_per_sf = UParameter(11); %[$/sf]
Total_cost_rep = Cost_per_sf*1.15; %0.15* cost of membrane for replacement labor
Expected_membrane_life = UParameter(10); %[yrs]
% N_modules = Mod_per_cas*N_cassettes;
Total_SA = Contact_SA*Cas_per_tank*N_train; %[ft2]
Membrane_replacement_cost = Total_SA*Total_cost_rep; %[$]
Yearly_cost_replacement = Membrane_replacement_cost*(interest/(((1+interest)^Expected_membrane_life)-1)); % [$/yr]
Output_LCC.Memb_replacement=Yearly_cost_replacement;
total_OM=total_OM+Yearly_cost_replacement;

%DM costs
if ih==1
    Cost_DM = 10000*n_DM; %Assumes these don''t need to be replaced
    Output_LCC.Degassing_memb=Cost_DM;
    total_const = total_const + Cost_DM;
else
    Output_LCC.Degassing_memb=0;
end

%CHP costs [From Woer et al., 2012 (Range given, took midpoint)]
% if strcmp(step_I,'IC')
%     Cost_CHP_const = P_offset * 1032.5; %[$]
%     Cost_CHP_OM = P_offset * 365 * 24*.0175; %[$/yr]
% elseif strcmp(step_I,'CG')
%     Cost_CHP_const = P_offset * 1550;
%     Cost_CHP_OM = P_offset * 365 * 24*.011;
% elseif strcmp(step_I,'micro')
Cost_CHP_const = P_offset * 1225;
Cost_CHP_OM = P_offset * 365 * 24*.0185;
% else
%     Cost_CHP_const = P_offset * 4540;
%     Cost_CHP_OM = P_offset * 365 * 24*.0115;
% end
Output_LCC.CHP_Const=Cost_CHP_const;
Output_LCC.CHP_OM=Cost_CHP_OM;
total_const = total_const+Cost_CHP_const;
total_OM = total_OM+Cost_CHP_OM;

% Net present value of O&M costs
Present_OM = total_OM*(((1+interest)^Year-1)/(interest*(1+interest)^Year)); %From Engineering Economic Analysis, 11th edition, Newnan et al.

%Total Cost
total=Present_OM+total_const;
Output_LCC.Total_Const=total_const;
Output_LCC.Total_Present_OM=Present_OM;
Output_LCC.Total_Yearly_OM=total_OM;
Output_LCC.Total=total;
Output_LCC.Percent_Const=total_const./total;
Output_LCC.Percent_OM=Present_OM./total;
