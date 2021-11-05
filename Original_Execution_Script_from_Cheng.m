clc; clear

%% Decision variables
J = 8.5; % [L/m^2-hr]

TMP_high = 25; % [psi]

TMP_low = 20; % [psi]

OL_AF = 5; % [g-COD/L-d] or [kg/m^3-day]

HL_AF = 10; % [m/hr]

OL_AER = 2; % [g-COD/L-d] or [kg/m^3-day]

HL_AER = 2; % [m/hr]

HRT = 10; % [hr] 

N_train = 7;

IRR = 4;

% Assign values to decision variables
%if strcmp(step_A, 'CSTR') && strcmp(step_B, 'Submerged') && strcmp(step_A1, 'none')
     DVariable_1 = [J, TMP_high, N_train, IRR, HRT];
%elseif strcmp(step_A, 'CSTR') && strcmp(step_B, 'Submerged') && strcmp(step_A1, '+GAC')
     DVariable_2 = [J, TMP_low, N_train, IRR, HRT];
%elseif strcmp(step_A, 'CSTR') && strcmp(step_B, 'Cross-flow') && strcmp(step_A1, 'none')
     DVariable_3 = [J, TMP_high, N_train, IRR, HRT];
%elseif strcmp(step_A, 'CSTR') && strcmp(step_B, 'Cross-flow') && strcmp(step_A1, '+GAC')
     DVariable_6 = [J, TMP_low, N_train, IRR, HRT];    
%elseif strcmp(step_A, 'AF') && strcmp(step_A1, 'none') 
     DVariable_4 = [J, TMP_high, OL_AF, HL_AF, IRR];
%elseif strcmp(step_A, 'AF') && strcmp(step_A1, '+GAC')
     DVariable_7 = [J, TMP_low, OL_AF, HL_AF, IRR];
%elseif strcmp(step_A, 'AF') && strcmp(step_A1, '+AER')
     DVariable_5 = [J, TMP_low, OL_AF, HL_AF, IRR, OL_AER, HL_AER];
%end

%% Uncertainty parameters
Q_mgd = 20; % [mgd]

S_SO = 300; % [mg-COD/L] or [g/m^3] 

X_SO = 100; % [mg-COD/L] or [g/m^3] 

UParameter = [Q_mgd, S_SO, X_SO];

%% System 1: CSTR + submerged (hollow fiber)
disp('     System 1'); 

%% Inventory Analysis
[INV_CON_1, INV_OP_1, DEmis_1, Power_pct_1, E_input_kWh_1, E_offset_kWh_1, V_treated_1, Output_cost_1] = LCI_anMBR('CSTR', 'none', 'Submerged', 'HF', 'PET', DVariable_2, UParameter);

%% Impact Assessment
[IMPACT_CON_1, IMPACT_OP_1, IMPACT_DE_1, IMPACT_avoided_1] = Impact_Assessment (INV_CON_1, INV_OP_1, DEmis_1, E_offset_kWh_1, V_treated_1);

% Total Impact (construction + operation)
CON_TOT_1 = sum(IMPACT_CON_1);
OP_TOT_1 = sum(IMPACT_OP_1); 
IMPACT_TOT_1 = sum(IMPACT_CON_1) + sum(IMPACT_OP_1) + sum(IMPACT_DE_1);

%% Impact by Sources
% Construction phase
Concrete_pct_1 = (IMPACT_CON_1(1,:)) ./ IMPACT_TOT_1;
Steel_pct_1 = (IMPACT_CON_1(2,:) + IMPACT_CON_1(8,:)) ./ IMPACT_TOT_1;
Excavation_pct_1 = (IMPACT_CON_1(18,:)) ./ IMPACT_TOT_1;
CON_MISC_pct_1 = (sum(IMPACT_CON_1(3:7,:)) + sum(IMPACT_CON_1(9:17,:)) + sum(IMPACT_CON_1(19:39,:))) ./ IMPACT_TOT_1;


% Operation phase
Pumping_PERM_pct_1  = (Power_pct_1(1) * sum(IMPACT_OP_1(7:14,:))) ./ IMPACT_TOT_1;
Pumping_IR_pct_1 = (Power_pct_1(2) * sum(IMPACT_OP_1(7:14,:))) ./ IMPACT_TOT_1;
Sparging_pct_1 = (Power_pct_1(3) * sum(IMPACT_OP_1(7:14,:))) ./ IMPACT_TOT_1;
Pumping_CHEM_pct_1 = (Power_pct_1(4) * sum(IMPACT_OP_1(7:14,:))) ./ IMPACT_TOT_1;
Pumping_GBT_pct_1 = (Power_pct_1(5) * sum(IMPACT_OP_1(7:14,:))) ./ IMPACT_TOT_1;
OP_MISC_pct_1 = (sum(IMPACT_OP_1(1:6,:)) + sum(IMPACT_OP_1(15:20,:))) ./ IMPACT_TOT_1;

Energy_TOT_1 = sum(IMPACT_OP_1(7:14,:));

% Direct emission
DE_pct_1 = sum(IMPACT_DE_1) ./ IMPACT_TOT_1;

% Avoided Energy
Avoided_pct_1 = sum(IMPACT_avoided_1) ./ IMPACT_TOT_1;

Impact_source_pct_1 = [Concrete_pct_1; Steel_pct_1; Excavation_pct_1; CON_MISC_pct_1; ...
    Pumping_PERM_pct_1; Pumping_IR_pct_1; Sparging_pct_1; Pumping_CHEM_pct_1; Pumping_GBT_pct_1; OP_MISC_pct_1; DE_pct_1; Avoided_pct_1];

%Labeled Output
Labeled_Concrete_pct_1 = horzcat('Concrete',num2cell(Concrete_pct_1));
Labeled_Steel_pct_1 = horzcat('Steel',num2cell(Steel_pct_1));
Labeled_Excavation_pct_1 = horzcat('Excavation',num2cell(Excavation_pct_1));
Labeled_CON_MISC_pct_1 = horzcat('Construction Misc.',num2cell(CON_MISC_pct_1));
Labeled_Pumping_PERM_pct_1 = horzcat('Permeate Pumping',num2cell(Pumping_PERM_pct_1));
Labeled_Pumping_IR_pct_1 = horzcat('Recirculation Pumping',num2cell(Pumping_IR_pct_1));
Labeled_Sparging_pct_1 = horzcat('Gas Sparging',num2cell(Sparging_pct_1));
Labeled_Pumping_CHEM_pct_1 = horzcat('Chemical Pumping',num2cell(Pumping_CHEM_pct_1));
Labeled_Pumping_GBT_pct_1 = horzcat('Gravity Belt Thickeners',num2cell(Pumping_GBT_pct_1));
Labeled_OP_MISC_pct_1 = horzcat('Operation Misc.',num2cell(OP_MISC_pct_1));
Labeled_DE_pct_1 = horzcat('Direct Emissions',num2cell(DE_pct_1));
Labeled_Avoided_pct_1 = horzcat('Avoided Energy',num2cell(Avoided_pct_1));

Labeled_Impact_source_pct_1 = {'Type','OZO','GWM','SMG','ACD','EUT', 'CAR','NONCAR','RSP','ECOTOX'};
Labeled_Impact_source_pct_10=vertcat(Labeled_Concrete_pct_1,Labeled_Steel_pct_1,Labeled_Excavation_pct_1,...
    Labeled_CON_MISC_pct_1,Labeled_Pumping_PERM_pct_1,Labeled_Pumping_IR_pct_1,Labeled_Pumping_CHEM_pct_1,Labeled_Pumping_GBT_pct_1,...
    Labeled_OP_MISC_pct_1,Labeled_DE_pct_1,Labeled_Avoided_pct_1);
Labeled_Impact_source_pct_11=array2table(Labeled_Impact_source_pct_10,'VariableNames',Labeled_Impact_source_pct_1);
Labeled_Impact_source_pct_11.Properties.RowNames = Labeled_Impact_source_pct_11.Type;
Labeled_Impact_source_pct_11.Type = [];

%% Plot (1) Impact by sources 
barExtended(Impact_source_pct_1');
Impact_Categ = {'OZO';'GWM';'SMG';'ACD';...
    'EUT'; 'CAR';'NON-CAR';'RSP';...
    'ECOTOX'};
set(gca,'xticklabel',Impact_Categ);
legend('Concrete', 'Steel', 'Excavation', 'Construction Misc.', ...
    'Permeate Pumping', 'Recirculation Pumping', 'Gas Sparging', 'Chemical Pumping', 'Gravity belt thickeners', ...
    'Operation Misc.', 'Direction Emissions', 'Avoided Energy', 'Location', 'NorthEastOutside');


%% System 2: CSTR + submerged (hollow fiber) + GAC
disp('     System 2'); 
%% LCI
[INV_CON_2, INV_OP_2, DEmis_2, Power_pct_2, E_input_kWh_2, E_offset_kWh_2, V_treated_2] = LCI_anMBR('CSTR', '+GAC', 'Submerged', 'HF', 'PET', DVariable_2, UParameter);


%% LCI --> Impact
[IMPACT_CON_2, IMPACT_OP_2, IMPACT_DE_2, IMPACT_avoided_2] = Impact_Assessment (INV_CON_2, INV_OP_2, DEmis_2, E_offset_kWh_2, V_treated_2);

% Total Impact (construction + operation)
CON_TOT_2 = sum(IMPACT_CON_2);
OP_TOT_2 = sum(IMPACT_OP_2);  
IMPACT_TOT_2 = sum(IMPACT_CON_2) + sum(IMPACT_OP_2) + sum(IMPACT_DE_2);

%% Impact by Sources
% Construction phase
Concrete_pct_2 = (IMPACT_CON_2(1,:)) ./ IMPACT_TOT_2;
Steel_pct_2 = (IMPACT_CON_2(2,:) + IMPACT_CON_2(8,:)) ./ IMPACT_TOT_2;
Excavation_pct_2 = (IMPACT_CON_2(18,:)) ./ IMPACT_TOT_2;
CON_MISC_pct_2 = (sum(IMPACT_CON_2(3:7,:)) + sum(IMPACT_CON_2(9:17,:)) + sum(IMPACT_CON_2(19:39,:))) ./ IMPACT_TOT_2;

% Operation phase
Pumping_PERM_pct_2  = (Power_pct_2(1) * sum(IMPACT_OP_2(7:14,:))) ./ IMPACT_TOT_2;
Pumping_IR_pct_2 = (Power_pct_2(2) * sum(IMPACT_OP_2(7:14,:))) ./ IMPACT_TOT_2;
Pumping_CHEM_pct_2 = (Power_pct_2(3) * sum(IMPACT_OP_2(7:14,:))) ./ IMPACT_TOT_2;
Pumping_GBT_pct_2 = (Power_pct_2(4) * sum(IMPACT_OP_2(7:14,:))) ./ IMPACT_TOT_2;
OP_MISC_pct_2 = (sum(IMPACT_OP_2(1:6,:)) + sum(IMPACT_OP_2(15:20,:))) ./ IMPACT_TOT_2;

Energy_TOT_2 = sum(IMPACT_OP_2(7:14,:));

% Direct emission
DE_pct_2 = sum(IMPACT_DE_2) ./ IMPACT_TOT_2;

% Avoided Energy
Avoided_pct_2 = sum(IMPACT_avoided_2) ./ IMPACT_TOT_2;

Impact_source_pct_2 = [Concrete_pct_2; Steel_pct_2; Excavation_pct_2; CON_MISC_pct_2; ...
    Pumping_PERM_pct_2; Pumping_IR_pct_2; Pumping_CHEM_pct_2; Pumping_GBT_pct_2; OP_MISC_pct_2; DE_pct_2; Avoided_pct_2];

%Labeled Output
Labeled_Concrete_pct_2 = horzcat('Concrete',num2cell(Concrete_pct_2));
Labeled_Steel_pct_2 = horzcat('Steel',num2cell(Steel_pct_2));
Labeled_Excavation_pct_2 = horzcat('Excavation',num2cell(Excavation_pct_2));
Labeled_CON_MISC_pct_2 = horzcat('Construction Misc.',num2cell(CON_MISC_pct_2));
Labeled_Pumping_PERM_pct_2 = horzcat('Permeate Pumping',num2cell(Pumping_PERM_pct_2));
Labeled_Pumping_IR_pct_2 = horzcat('Recirculation Pumping',num2cell(Pumping_IR_pct_2));
Labeled_Pumping_CHEM_pct_2 = horzcat('Chemical Pumping',num2cell(Pumping_CHEM_pct_2));
Labeled_Pumping_GBT_pct_2 = horzcat('Gravity Belt Thickeners',num2cell(Pumping_GBT_pct_2));
Labeled_OP_MISC_pct_2 = horzcat('Operation Misc.',num2cell(OP_MISC_pct_2));
Labeled_DE_pct_2 = horzcat('Direct Emissions',num2cell(DE_pct_2));
Labeled_Avoided_pct_2 = horzcat('Avoided Energy',num2cell(Avoided_pct_2));

Labeled_Impact_source_pct_2 = {'CSTR,Submerged with GAC','OZO','GWM','SMG','ACD','EUT', 'CAR','NON-CAR','RSP','ECOTOX'};
Labeled_Impact_source_pct_2=vertcat(Labeled_Impact_source_pct_2,Labeled_Concrete_pct_2,Labeled_Steel_pct_2,Labeled_Excavation_pct_2,...
    Labeled_CON_MISC_pct_2,Labeled_Pumping_PERM_pct_2,Labeled_Pumping_IR_pct_2,Labeled_Pumping_CHEM_pct_2,Labeled_Pumping_GBT_pct_2,...
    Labeled_OP_MISC_pct_2,Labeled_DE_pct_2,Labeled_Avoided_pct_2);
%% Plot (2) Impact by sources 
barExtended(Impact_source_pct_2');
Impact_Categ = {'OZO';'GWM';'SMG';'ACD';...
    'EUT'; 'CAR';'NON-CAR';'RSP';...
    'ECOTOX'};
set(gca,'xticklabel',Impact_Categ);
legend('Concrete', 'Steel', 'Excavation', 'Construction Misc.', ...
    'Permeate Pumping', 'Recirculation Pumping', 'Chemical Pumping', 'Gravity belt thickeners', 'Operation Misc.', 'Direction Emissions', 'Avoided Energy', 'Location', 'NorthEastOutside');




%% System 3: CSTR + Cross-flow (Multi-tube)
disp('     System 3');
%% Inventory Analysis
[INV_CON_3, INV_OP_3, DEmis_3, Power_pct_3, E_input_kWh_3, E_offset_kWh_3, V_treated_3, Output_cost_3] = LCI_anMBR('CSTR', 'none', 'Cross-flow', 'MT', 'PET', DVariable_3, UParameter);

%% Impact Assessment
[IMPACT_CON_3, IMPACT_OP_3, IMPACT_DE_3, IMPACT_avoided_3] = Impact_Assessment (INV_CON_3, INV_OP_3, DEmis_3, E_offset_kWh_3, V_treated_3);

% Total Impact (construction + operation)
CON_TOT_3 = sum(IMPACT_CON_3);
OP_TOT_3 = sum(IMPACT_OP_3); 
IMPACT_TOT_3 = sum(IMPACT_CON_3) + sum(IMPACT_OP_3) + sum(IMPACT_DE_3);

%% Impact by Sources
% Construction phase
Concrete_pct_3 = (IMPACT_CON_3(1,:)) ./ IMPACT_TOT_3;
Steel_pct_3 = (IMPACT_CON_3(2,:) + IMPACT_CON_3(8,:)) ./ IMPACT_TOT_3;
Excavation_pct_3 = (IMPACT_CON_3(18,:)) ./ IMPACT_TOT_3;
CON_MISC_pct_3 = (sum(IMPACT_CON_3(3:7,:)) + sum(IMPACT_CON_3(9:17,:)) + sum(IMPACT_CON_3(19:39,:))) ./ IMPACT_TOT_3;

% Operation phase
Pumping_PERM_pct_3  = (Power_pct_3(1) * sum(IMPACT_OP_3(7:14,:))) ./ IMPACT_TOT_3;
Pumping_IR_pct_3 = (Power_pct_3(2) * sum(IMPACT_OP_3(7:14,:))) ./ IMPACT_TOT_3;
Pumping_R_pct_3 = (Power_pct_3(3) * sum(IMPACT_OP_3(7:14,:))) ./ IMPACT_TOT_3;
Pumping_CHEM_pct_3 = (Power_pct_3(4) * sum(IMPACT_OP_3(7:14,:))) ./ IMPACT_TOT_3;
Pumping_GBT_pct_3 = (Power_pct_3(5) * sum(IMPACT_OP_3(7:14,:))) ./ IMPACT_TOT_3;
OP_MISC_pct_3 = (sum(IMPACT_OP_3(1:6,:)) + sum(IMPACT_OP_3(15:20,:))) ./ IMPACT_TOT_3;

Energy_TOT_3 = sum(IMPACT_OP_3(7:14,:));

% Direct emission
DE_pct_3 = sum(IMPACT_DE_3) ./ IMPACT_TOT_3;

% Avoided Energy
Avoided_pct_3 = sum(IMPACT_avoided_3) ./ IMPACT_TOT_3;

Impact_source_pct_3 = [Concrete_pct_3; Steel_pct_3; Excavation_pct_3; CON_MISC_pct_3; ...
    Pumping_PERM_pct_3; Pumping_IR_pct_3; Pumping_R_pct_3; Pumping_CHEM_pct_3; Pumping_GBT_pct_3; OP_MISC_pct_3; DE_pct_3; Avoided_pct_3];

%Labeled Output
Labeled_Concrete_pct_3 = horzcat('Concrete',num2cell(Concrete_pct_3));
Labeled_Steel_pct_3 = horzcat('Steel',num2cell(Steel_pct_3));
Labeled_Excavation_pct_3 = horzcat('Excavation',num2cell(Excavation_pct_3));
Labeled_CON_MISC_pct_3 = horzcat('Construction Misc.',num2cell(CON_MISC_pct_3));
Labeled_Pumping_PERM_pct_3 = horzcat('Permeate Pumping',num2cell(Pumping_PERM_pct_3));
Labeled_Pumping_IR_pct_3 = horzcat('Recirculation Pumping',num2cell(Pumping_IR_pct_3));
Labeled_Pumping_R_pct_3 = horzcat('Gas Sparging',num2cell(Pumping_R_pct_3));
Labeled_Pumping_CHEM_pct_3 = horzcat('Chemical Pumping',num2cell(Pumping_CHEM_pct_3));
Labeled_Pumping_GBT_pct_3 = horzcat('Gravity Belt Thickeners',num2cell(Pumping_GBT_pct_3));
Labeled_OP_MISC_pct_3 = horzcat('Operation Misc.',num2cell(OP_MISC_pct_3));
Labeled_DE_pct_3 = horzcat('Direct Emissions',num2cell(DE_pct_3));
Labeled_Avoided_pct_3 = horzcat('Avoided Energy',num2cell(Avoided_pct_3));

Labeled_Impact_source_pct_3 = {'CSTR,Cross-flow (Multi-tube)','OZO','GWM','SMG','ACD','EUT', 'CAR','NON-CAR','RSP','ECOTOX'};
Labeled_Impact_source_pct_3=vertcat(Labeled_Impact_source_pct_3,Labeled_Concrete_pct_3,Labeled_Steel_pct_3,Labeled_Excavation_pct_3,...
    Labeled_CON_MISC_pct_3,Labeled_Pumping_PERM_pct_3,Labeled_Pumping_IR_pct_3,Labeled_Pumping_R_pct_3,Labeled_Pumping_CHEM_pct_3,Labeled_Pumping_GBT_pct_3,...
    Labeled_OP_MISC_pct_3,Labeled_DE_pct_3,Labeled_Avoided_pct_3);

%% Plot (3) Impact by sources 
barExtended(Impact_source_pct_3');
Impact_Categ = {'OZO';'GWM';'SMG';'ACD';...
    'EUT'; 'CAR';'NON-CAR';'RSP';...
    'ECOTOX'};
set(gca,'xticklabel',Impact_Categ);
legend('Concrete', 'Steel', 'Excavation', 'Construction Misc.', ...
    'Permeate Pumping', 'Recirculation Pumping', 'Retentate Pumping', 'Chemical Pumping', 'Gravity belt thickeners', 'Operation Misc.', 'Direction Emissions', 'Avoided Energy', 'Location', 'NorthEastOutside');




%% System 4: AF + Crossflow (Multi-tube)
disp('     System 4');
%% LCI
[INV_CON_4, INV_OP_4, DEmis_4, Power_pct_4, E_input_kWh_4, E_offset_kWh_4, V_treated_4, Output_cost_4] = LCI_anMBR('AF', 'none', 'Cross-flow', 'MT', 'PET', DVariable_4, UParameter);

%% LCI --> Impact
[IMPACT_CON_4, IMPACT_OP_4, IMPACT_DE_4, IMPACT_avoided_4] = Impact_Assessment (INV_CON_4, INV_OP_4, DEmis_4, E_offset_kWh_4, V_treated_4);

% Total Impact (construction + operation)
CON_TOT_4 = sum(IMPACT_CON_4);
OP_TOT_4 = sum(IMPACT_OP_4); 
IMPACT_TOT_4 = sum(IMPACT_CON_4) + sum(IMPACT_OP_4) + sum(IMPACT_DE_4);

%% Impact by Sources
% Construction phase
Concrete_pct_4 = (IMPACT_CON_4(1,:)) ./ IMPACT_TOT_4;
Steel_pct_4 = (IMPACT_CON_4(2,:) + IMPACT_CON_4(8,:)) ./ IMPACT_TOT_4;
Excavation_pct_4 = (IMPACT_CON_4(18,:)) ./ IMPACT_TOT_4;
CON_MISC_pct_4 = (sum(IMPACT_CON_4(3:7,:)) + sum(IMPACT_CON_4(9:17,:)) + sum(IMPACT_CON_4(19:39,:))) ./ IMPACT_TOT_4;

% Operation phase
Pumping_LIFT_pct_4  = (Power_pct_4(1) * sum(IMPACT_OP_4(7:14,:))) ./ IMPACT_TOT_4;
Pumping_PERM_pct_4  = (Power_pct_4(2) * sum(IMPACT_OP_4(7:14,:))) ./ IMPACT_TOT_4;
Pumping_IR_pct_4 = (Power_pct_4(3) * sum(IMPACT_OP_4(7:14,:))) ./ IMPACT_TOT_4;
Pumping_R_pct_4 = (Power_pct_4(4) * sum(IMPACT_OP_4(7:14,:))) ./ IMPACT_TOT_4;
Pumping_CHEM_pct_4 = (Power_pct_4(5) * sum(IMPACT_OP_4(7:14,:))) ./ IMPACT_TOT_4;
Pumping_GBT_pct_4 = (Power_pct_4(6) * sum(IMPACT_OP_4(7:14,:))) ./ IMPACT_TOT_4;
OP_MISC_pct_4 = (sum(IMPACT_OP_4(1:6,:)) + sum(IMPACT_OP_4(15:20,:))) ./ IMPACT_TOT_4;

Energy_TOT_4 = sum(IMPACT_OP_4(7:14,:));

% Direct emission
DE_pct_4 = sum(IMPACT_DE_4) ./ IMPACT_TOT_4;

% Avoided Energy
Avoided_pct_4 = sum(IMPACT_avoided_4) ./ IMPACT_TOT_4;

Impact_source_pct_4 = [Concrete_pct_4; Steel_pct_4; Excavation_pct_4; CON_MISC_pct_4; ...
    Pumping_LIFT_pct_4; Pumping_PERM_pct_4; Pumping_IR_pct_4; Pumping_R_pct_4; Pumping_CHEM_pct_4; Pumping_GBT_pct_4; OP_MISC_pct_4; DE_pct_4; Avoided_pct_4];

%Labeled Output
Labeled_Concrete_pct_4 = horzcat('Concrete',num2cell(Concrete_pct_4));
Labeled_Steel_pct_4 = horzcat('Steel',num2cell(Steel_pct_4));
Labeled_Excavation_pct_4 = horzcat('Excavation',num2cell(Excavation_pct_4));
Labeled_CON_MISC_pct_4 = horzcat('Construction Misc.',num2cell(CON_MISC_pct_4));
Labeled_Pumping_LIFT_4 = horzcat('Lift Pumping',num2cell(Pumping_LIFT_pct_4));
Labeled_Pumping_PERM_pct_4 = horzcat('Permeate Pumping',num2cell(Pumping_PERM_pct_4));
Labeled_Pumping_IR_pct_4 = horzcat('Recirculation Pumping',num2cell(Pumping_IR_pct_4));
Labeled_Pumping_R_pct_4 = horzcat('Gas Sparging',num2cell(Pumping_R_pct_4));
Labeled_Pumping_CHEM_pct_4 = horzcat('Chemical Pumping',num2cell(Pumping_CHEM_pct_4));
Labeled_Pumping_GBT_pct_4 = horzcat('Gravity Belt Thickeners',num2cell(Pumping_GBT_pct_4));
Labeled_OP_MISC_pct_4 = horzcat('Operation Misc.',num2cell(OP_MISC_pct_4));
Labeled_DE_pct_4 = horzcat('Direct Emissions',num2cell(DE_pct_4));
Labeled_Avoided_pct_4 = horzcat('Avoided Energy',num2cell(Avoided_pct_4));

Labeled_Impact_source_pct_4 = {'AF,Cross-flow (Multi-tube)','OZO','GWM','SMG','ACD','EUT', 'CAR','NON-CAR','RSP','ECOTOX'};
Labeled_Impact_source_pct_4=vertcat(Labeled_Impact_source_pct_4,Labeled_Concrete_pct_4,Labeled_Steel_pct_4,Labeled_Excavation_pct_4,...
    Labeled_CON_MISC_pct_4,Labeled_Pumping_LIFT_4,Labeled_Pumping_PERM_pct_4,Labeled_Pumping_IR_pct_4,Labeled_Pumping_R_pct_4,Labeled_Pumping_CHEM_pct_4,Labeled_Pumping_GBT_pct_4,...
    Labeled_OP_MISC_pct_4,Labeled_DE_pct_4,Labeled_Avoided_pct_4);


%% Plot (4) Impact by sources 
barExtended(Impact_source_pct_4');
Impact_Categ = {'OZO';'GWM';'SMG';'ACD';...
    'EUT'; 'CAR';'NON-CAR';'RSP';...
    'ECOTOX'};
set(gca,'xticklabel',Impact_Categ);
legend('Concrete', 'Steel', 'Excavation', 'Construction Misc.', ...
    'Lift Pumping', 'Permeate Pumping', 'Recirculation Pumping', 'Retentate Pumping', 'Chemical Pumping', 'Gravity belt thickeners', 'Operation Misc.', 'Direction Emissions', 'Avoided Energy', 'Location', 'NorthEastOutside');



%% System 5: AF + AER + Crossflow (Multi-tube)
disp('     System 5');
%% LCI
[INV_CON_5, INV_OP_5, DEmis_5, Power_pct_5, E_input_kWh_5, E_offset_kWh_5, V_treated_5, Output_cost_5] = LCI_anMBR('AF', '+AER', 'Cross-flow', 'MT', 'PET', DVariable_5, UParameter);

%% LCI --> Impact
[IMPACT_CON_5, IMPACT_OP_5, IMPACT_DE_5, IMPACT_avoided_5] = Impact_Assessment (INV_CON_5, INV_OP_5, DEmis_5, E_offset_kWh_5, V_treated_5);

% Total Impact (construction + operation)
CON_TOT_5 = sum(IMPACT_CON_5);
OP_TOT_5 = sum(IMPACT_OP_5); 
IMPACT_TOT_5 = sum(IMPACT_CON_5) + sum(IMPACT_OP_5) + sum(IMPACT_DE_5);

%% Impact by Sources
% Construction phase
Concrete_pct_5 = (IMPACT_CON_5(1,:)) ./ IMPACT_TOT_5;
Steel_pct_5 = (IMPACT_CON_5(2,:) + IMPACT_CON_5(8,:)) ./ IMPACT_TOT_5;
Excavation_pct_5 = (IMPACT_CON_5(18,:)) ./ IMPACT_TOT_5;
CON_MISC_pct_5 = (sum(IMPACT_CON_5(3:7,:)) + sum(IMPACT_CON_5(9:17,:)) + sum(IMPACT_CON_5(19:39,:))) ./ IMPACT_TOT_5;

% Operation phase
Pumping_LIFT_pct_5  = (Power_pct_5(1) * sum(IMPACT_OP_5(7:14,:))) ./ IMPACT_TOT_5;
Pumping_PERM_pct_5  = (Power_pct_5(2) * sum(IMPACT_OP_5(7:14,:))) ./ IMPACT_TOT_5;
Pumping_IR_pct_5 = (Power_pct_5(3) * sum(IMPACT_OP_5(7:14,:))) ./ IMPACT_TOT_5;
Pumping_R_pct_5 = (Power_pct_5(4) * sum(IMPACT_OP_5(7:14,:))) ./ IMPACT_TOT_5;
Pumping_CHEM_pct_5 = (Power_pct_5(5) * sum(IMPACT_OP_5(7:14,:))) ./ IMPACT_TOT_5;
Pumping_GBT_pct_5 = (Power_pct_5(6) * sum(IMPACT_OP_5(7:14,:))) ./ IMPACT_TOT_5;
OP_MISC_pct_5 = (sum(IMPACT_OP_5(1:6,:)) + sum(IMPACT_OP_5(15:20,:))) ./ IMPACT_TOT_5;

Energy_TOT_5 = sum(IMPACT_OP_5(7:14,:));

% Direct emission
DE_pct_5 = sum(IMPACT_DE_5) ./ IMPACT_TOT_5;

% Avoided Energy
Avoided_pct_5 = sum(IMPACT_avoided_5) ./ IMPACT_TOT_5;

Impact_source_pct_5 = [Concrete_pct_5; Steel_pct_5; Excavation_pct_5; CON_MISC_pct_5; ...
    Pumping_LIFT_pct_5; Pumping_PERM_pct_5; Pumping_IR_pct_5; Pumping_R_pct_5; Pumping_CHEM_pct_5; Pumping_GBT_pct_5; OP_MISC_pct_5; DE_pct_5; Avoided_pct_5];

%Labeled Output
Labeled_Concrete_pct_5 = horzcat('Concrete',num2cell(Concrete_pct_5));
Labeled_Steel_pct_5 = horzcat('Steel',num2cell(Steel_pct_5));
Labeled_Excavation_pct_5 = horzcat('Excavation',num2cell(Excavation_pct_5));
Labeled_CON_MISC_pct_5 = horzcat('Construction Misc.',num2cell(CON_MISC_pct_5));
Labeled_Pumping_LIFT_5 = horzcat('Lift Pumping',num2cell(Pumping_LIFT_pct_5));
Labeled_Pumping_PERM_pct_5 = horzcat('Permeate Pumping',num2cell(Pumping_PERM_pct_5));
Labeled_Pumping_IR_pct_5 = horzcat('Recirculation Pumping',num2cell(Pumping_IR_pct_5));
Labeled_Pumping_R_pct_5 = horzcat('Gas Sparging',num2cell(Pumping_R_pct_5));
Labeled_Pumping_CHEM_pct_5 = horzcat('Chemical Pumping',num2cell(Pumping_CHEM_pct_5));
Labeled_Pumping_GBT_pct_5 = horzcat('Gravity Belt Thickeners',num2cell(Pumping_GBT_pct_5));
Labeled_OP_MISC_pct_5 = horzcat('Operation Misc.',num2cell(OP_MISC_pct_5));
Labeled_DE_pct_5 = horzcat('Direct Emissions',num2cell(DE_pct_5));
Labeled_Avoided_pct_5 = horzcat('Avoided Energy',num2cell(Avoided_pct_5));

Labeled_Impact_source_pct_5 = {'AF + AER + Crossflow (Multi-tube)','OZO','GWM','SMG','ACD','EUT', 'CAR','NON-CAR','RSP','ECOTOX'};
Labeled_Impact_source_pct_5=vertcat(Labeled_Impact_source_pct_5,Labeled_Concrete_pct_5,Labeled_Steel_pct_5,Labeled_Excavation_pct_5,...
    Labeled_CON_MISC_pct_5,Labeled_Pumping_LIFT_5,Labeled_Pumping_PERM_pct_5,Labeled_Pumping_IR_pct_5,Labeled_Pumping_R_pct_5,Labeled_Pumping_CHEM_pct_5,Labeled_Pumping_GBT_pct_5,...
    Labeled_OP_MISC_pct_5,Labeled_DE_pct_5,Labeled_Avoided_pct_5);

%% Plot (5): Impact by sources 
barExtended(Impact_source_pct_5');
Impact_Categ = {'OZO';'GWM';'SMG';'ACD';...
    'EUT'; 'CAR';'NON-CAR';'RSP';...
    'ECOTOX'};
set(gca,'xticklabel',Impact_Categ);
legend('Concrete', 'Steel', 'Excavation', 'Construction Misc.', ...
    'Lift Pumping', 'Permeate Pumping', 'Recirculation Pumping', 'Retentate Pumping', 'Chemical Pumping', 'Gravity belt thickeners', 'Operation Misc.', 'Direction Emissions', 'Avoided Energy', 'Location', 'NorthEastOutside');

%% Plot (6) Total Impact 
figure
% bar([IMPACT_TOT_1 ./IMPACT_TOT_1; IMPACT_TOT_2 ./IMPACT_TOT_1; IMPACT_TOT_3 ./IMPACT_TOT_1; IMPACT_TOT_4 ./IMPACT_TOT_1; IMPACT_TOT_5 ./IMPACT_TOT_1]');
% hold
bar([sum(IMPACT_avoided_1) ./IMPACT_TOT_1; sum(IMPACT_avoided_2) ./IMPACT_TOT_1; sum(IMPACT_avoided_3) ./IMPACT_TOT_1; sum(IMPACT_avoided_4) ./IMPACT_TOT_1; sum(IMPACT_avoided_5) ./IMPACT_TOT_1]');
Impact_Categ = {'OZO';'GWM';'SMG';'ACD';...
    'EUT'; 'CAR';'NON-CAR';'RSP';...
    'ECOTOX'};
set(gca,'xticklabel',Impact_Categ,'fontSize',10);
legend('System 1', 'System 2', 'System 3', 'System 4', 'System 5', 5, 'Location', 'NorthEastOutside')
hold off

%% Plot (7) Construction Inventory comparison
figure
bar([CON_TOT_1 ./CON_TOT_1; CON_TOT_2 ./CON_TOT_1; CON_TOT_3 ./CON_TOT_1; CON_TOT_4 ./CON_TOT_1; CON_TOT_5 ./CON_TOT_1]');
Impact_Categ = {'OZO';'GWM';'SMG';'ACD';...
    'EUT'; 'CAR';'NON-CAR';'RSP';...
    'ECOTOX'};
set(gca,'xticklabel',Impact_Categ,'fontSize',10);
legend('System 1', 'System 2', 'System 3', 'System 4', 'System 5', 5, 'Location', 'NorthEastOutside')
hold off


%% Plot (8) Operational Inventory comparison
figure
bar([OP_TOT_1 ./OP_TOT_1; OP_TOT_2 ./OP_TOT_1; OP_TOT_3 ./OP_TOT_1; OP_TOT_4 ./OP_TOT_1; OP_TOT_5 ./OP_TOT_1]');
Impact_Categ = {'OZO';'GWM';'SMG';'ACD';...
    'EUT'; 'CAR';'NON-CAR';'RSP';...
    'ECOTOX'};
set(gca,'xticklabel',Impact_Categ,'fontSize',10);
legend('System 1', 'System 2', 'System 3', 'System 4', 'System 5', 5, 'Location', 'NorthEastOutside')
hold off

%% Plot (9) Energy consumption 
figure
bar([Energy_TOT_1 ./Energy_TOT_1; Energy_TOT_2 ./Energy_TOT_1; Energy_TOT_3 ./Energy_TOT_1; Energy_TOT_4 ./Energy_TOT_1; Energy_TOT_5 ./Energy_TOT_1]');
Impact_Categ = {'OZO';'GWM';'SMG';'ACD';...
    'EUT'; 'CAR';'NON-CAR';'RSP';...
    'ECOTOX'};
set(gca,'xticklabel',Impact_Categ,'fontSize',10);
legend('System 1', 'System 2', 'System 3', 'System 4', 'System 5', 5, 'Location', 'NorthEastOutside')
hold off

%% Final Labeled Outputs
Labeled_output = vertcat(Labeled_Impact_source_pct,Labeled_Impact_source_pct_2,Labeled_Impact_source_pct_3,Labeled_Impact_source_pct_4,Labeled_Impact_source_pct_5);