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

IRR = 1;


%% Uncertainty parameters
Q_mgd = 20; % [mgd]

S_SO = 300; % [mg-COD/L] or [g/m^3] 

X_SO = 100; % [mg-COD/L] or [g/m^3] 

UParameter = [Q_mgd, S_SO, X_SO];

%% Possible Systems
% Reactor type: step_A = 'CSTR' or 'AF'
choice_A = {'CSTR', 'AF'};
    % If there's polishing aerobic filter: step_A1 = '+AER'
    % If there's GAC in reactor: step_A1 = '+GAC'
    % Otherwise, step_A1 = 'none'
choice_A1 = {'+AER', '+GAC', 'none'};
    
% Configuration: step_B = 'Submerged' or 'Cross-flow'
choice_B = {'Submerged', 'Cross-flow'};

% Membrane type: step_C = 'HF' or 'FS' or 'MT' 
choice_C = {'HF', 'FS', 'MT'};

% Membrane material: step_D = 'PET' or 'PTFE'
choice_D = {'PET', 'PTFE'};

% Chemical cleaning: no choice
% Physical cleaning: 
% Soluble methane management: no choice
% Methane processing: no choice

% Initialize index
% Output A - REACTOR TYPE
i_CSTR = 1;
i_AF = 1;
i_Submerged = 1;
i_Xflow = 1;
i_HF = 1;
i_FS = 1;
i_MT = 1;
i_PET = 1; 
i_PTFE = 1;

for A = 1:2
    for A1 = 1:3
        for B = 1:2
            for C = 1:3
                for D = 1:2
% Assign values to each design choice                    
step_A = choice_A(A);
step_A1 = choice_A1(A1);
step_B = choice_B(B);
step_C = choice_C(C);
step_D = choice_D(D);

% Uncompatible design decisions
if strcmp(step_A, 'CSTR') && strcmp(step_A1, '+AER')
    continue
elseif strcmp(step_A, 'AF') && strcmp(step_B, 'Submerged')
    continue
elseif strcmp(step_B, 'Submerged') && strcmp(step_C, 'MT')
    continue
elseif strcmp(step_B, 'Cross-flow') && strcmp(step_C, 'HF')    
    continue
end

% Assign values to decision variables
if strcmp(step_A, 'CSTR') && strcmp(step_B, 'Submerged') && strcmp(step_A1, 'none')
    DVariable = [J, TMP_high, N_train, IRR];
elseif strcmp(step_A, 'CSTR') && strcmp(step_B, 'Submerged') && strcmp(step_A1, '+GAC')
    DVariable = [J, TMP_low, N_train, IRR];
elseif strcmp(step_A, 'CSTR') && strcmp(step_B, 'Cross-flow') && strcmp(step_A1, 'none')
    DVariable = [J, TMP_high, N_train, IRR, HRT];
elseif strcmp(step_A, 'CSTR') && strcmp(step_B, 'Cross-flow') && strcmp(step_A1, '+GAC')
    DVariable = [J, TMP_low, N_train, IRR, HRT];    
elseif strcmp(step_A, 'AF') && strcmp(step_A1, 'none') 
    DVariable = [J, TMP_high, OL_AF, HL_AF, IRR];
elseif strcmp(step_A, 'AF') && strcmp(step_A1, '+GAC')
    DVariable = [J, TMP_low, OL_AF, HL_AF, IRR];
elseif strcmp(step_A, 'AF') && strcmp(step_A1, '+AER')
    DVariable = [J, TMP_low, OL_AF, HL_AF, IRR, OL_AER, HL_AER];
end
           
%% Inventory Analysis
[INV_CON, INV_OP, DEmis, Power_pct, E_input_kWh, E_offset_kWh, V_treated, Output_cost] = LCI_anMBR(step_A, step_A1, step_B, step_C, step_D, DVariable, UParameter);

%% Impact Assessment
[IMPACT_CON, IMPACT_OP, IMPACT_DE, IMPACT_avoided] = Impact_Assessment (INV_CON, INV_OP, DEmis, E_offset_kWh, V_treated);

% Total Impact (construction + operation)
IMPACT_TOT = sum(IMPACT_CON) + sum(IMPACT_OP) + sum(IMPACT_DE);


%% Output (only output global warming potential)
% Output A - REACTOR TYPE
if strcmp(step_A, 'CSTR')
    A_CSTR(i_CSTR) = IMPACT_TOT(2);
    i_CSTR = i_CSTR + 1;
elseif strcmp(step_A, 'AF')
    A_AF(i_AF) = IMPACT_TOT(2);
    i_AF = i_AF + 1;
end

% Output B - CONFIGURATION	
if strcmp(step_B, 'Submerged')
    B_Submerged(i_Submerged) = IMPACT_TOT(2);
    i_Submerged = i_Submerged + 1;
elseif strcmp(step_B, 'Cross-flow') 
    B_Xflow(i_Xflow) = IMPACT_TOT(2);
    i_Xflow = i_Xflow + 1;
end

% Output C - MEMBRANE TYPE	
if strcmp(step_C, 'HF')
    C_HF(i_HF) = IMPACT_TOT(2);
    i_HF = i_HF + 1;
elseif strcmp(step_C, 'FS')
    C_FS(i_FS) = IMPACT_TOT(2);
    i_FS = i_FS + 1;
elseif strcmp(step_C, 'MT')
    C_MT(i_MT) = IMPACT_TOT(2);
    i_MT = i_MT + 1;
end

% Output D - MEMBRANE MATERIAL
if strcmp(step_D, 'PET')
    D_PET(i_PET) = IMPACT_TOT(2);
    i_PET = i_PET + 1; 
elseif strcmp(step_D, 'PTFE') 
    D_PTFE(i_PTFE) = IMPACT_TOT(2);
    i_PTFE = i_PTFE + 1;
end
% Output E - OPERATION MODE	
% OutputF - PHYSICAL CLEANING	
% Output G - CHEMICAL CLEANING 	
% Output H - SOLUBLE METHANE MANAGEMENT	
% Output I - METHANE PROCESSING
       
                end
            end
        end
    end
end

%% Write to Excel
xlswrite('A_CSTR.xls', A_CSTR', 'Sheet 1', 'A1')
xlswrite('A_AF.xls', A_AF', 'Sheet 1', 'A1')
xlswrite('B_Submerged.xls', B_Submerged', 'Sheet 1', 'A1')
xlswrite('B_Xflow.xls', B_Xflow', 'Sheet 1', 'A1')
xlswrite('C_HF.xls', C_HF', 'Sheet 1', 'A1')
xlswrite('C_FS.xls', C_FS', 'Sheet 1', 'A1')
xlswrite('C_MT.xls', C_MT', 'Sheet 1', 'A1')
xlswrite('D_PET.xls', D_PET', 'Sheet 1', 'A1')
xlswrite('D_PTFE.xls', D_PTFE', 'Sheet 1', 'A1')


