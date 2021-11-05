function [D_train, L_CSTR, W_train, W_N_trains, W_PB, L_BB, VWC, VSC, VEX] = CSTR_Submerged(Q_mgd, HRT, N_train, Cas_per_tank,HRT_membrane_tank,L_membrane_tank, ia1,HRT_for_GAC)
% Input: 
    % Influent flow rate, Q_cfh [cfh]
    % Number of trains, N_train
    % Number of membrane cassettes per train, N_cassette_pt

% --> CSTR Design 
W_train = 21; % [ft] Width of one train
D_train = 12; % [ft] Depth of one train
Q_cfh = Q_mgd * 133681 / 24; % [ft^3/hr] unit conversion
if ia1==2
    HRT_CSTR = HRT_for_GAC-HRT_membrane_tank;
else
    HRT_CSTR = HRT-HRT_membrane_tank;
end
if HRT_CSTR < 0
    HRT_CSTR=0;
end
L_CSTR = Q_cfh * HRT_CSTR / (N_train * W_train * D_train); % [ft] Length of CSTR

% Width of pump building, W_PB [ft] (based on Hazen & Sawyer data)
if Cas_per_tank >= 0 && Cas_per_tank <= 10 
    W_PB = 27 + 4/12;
elseif Cas_per_tank>= 11 && Cas_per_tank <= 16
    W_PB = 29 + 6/12;
elseif Cas_per_tank >= 17 && Cas_per_tank <= 22
    W_PB = 31 + 8/12;
elseif Cas_per_tank >= 23 && Cas_per_tank <= 28
    W_PB = 35;
elseif Cas_per_tank >= 29 
    W_PB = 38 + 4/12;
else
    W_PB = 0;
end

% Width of blower building, W_BB [ft]
if Cas_per_tank <= 18
    W_BB = 18 + 8/12;  
else
    W_BB = 22;
end

% Length of blower building, L_BB [ft]
if Cas_per_tank <=18
    L_BB = 69 + 6/12;
else
    L_BB = 76 + 8/12;
end

% --> Concrete 
% (assume walls are built on slabs)
% Concrete wall thickness [ft]
if D_train < 12
    t_wall = 1; % Minimum wall thickness = 12 inches
else
    t_wall = 1 + (D_train - 12)/12; % Adding an inch for every foot of depth over 12 ft
end
    
% Concrete slab thickness [ft]
t_slab = t_wall + 2/12; % Slab thickness = wall thickness + 2 inches

% Concrete [ft^3] - Distribution Channel
L_dist = 4.5; % [ft] Width of distribution channel
W_N_trains = (W_train + 2 * t_wall) * N_train - t_wall * (N_train - 1);
VWC_dist = (D_train+2) * t_wall * (2 * W_N_trains + 2 * L_dist); %2 ft freeboard added to depth of train
VSC_dist = W_N_trains * (L_dist + 2 * t_wall)*t_wall+ W_N_trains * (L_dist + 2 * t_wall) * t_slab;

% Concrete [ft^3] - CSTR Trains
VWC_CSTR = (D_train+2) * t_wall * (N_train + 1) * L_CSTR;
VSC_CSTR = (W_N_trains) * L_CSTR*t_wall+(W_N_trains) * t_slab * L_CSTR;

% Concrete [ft^3] - Membrane tanks
VWC_membrane_tank = (D_train+2) * t_wall * (N_train + 1) * L_membrane_tank;
VSC_membrane_tank = (W_N_trains) * L_membrane_tank*t_wall+(W_N_trains) * t_slab * L_membrane_tank;

% Concrete [ft^3] - Effluent Channel 
L_eff = 4.5; % [ft] Width of effluent channel 
VWC_eff = (D_train+2) * t_wall * (2 * W_N_trains +  2 * L_eff);
VSC_eff = W_N_trains * (L_eff + 2 * t_wall)*t_wall+ W_N_trains * (L_eff + 2 * t_wall) * t_slab;

% Pump/Blower Building
VWC_PBB = (D_train+2) * t_wall * (2 * W_N_trains + 2 * W_PB + 2 * W_BB);
VSC_PBB = W_N_trains * (W_PB + t_wall + W_BB)*t_wall+W_N_trains * (W_PB + t_wall + W_BB) * t_slab;

% Wet Well (for mix liquor storage)
D_WW = 12; % [ft] Depth of wet well
W_WW = 8; % [ft] Width of wet well
L_WW = 8; % [ft] Length of wet well
VWC_WW = D_WW * 2*(L_WW * t_wall + (W_WW+ (2 * t_wall)) * t_wall);
VSC_WW = t_slab * (L_WW + 2 * t_wall) * (W_WW + 2 * t_wall)+t_wall *(L_WW + 2 * t_wall) * (W_WW + 2 * t_wall);


% Total Volume of Wall Concrete [f^3]
VWC = VWC_dist + VWC_CSTR + VWC_membrane_tank + VWC_eff + VWC_PBB + VWC_WW;

% Total Volume of Slab Concrete [f^3]
VSC = VSC_dist + VSC_CSTR + VSC_membrane_tank + VSC_eff + VSC_PBB + VSC_WW;


% --> Excavation 
SL = 1.5; % Slope = horizontal/vertical
CA = 3; % [ft] Construction Access

% Volume of excavation of membrane trains [ft^3]
Area_B_train = (L_dist + L_membrane_tank + L_eff + 2 * CA+L_CSTR) * (W_N_trains + 2 * CA); % top area
Area_T_train = (L_dist + L_membrane_tank + L_eff + 2 * CA + D_train * SL+L_CSTR) * (W_N_trains + 2 * CA + D_train * SL); % bottom area
VEX_train = 0.5 * (Area_B_train + Area_T_train) * D_train;

% Volume of excavation of pump/blowr building [ft^3]
Area_B_PB = (W_PB + W_BB + 2 * CA) * (W_N_trains + 2 * CA);
Area_T_P = (W_PB + W_BB + 2 * CA + D_train * SL) * (W_N_trains + 2 * CA + D_train * SL);
VEX_PB = 0.5 * (Area_B_PB + Area_T_P) * D_train;

% Total Volume of Excavation [ft^3]
VEX = VEX_train + VEX_PB;

end