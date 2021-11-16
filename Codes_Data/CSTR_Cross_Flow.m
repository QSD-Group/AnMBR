function [D_train, L_CSTR, W_train, W_N_trains, W_PB, L_BB, VWC, VSC, VEX] = CSTR_Cross_Flow(Q_mgd, HRT, N_train,L_membrane_tank)
% Input:
% Influent flow rate, Q_mgd [mgd]
% Hydraulic retention time, HRT [hr]
% Number of trains, N_train


% --> CSTR Design
W_train = 21; % [ft] Width of one train
D_train = 12; % [ft] Depth of one train
Q_cfh = Q_mgd * 133681 / 24; % [ft^3/hr] unit conversion
L_CSTR = Q_cfh * HRT / (N_train * W_train * D_train); % [ft] Length of CSTR

% Width of pump building, W_PB [ft] (based on Hazen & Sawyer data)
N_eq = L_CSTR / ((1 + 8/12) + (3 + 4/12)); % based on CSTR with submerged membrane
if N_eq >= 0 && N_eq <= 10
    W_PB = 27 + 4/12;
elseif N_eq>= 11 && N_eq <= 16
    W_PB = 29 + 6/12;
elseif N_eq >= 17 && N_eq <= 22
    W_PB = 31 + 8/12;
elseif N_eq >= 23 && N_eq <= 28
    W_PB = 35;
elseif N_eq >= 29
    W_PB = 38 + 4/12;
else
    W_PB = 0;
end

% Width of blower building, W_BB [ft]
if N_eq <= 18
    W_BB = L_membrane_tank;  %No blower building required for X-flow units, assume all membranes are indoors
else
    W_BB = L_membrane_tank;
end

% Length of blower building, L_BB [ft]
if N_eq <=18
    L_BB = 0;  %No blower building required for X-flow units
else
    L_BB = 0;
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

% Concrete [ft^3] - Effluent Channel
L_eff = 4.5; % [ft] Width of effluent channel
VWC_eff = (D_train+2) * t_wall * (2 * W_N_trains +  2 * L_eff);
VSC_eff = W_N_trains * (L_eff + 2 * t_wall)*t_wall+ W_N_trains * (L_eff + 2 * t_wall) * t_slab;

% Pump/Blower Building
VWC_PBB = (D_train+2) * t_wall * (2 * W_N_trains + 2 * W_PB + 2 * W_BB);
VSC_PBB = W_N_trains * (W_PB + t_wall + W_BB)*t_wall+W_N_trains * (W_PB + t_wall + W_BB) * t_slab;

% Total Volume of Wall Concrete [f^3]
VWC = VWC_dist + VWC_CSTR + VWC_eff + VWC_PBB;

% Total Volume of Slab Concrete [f^3]
VSC = VSC_dist + VSC_CSTR + VSC_eff + VSC_PBB;

% --> Excavation
SL = 1.5; % Slope = horizontal/vertical
CA = 3; % [ft] Construction Access

% Volume of excavation of membrane trains [ft^3]
Area_B_train = (L_dist + L_CSTR + L_eff + 2 * CA) * (W_N_trains + 2 * CA); % top area
Area_T_train = (L_dist + L_CSTR + L_eff + 2 * CA + D_train * SL) * (W_N_trains + 2 * CA + D_train * SL); % bottom area
VEX_train = 0.5 * (Area_B_train + Area_T_train) * D_train;

% Volume of excavation of pump/blowr building [ft^3]
Area_B_PB = (W_PB + W_BB + 2 * CA) * (W_N_trains + 2 * CA);
Area_T_P = (W_PB + W_BB + 2 * CA + D_train * SL) * (W_N_trains + 2 * CA + D_train * SL);
VEX_PB = 0.5 * (Area_B_PB + Area_T_P) * D_train;

% Total Volume of Excavation [ft^3]
VEX = VEX_train + VEX_PB;

end