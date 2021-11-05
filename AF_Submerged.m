function [V_m_AF_ft, D_AF_ft, Dia_AF_ft, N_AF, VWC, VSC, VEX, W_PB, L_BB,W_N_trains] = AF_Submerged(Q_mgd, S_SO, X_SO, OL_AF, HL_AF, R_AF,Cas_per_tank,N_train,L_membrane_tank)
% Input:
% Influent flow rate, Q_mgd [mgd]
% Influent readily biodegradable (soluble) substrate concentration, S_SO [mg-BOD5/L or g/m^3]
% Influent slowly biodegradable (particulate) substrate concentration, X_SO [mg-BOD5/L or g/m^3]
% Organic loading rate, OL_AF [kg-BOD5/m^3-day]
% Hydraulic loading rate, HL_AF [m^3/m^2-hr]
% Recirculation ratio, R_AF


%% Anaerobic Filter Design
N_AF = 2;  % Initialize number of ANA filters
W_train = 21; % [ft] Width of one train

% Volume of packing media in each filter, V_m_AF [m^3] (1000 = 'g' to 'kg' unit conversion factor)
Q_cmd = Q_mgd * 3785.41178; % [m^3/day]
V_m_AF = (Q_cmd / N_AF) * (S_SO + X_SO) / OL_AF / 1000;

% X-sectional area of each filter, A_AF [m^2]
Q_cmh = Q_cmd / 24; % [m^3/hour]
A_AF = Q_cmh * (1 + R_AF) / N_AF / HL_AF;

% Diameter of each filter, Dia_AF [m]
Dia_AF = (4 * A_AF / pi) ^ 0.5;

% Depth of each ANA filter, D_AF [m]
D_AF = V_m_AF / A_AF;

while D_AF > 6 % Maximum depth assumption [m]
    R_AF = R_AF + 0.1;
    A_AF = Q_cmh * (1 + R_AF) / N_AF / HL_AF;
    Dia_AF = (4 * A_AF / pi) ^ 0.5;
    % Check if more than 1 filter is needed
    while Dia_AF > 12 % Maximum diameter assumption [m]
        N_AF = N_AF + 1;
        A_AF = Q_cmh * (1 + R_AF) / N_AF / HL_AF;
        Dia_AF = (4 * A_AF / pi) ^ 0.5;
    end
    V_m_AF = (Q_cmd / N_AF) * (S_SO + X_SO) / OL_AF / 1000;
    D_AF = V_m_AF / A_AF;
end

% Unit conversion ('m' to 'ft')
Dia_AF_ft = Dia_AF * 3.28084; % [ft]
D_AF_ft = D_AF * 3.28084; % [ft]
V_m_AF_ft = V_m_AF * 3.28084^3; % [ft^3]

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
% Concrete of anaerobic filter
FB_AF = 3; % [ft] freeboard
t_wall_AF = 1 + (D_AF - 12)/12; % [ft] wall thickness
t_slab_AF = t_wall_AF + 2/12; % [ft] slab thickness

W_N_trains = (W_train + 2 * t_wall_AF) * N_train - t_wall_AF * (N_train - 1);

% External wall (wall concrete), VWC_E_AF [ft^3]
VWC_E_AF = t_wall_AF * pi * Dia_AF_ft * (D_AF_ft + FB_AF);

% Floor (slab concrete), VSC_F_AF [ft^3]
VSC_F_AF =(pi/4) * Dia_AF_ft^2+ t_slab_AF * (pi/4) * Dia_AF_ft^2;

VWC_AF = N_AF * VWC_E_AF; % [ft^3] Volume of wall concrete
VSC_AF = N_AF * VSC_F_AF; % [ft^3] Volume of slab concrete

% Concrete [ft^3] - Membrane tanks
VWC_membrane_tank = (D_AF+2) * t_wall_AF * (N_train + 1) * L_membrane_tank;
VSC_membrane_tank = (W_N_trains) * L_membrane_tank*t_wall_AF+(W_N_trains) * t_slab_AF * L_membrane_tank;

% Pump/Blower Building
VWC_PBB = (D_AF+2) * t_wall_AF * (2 * W_N_trains + 2 * W_PB + 2 * W_BB);
VSC_PBB = W_N_trains * (W_PB + t_wall_AF + W_BB)*t_wall_AF+W_N_trains * (W_PB + t_wall_AF + W_BB) * t_slab_AF;

% Wet Well (for mix liquor storage)
D_WW = 12; % [ft] Depth of wet well
W_WW = 8; % [ft] Width of wet well
L_WW = 8; % [ft] Length of wet well
VWC_WW = D_WW * 2*(L_WW * t_wall_AF + (W_WW+ (2 * t_wall_AF)) * t_wall_AF);
VSC_WW = t_slab_AF * (L_WW + 2 * t_wall_AF) * (W_WW + 2 * t_wall_AF)+t_wall_AF *(L_WW + 2 * t_wall_AF) * (W_WW + 2 * t_wall_AF);

% Total Volume of Wall Concrete [f^3]
VWC = VWC_AF + VWC_membrane_tank + VWC_PBB + VWC_WW;

% Total Volume of Slab Concrete [f^3]
VSC = VSC_AF + VSC_membrane_tank + VSC_PBB + VSC_WW;

% --> Excavation
SL = 1.5; % Slope = horizontal/vertical
CA = 3; % [ft] Construction Access

% Volume of excavation of membrane trains [ft^3]
Area_B_train = (L_membrane_tank + 2 * CA) * (W_N_trains + 2 * CA); % top area
Area_T_train = (L_membrane_tank + 2 * CA + D_AF * SL) * (W_N_trains + 2 * CA + D_AF * SL); % bottom area
VEX_train = 0.5 * (Area_B_train + Area_T_train) * D_AF;

% Volume of excavation of pump/blowr building [ft^3]
Area_B_PB = (W_PB + W_BB + 2 * CA) * (W_N_trains + 2 * CA);
Area_T_P = (W_PB + W_BB + 2 * CA + D_AF * SL) * (W_N_trains + 2 * CA + D_AF * SL);
VEX_PB = 0.5 * (Area_B_PB + Area_T_P) * D_AF;

% Total Volume of Excavation [ft^3]
VEX = VEX_train + VEX_PB;


end
