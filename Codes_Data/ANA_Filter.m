function [V_m_AF_ft, D_AF_ft, Dia_AF_ft, N_AF, VWC_AF, VSC_AF, VEX] = ANA_Filter(Q_mgd, S_SO, X_SO, OL_AF, HL_AF, R_AF)
% Input: 
    % Influent flow rate, Q_mgd [mgd]
    % Influent readily biodegradable (soluble) substrate concentration, S_SO [mg-BOD5/L or g/m^3] 
    % Influent slowly biodegradable (particulate) substrate concentration, X_SO [mg-BOD5/L or g/m^3] 
    % Organic loading rate, OL_AF [kg-BOD5/m^3-day]
    % Hydraulic loading rate, HL_AF [m^3/m^2-hr]
    % Recirculation ratio, R_AF


%% Anaerobic Filter Design
N_AF = 2;  % Initialize number of ANA filters

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

% --> Concrete 
% Concrete of anaerobic filter 
FB_AF = 3; % [ft] freeboard
t_wall_AF = 6/12; % [ft] wall thickness
t_slab_AF = 8/12; % [ft] slab thickness

% External wall (wall concrete), VWC_E_AF [ft^3]
VWC_E_AF = t_wall_AF * pi * Dia_AF_ft * (D_AF_ft + FB_AF);

% Floor (slab concrete), VSC_F_AF [ft^3]
VSC_F_AF =(pi/4) * Dia_AF_ft^2+ t_slab_AF * (pi/4) * Dia_AF_ft^2;

VWC_AF = N_AF * VWC_E_AF; % [ft^3] Volume of wall concrete
VSC_AF = N_AF * VSC_F_AF; % [ft^3] Volume of slab concrete


% --> Excavation 
SL = 1.5; % Slope = horizontal/vertical
CA = 3; % [ft] Construction Access

%   Excavation of Pump Building
PBL = 50; % [ft] Pump Building Length
PBW = 30; % [ft] Pump Building Width
PBD = 10; % [ft] Pump Building Depth
Area_B_P = (PBL + 2 * CA) * (PBW + 2 * CA); % [ft^2] Bottom Area of frustum  
Area_T_P = (PBL + 2 * CA + PBW * SL) * (PBW + 2 * CA + PBD * SL); % [ft^2] Top Area of frustum  
VEX_PB = 0.5 * (Area_B_P + Area_T_P) * PBD; % [ft^2] Volume of excavaion of Pump Building

VEX = VEX_PB; % [ft^3] Total volume of excavation 


end
