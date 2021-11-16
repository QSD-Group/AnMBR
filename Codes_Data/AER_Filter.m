function [V_m_AER_ft, D_AER_ft, Dia_AER_ft, N_AER, VWC_AER, VSC_AER, VEX] = AER_Filter(Q_mgd, S_SO, X_SO, OL_AER, HL_AER)
% Input:
% Influent flow rate, Q_mgd [mgd]
% Influent readily biodegradable (soluble) substrate concentration, S_SO [mg-BOD5/L or g/m^3]
% Influent slowly biodegradable (particulate) substrate concentration, X_SO [mg-BOD5/L or g/m^3]
% Organic loading rate, OL_AER [kg-BOD5/m^3-day]
% Hydraulic loading rate, HL_AER [m^3/m^2-hr]


%% Aerobic Polishing Filter Design
N_AER = 2;  % Initialize number of AER filters

Removal_ANA = 0.95;

%   Influent Water Parameters
S_SO_AER = S_SO * (1-Removal_ANA); % [mg-BOD5/L] redily biodegradable substrate (Assume 95% removal in ANA filter)
X_SO_AER = X_SO * (1-Removal_ANA); % [mg/L] Slowly biodegradable (particulate) substrate concentration (Assume 95% removal in ANA filter)


% Volume of packing media in each filter, V_m_AER [m^3] (1000 = 'g' to 'kg' unit conversion factor)
Q_cmd = Q_mgd * 3785.41178; % [m^3/day]
V_m_AER = (Q_cmd / N_AER) * (S_SO_AER + X_SO_AER) / OL_AER / 1000;

% X-sectional area of each filter, A_AER [m^2]
Q_cmh = Q_cmd / 24; % [m^3/hour]
A_AER = Q_cmh / N_AER / HL_AER;

% Diameter of each filter, Dia_AER [m]
Dia_AER = (4 * A_AER / pi) ^ 0.5;

% Depth of each ANA filter, D_AER [m]
D_AER = V_m_AER / A_AER;


% Check if more than 1 filter is needed
while Dia_AER > 12 % Maximum diameter assumption [m]
    N_AER = N_AER + 1;
    V_m_AER = (Q_cmd / N_AER) * (S_SO + X_SO) / OL_AER / 1000;
    A_AER = Q_cmh / N_AER / HL_AER;
    Dia_AER = (4 * A_AER / pi) ^ 0.5;
    D_AER = V_m_AER / A_AER;
end

% Unit conversion ('m' to 'ft')
Dia_AER_ft = Dia_AER * 3.28084; % [ft]
D_AER_ft = D_AER * 3.28084; % [ft]
V_m_AER_ft = V_m_AER * 3.28084^3; % [ft^3]

% --> Concrete
% Concrete of anaerobic filter
FB_AER = 3; % [ft] freeboard
t_wall_AER = 8/12; % [ft] wall thickness
t_slab_AER = 1; % [ft] slab thickness

% External wall (wall concrete), VWC_E_AER [ft^3]
VWC_E_AER = t_wall_AER * pi * Dia_AER_ft * (D_AER_ft + FB_AER);

% Floor (slab concrete), VSC_F_AER [ft^3]
VSC_F_AER = (pi/4) * Dia_AER_ft^2+t_slab_AER * (pi/4) * Dia_AER_ft^2;

VWC_AER = N_AER * VWC_E_AER; % [ft^3] Volume of wall concrete
VSC_AER = N_AER * VSC_F_AER; % [ft^3] Volume of slab concrete


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
