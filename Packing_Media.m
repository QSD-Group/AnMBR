function [M_LDPE_kg, M_HDPE_kg] = Packing_Media (V_m) 
% Input:
% Volume of packing media, V_m [m^3]

% Assume packing media is 50% LDPE and 50% HDPE
Ro_LDPE = 925; % [kg/m^3] Density of LDPE
Ro_HDPE = 950; % [kg/m^3] Density of HDPE
    % LDPE density = 0.910-0.940 g/cm3 or 910-940 kg/m^3
    % HDPE density = 0.93-0.97 g/cm3 or 930-970 kg/m^3
Frac_void = 0.9; % Void Fraction, usually 85% - 95% for plastic packing media

M_LDPE_kg = 0.5 * (1-Frac_void) * Ro_LDPE * V_m; % [kg] Mass of LDPE
M_HDPE_kg = 0.5 * (1-Frac_void) * Ro_HDPE * V_m; % [kg] Mass of HDPE
end