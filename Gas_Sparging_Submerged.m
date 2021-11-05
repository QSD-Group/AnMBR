function [Q_gas_cfm, P_input_blower, M_SS_gas, OD_gh, OD_gsm] = Gas_Sparging_Submerged(N_train, N_cassette_pt, A_cassette, SGD, L_train, W_PB, L_BB, freq)
% Input:
    % Number of trains, N_train
    % Number cassettes per train, N_cassette_pt
    % Surface area of each cassette, A_cassette [ft^3]
    % Specific gas demand, SGD [Nm^3 gas/m^2 membrane area-h]
    % Length of each train, L_train [ft]
    % Width of pump building, W_PB [ft]
    % Length of Blower building, L_BB [ft]
    % Sparging frequency, freq

% NB = ceil(N_train / 2); % # of blowers in duty (assume 1 blowers per 2 trains) !!Not actually used in LCI code!!

% Air Header
L_gh = L_train; % [ft] single length of air header 

% Air Supply Manifold
L_gsm = (21 + 2)*N_train - (N_train - 1)*1 + W_PB + L_BB; % [ft] length of air supply manifold

% Gas Requirement 
Q_gas_cfm_pt = (SGD * 3.2808399/60) * N_cassette_pt * A_cassette; % [ft^3/min] gas requirement per train
Q_gas_cfm = Q_gas_cfm_pt * N_train; % [ft^3/min] total gas requirement 
Q_gas_cfs_pt = Q_gas_cfm_pt * 60; % [ft^3/s] gas requirement per train
Q_gas_cfs = Q_gas_cfs_pt * N_train; % [ft^3/s] total gas requirement 


% Air Header
VEL_gh = 70; % [ft/s] air velocity in air headers (Median of values, Metcalf & Eddy, 2013, Table 5-29)
[OD_gh, t_gh, ID_gh] = pipe(Q_gas_cfs_pt, VEL_gh);

% Air Supply Manifold
VEL_gsm = 70; % [ft/s] air velocity in air supply manifold (Metcalf & Eddy, 2013, Table 5-29)
[OD_gsm, t_gsm, ID_gsm] = pipe(Q_gas_cfs, VEL_gsm);

% Pipe material (assume stainless steel, density = 0.29 lbs/in^3)
V_gh = N_train * pi/4 * ((OD_gh)^2 - (ID_gh)^2) * (L_gh * 12);
V_gsm = pi/4 * ((OD_gsm)^2 - (ID_gsm)^2) * (L_gsm * 12);
M_SS_gh = 0.29 * V_gh * 0.453592; % [kg] mass of stainless steel 
M_SS_gsm = 0.29 * V_gsm * 0.453592; % [kg] mass of stainless steel 
M_SS_blw = 1000; % [kg] mass of stainless steel for blower, assuming blower weighs 1000 kg, half is composed of SS and 15-year useful life (so doubled to 1000)
M_SS_gas = M_SS_gh + M_SS_gsm + M_SS_blw; % [kg] mass of stainless steel 
                                       
TDH_blower_psig = 6; % [psig] estimated TDH (estimation based on Hazen&Sawyer spreadsheet)
Eff_blower = 0.7; % blower efficiency
Eff_motor_BL = 0.7; % motor efficiency
BHP_blower_tot = (Q_gas_cfm * 0.23) * (((14.7 + TDH_blower_psig)/14.7)^0.283 - 1.0) / Eff_blower; 
P_input_blower = BHP_blower_tot * 0.746 / Eff_motor_BL * freq; % [kW] Power input to motor 

end