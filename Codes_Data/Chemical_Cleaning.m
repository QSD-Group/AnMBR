function [M_NaOCl_kg, Q_NaOCl_weekly, M_CA_kg, Q_CA_weekly] = Chemical_Cleaning (Q_mgd, Year)
% Input:
% Inflent water flow rate, Q_mgd [mgd]
% Operation time, Year [years]

% --> Clean In Place (CIP): Weekly cleaning with 500 mg/L NaOCl and 2000 mg/L citric acid
% --> Clean Out of Place (COP): Biannual cleaning with ? mg/L NaOCl and ? mg/L citric acid

% Sodium Hypochlorite (12.5% solution, 15% by volume)
Dose_NaOCl = 2200; % [gal/yr/mgd] NaOCl Usage Rate
Q_NaOCl_annual = Dose_NaOCl * Q_mgd; % [gal/yr] NaOCl annual flow rate
Q_NaOCl_weekly = Q_NaOCl_annual / 52; % [gal/week] NaOCl weekly flow rate
M_NaOCl_kg = Q_NaOCl_annual * 3.78541 * (12.5/15) * Year; % [kg] Mass of NaClO consumption over N years
% 12.5% by weight = 12.5 g solute/100 mL solution
% 15% by volume = 15 mL solute/100 mL solution
% (12.5 kg/15 L)
% 1 gal = 3.78541 L


% Citric Acid (100% solution, 13.8 lb/gal)
Dose_CA = 600; % [gal/yr/mgd] Citric acid Usage Rate
Q_CA_annual = Dose_CA * Q_mgd; % [gal/yr] Citric acid annual flow rate
Q_CA_weekly = Q_CA_annual / 52; % [gal/week] Citric acid weekly flow rate
M_CA_kg = Q_CA_annual * 13.8 * Year * 0.453592; % [kg] Mass of Citric acid consumption over N years
% 1 lb = 0.453592 kg

end