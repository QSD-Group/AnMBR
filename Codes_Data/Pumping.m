function [P_input, M_SS, OD_s, OD_d] = Pumping(Q_mgd, NP, H_ts, L_s, L_d, H_p)
% Input Parameters:
% Total pumping flow rate, Q [gmd]
% Number of pumps, NP
% Total static head, H_ts [ft]
% Length of pipes on suction side, L_s [ft]
% Length of pipes on discharge side, L_d [ft]
% Pressure head, H_p [ft]


VEL = 3; % [ft/s] water velocity in pipe
C = 110; % Hazen- Williams coefficient (for S.S. steel)
Q_cfs = Q_mgd * 1.54722865; % [ft^3/s] 
Q_gpm = Q_mgd * 694.444444; % [million gallon/day]


% --> Suction side 
% Suction pipe (permeate headers)
[OD_s, t_s, ID_s] = pipe(Q_cfs/NP, VEL); % [in, in, in]

% Suction Friction Head, H_sf [ft]     
H_sf = 3.02 * L_s * VEL ^1.85 * C^(-1.85) * (ID_s/12)^(-1.17);   


% --> Discharge side 
% Discharge pipe (Permeate collector)
[OD_d, t_d, ID_d] = pipe(Q_cfs, VEL); % [in, in, in]

% Discharge Friction Head, H_df [ft]     
H_df = 3.02 * L_d * VEL ^1.85 * C ^(-1.85) * (ID_d/12)^(-1.17);  


% --> Pipe material (assume stainless steel, density = 0.29 lbs/in^3)  
V_s = NP * pi/4 * ((OD_s)^2 - (ID_s)^2) * (L_s * 12);
V_d = pi/4 * ((OD_d)^2 - (ID_d)^2) * (L_d * 12);
M_SS_pipe = 0.29 * (V_s + V_d) * 0.453592; % [kg] mass of stainless steel (all pipes)
M_SS_pump_ea = 725 * 0.5;  % [kg] mass of stainless steel (each pump) (for a 300 - 1000 gpm pump --> reference http://www.godwinpumps.com/images/uploads/ProductCatalog_Nov_2011_spread2.pdf)
                            % assume 50% of the product weight is stainless steel 
M_SS_pump = M_SS_pump_ea * NP*2;  % [kg] mass of stainless steel (all pumps), assumes 15-year useful life (i.e., need two over a 30-year life for the plant)
M_SS = M_SS_pipe + M_SS_pump; % [kg] total mass of stainless steel required for pumping


% --> Electrical Energy
TDH = H_p + H_ts + H_sf + H_df; % [ft] Total Dynamic Head 
Eff_pump = 0.8; % assume pump efficiency is 80%
Eff_motor = 0.7; % assume motor efficiency is 70% 
BHP = (TDH * Q_gpm) / 3960 / Eff_pump; 
P_input = BHP * 0.746 / Eff_motor; % [Kw] Power input to motor 

end