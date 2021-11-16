function [P_input_sludge, M_SS_GBT, M_COD_WAS] = Sludge_Handling (Q_CH4_1, Q_WAS_mgd, M_COD_removed) 
% Input:
% Methane production rate, Q_CH4 [m^3-CH4/day]
% Sludge wastage rate, Q_WAS_mgd [mgd]
% Mass of COD removed, M_COD_removed [kg-COD/day]

% 1) Dewatering design steps:
% Assumptions:
% Jeremy's dissertation, Chapter 5
% Gravity belt thickeners (GBTs) were designed to operate at the same 
% frequency as the centrifuges (22.7 hours a day, 7 days a week). It is assumed that no 
% solids storage exists immediately upstream or downstream of the GBTs, requiring GBT 
% operation whenever secondary solids are wasted. GBTs were designed such that the 
% maximum hydraulic loading is no more than 150 gal·m-1·min-1 with all units in service and 
% no more than 200 gal·m-1·min-1 with one unit out of service. 

% 1a) Determine wastage rate, Q_WAS [gal/min]
% 2b) Determine the total width of GBTs, W_GBT_tot = Q_WAS / 150 gal·m-1·min-1
W_GBT_tot = Q_WAS_mgd * 10^6 / (24 * 60) / 150; % [m]
% 3c) Determine the number of GBTs needed, N_GBT = W_GBT_tot / 3 m 
N_GBT = ceil(W_GBT_tot / 3); %
% 4e) Energy consumption, from Jeremy's spreadsheet, find relations from figures
%WAS Pumping
[P_input_WAS, ~, ~, ~] = Pumping(Q_WAS_mgd, 1, 0, 50, 50, 0); %Q_mgd, NP, H_ts, L_s, L_d, H_p

% kWh/yr = 422831.7854 * MGD^0.9248
P_input_GBT = (422831.7854 * Q_WAS_mgd^0.9248) / (24*365); % [kW] This equation is from Jeremy's dissertation 
P_input_sludge = P_input_GBT+P_input_WAS;
% 5f) Material consumption, assume stainless steel
M_SS_GBT_ea = 0; % [kg] Mass of stainless steel per GBT; Changed by Brian and Jeremy
    %(mass based on http://haibar.en.alibaba.com/product/703897144-212173963/Gravity_Belt_Thickener_for_Filter_Press_and_Decanter_Centrifuge_HBT_1500.html)
M_SS_GBT = M_SS_GBT_ea * N_GBT; % [kg] Total mass of stainless steel 

% 2) Landfilling design steps:
% 2a) Calculate M_COD_biogas
% CH4 + 2O2 —> CO2 + 2H2O
% 1 mol-CH4 ~ 2 mol-O2 ~ 64 g-O2
% knowing that 1 mol-CH4 = 22.4 L-CH4
% --> Determine M_COD_biogas = Q_CH4 * (64 g-COD/22.4 L-CH4);
M_COD_biogas = Q_CH4_1 * (0.064 / 0.0224); % [kg/day]
% 2b) Determine mass of COD wasted, M_COD_WAS
M_COD_WAS = M_COD_removed - M_COD_biogas; % [kg/day]
M_COD_WAS = M_COD_WAS/1.5/0.22; %1.5 is the ratio of COD/VSS for medium strength municipal strength WW (Henze, 2008, Biological WW Treatment)
    %0.22 assumes the use of belt filter presses with negligible energy consumption to get 22% solids.

end