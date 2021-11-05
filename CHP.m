function [P_offset,heat] = CHP (Q_CH4)
% Input Parameters:
% Q_biogas: gas flow [m^3-CH4/day]

% Combined heat and power (CHP)
EC_CH4 = 9.8; % [kWh/m3 CH4] Energetic content of methane
kJ_CH4 = 35846; % [kJ/m3 CH4]

%Efficiencies from Woer et al., 2012 (Range given, took midpoint)

% if strcmp(step_I,'IC')
%     eta_power = 0.36;
%     eta_heat = 0.42;
% elseif strcmp(step_I,'CG')
%     eta_power = 0.315;
%     eta_heat = 0.41;
% elseif strcmp(step_I,'micro')
%     eta_power = 0.28;
%     eta_heat = 0.335;
% else
%     eta_power = 0.405;
%     eta_heat = 0.35;
% end

eta_power=0.27; %Pretel et al., 2013
eta_heat = 0.335; %Pretel et al., 2013

P_offset = Q_CH4 * EC_CH4 * eta_power / 24; % [kW] Electricity offset power
heat = Q_CH4*kJ_CH4*eta_heat; %[kJ/d]

end