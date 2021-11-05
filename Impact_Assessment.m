function [traci_con, CED_con, traci_op, CED_op, traci_avoid, CED_avoid, IMPACT_DE] = Impact_Assessment (INV_CON, INV_OP, DEmis, E_offset_kWh, V_treated)
% function [traci_con, CML_con,CED_con, traci_op, CML_op,CED_op, traci_avoid, CML_avoid,CED_avoid, IMPACT_DE] = Impact_Assessment (INV_CON, INV_OP, DEmis, E_offset_kWh, V_treated)

% Input:
% Construction Inventories, INV_CON =
% [1-Concrete; 2-Reinforcing steel; 3-Tap water; 4-Aluminum; 5-Limestone;
% 6-Chromium Steel; 7-Flat glass; 8-Copper; 9-Sythetic rubber; 10-Rock wool;
% 11-Bitumen; 12-LDPE; 13-HDPE; 14-Excavation; 15-Operation; 16-Transport;
% 17-Extrusion; 18-Transport; 19-Electricity; 20-Organic chemicals;
% 21-Inorganic chemicals; 22-PET; 23-PTFE; 24-GAC; 25-CHP];
% Operational Inventories, INV_OP =
% [1-Acetic acid; 2-Methanol; 3-Iron; 4-Sodium hypochlorite;
% 5-Electricity; 6-Chlorine; 7-Methyl methcrylate; 8-Disposal; 9-Citric Acid]
% Direction emissions, DEmis (8x1 matrix)
% Output:

% Ouput:
% Impacts from construction phase, IMPACT_CON (normalized to per m^3 wastewater treated)
% Impacts from operational phase, IMPACT_OP (normalized to per m^3 wastewater treated)
% Impacts from direct emissions, IMPACT_DE (normalized to per m^3 wastewater treated)
% Avoided impacts from energy offset, IMPACT_avoided (normalized to per m^3 wastewater treated)

load simapro

traci_con=traci(1:38,:);
traci_op=traci([39:45 22 47 23:29],:);
traci_avoid=traci(23:29,:);
IMPACT_DE = DE;

CED_con=CED(1:38,:);
CED_op=CED([39:45 22 47 23:29],:);
CED_avoid=CED(23:29,:);
% CML_con=CML(1:38,:);
% CML_op=CML([39:45 22 47 23:29],:);
% CML_avoid=CML(23:29,:);

%% Construction Phase
UP_CON = zeros(38,1); % Initialize Ecoinvent Unit Process matrix

UP_CON(2) = INV_CON(1); % Concrete, normal {CH}| production | Alloc Def, U
UP_CON(3) = INV_CON(2); % Reinforcing steel {RER}| production | Alloc Def, U
UP_CON(4) = INV_CON(3); % Tap water, at user/RER U
UP_CON(5) = INV_CON(4); % Aluminium, production mix, at plant/RER U

% --------------------------- Limestone ---------------------------------
UP_CON(6) = INV_CON(5); % Lime {CH}| production, milled, loose | Alloc Def, U
UP_CON(7) = INV_CON(5); % Limestone, crushed, for mill {CH}| production | Alloc Def, U
UP_CON(8) = INV_CON(5); % Limestone, unprocessed {CH}| limestone quarry operation | Alloc Def, U

UP_CON(9) = INV_CON(6); % Steel, chromium steel 18/8, hot rolled {RER}| production | Alloc Def, U
UP_CON(10) = INV_CON(7); % Flat glass, uncoated {RER}| production | Alloc Def, U

% -------------------------- Copper -------------------------------------
UP_CON(11) = INV_CON(8); % Copper {RER}| production, primary | Alloc Def, U
UP_CON(12) = INV_CON(8); % Copper concentrate {RER}| copper mine operation | Alloc Def, U

UP_CON(13) = INV_CON(9); % Synthetic rubber {RER}| production | Alloc Def, U

% ------------------------- Rock wool -----------------------------------
UP_CON(14) = INV_CON(10); % Rock wool {CH}| production | Alloc Def, U
UP_CON(15) = INV_CON(10); % Rock wool, packed {CH}| production | Alloc Def, U

UP_CON(16) = INV_CON(11); % Bitumen, at refinery/RER U

% ------------------------ Polyethylene ---------------------------------
UP_CON(18) = 0.5 * INV_CON(12); % Polyethylene, high density, granulate {RER}| production | Alloc Def, U
UP_CON(17) = 0.5 * INV_CON(13); % Polyethylene, low density, granulate {RER}| production | Alloc Def, U

UP_CON(19) = INV_CON(14); % Excavation, hydraulic digger {RER}| processing | Alloc Def, U

% ------------------------ Operation ------------------------------------
% UP_CON(19) = INV_CON(15); % Operation, freight train/RER U
% UP_CON(20) = 0; % Operation, freight train, diesel/RER U

UP_CON(20) = INV_CON(16); % Transport, freight train {US}| diesel | Alloc Def, U
UP_CON(21) = INV_CON(17); % Extrusion, plastic pipes {RER}| production | Alloc Def, U
UP_CON(22) = INV_CON(18); % Transport, freight, lorry 16-32 metric ton, EURO5 {RER}| transport, freight, lorry 16-32 metric ton, EURO5 | Alloc Def, U

% ------------------------ Electricity ----------------------------------
% Electricity source percentages are based on US average in 2013
UP_CON(23) = Electricity{2,1} * INV_CON(19); % Electricity, natural gas, at power plant/US U
UP_CON(24) = Electricity{1,1} * INV_CON(19); % Electricity, hard coal, at power plant/US U
UP_CON(25) = Electricity{4,1} * INV_CON(19); % Electricity, hydropower, at pumped storage power plant/US U
UP_CON(26) = Electricity{7,1} * INV_CON(19); % Electricity, high voltage {GB}| electricity production, oil | Alloc Def, U
UP_CON(27) = Electricity{3,1} * INV_CON(19); % Electricity, nuclear, at power plant/US U
UP_CON(28) = Electricity{6,1} * INV_CON(19); % Electricity, at wind power plant 800kW/RER U
UP_CON(29) = Electricity{5,1} * INV_CON(19); % Electricity, high voltage {CH}| treatment of municipal solid waste, incineration | Alloc Def, U

UP_CON(30) = INV_CON(20); % Chemical, organic {GLO}| production | Alloc Def, U
UP_CON(31) = INV_CON(21); % Chemical, inorganic {GLO}| production | Alloc Def, U

UP_CON(32) = INV_CON(22); % Polyethylene terephthalate, granulate, amorphous {RER}| production | Alloc Def, U
UP_CON(33) = INV_CON(23); % Tetrafluoroethylene {RER}| production | Alloc Def, U
UP_CON(34) = INV_CON(24); % Carbon black {GLO}| production | Alloc Def, U
UP_CON(35) = INV_CON(25); % Mini CHP plant, common components for heat+electricity {CH}| construction | Alloc Def, U
UP_CON(1) = INV_CON(26); %Alumina, at plant/US
UP_CON(36) = INV_CON(27); %Sinter, iron {GLO}|production|Alloc Def, U
UP_CON(37) = INV_CON(28); %Polysulfone {GLO}| polysulfone production, for membrane filtration
UP_CON(38) = INV_CON(29); %Polyvinylfluoride, film {US}| production| Alloc Def, U

for i = 1:length(UP_CON)
    traci_con{i,:} = UP_CON(i).*traci_con{i,:} ./ V_treated;
    %     CML_con{i,:} = UP_CON(i).*CML_con{i,:} ./ V_treated;
    CED_con{i,:} = UP_CON(i).*CED_con{i,:} ./ V_treated;
end
% IMPACT_CON = UP_CON_matrix .* traci ./ V_treated; % Construction impact from each unit process

%% Operational Phase
UP_OP = zeros(16,1); % Initialized Ecoinvent Unit Process matrix
% -------------------- Acetic acid (pick one) ---------------------------
UP_OP(1) = INV_OP(1); % Acetic acid, without water, in 98% solution state {RER}| acetaldehyde oxidation | Alloc Def, U
% UP_OP(2) = 0; % Acetic acid, without water, in 98% solution state {RER}| oxidation of butane | Alloc Def, U
% UP_OP(3) = 0; % Acetic acid, without water, in 98% solution state {RER}| acetic acid production, product in 98% solution state | Alloc Def, U

UP_OP(2) = INV_OP(2); % Methanol {GLO}| production | Alloc Def, U
UP_OP(3) = INV_OP(3); % Iron (III) chloride, without water, in 40% solution state {CH}| iron (III) chloride production, product in 40% solution state | Alloc Def, U
UP_OP(4) = INV_OP(4); % Sodium hypochlorite, without water, in 15% solution state {RER}| sodium hypochlorite production, product in 15% solution state | Alloc Def, U
% ------------------------ Electricity ----------------------------------
% Electricity source percentages are based on US average in 2013
UP_OP(10) = Electricity{2,1} * INV_OP(5); % Electricity, natural gas, at power plant/US U
UP_OP(11) = Electricity{1,1} * INV_OP(5); % Electricity, hard coal, at power plant/US U
UP_OP(12) = Electricity{4,1} * INV_OP(5); % Electricity, hydropower, at pumped storage power plant/US U
UP_OP(13) = Electricity{7,1} * INV_OP(5); % Electricity, high voltage {GB}| electricity production, oil | Alloc Def, U
UP_OP(14) = Electricity{3,1} * INV_OP(5); % Electricity, nuclear, at power plant/US U
UP_OP(15) = Electricity{6,1} * INV_OP(5); % Electricity, at wind power plant 800kW/RER U
UP_OP(16) = Electricity{5,1} * INV_OP(5); % Electricity, high voltage {CH}| treatment of municipal solid waste, incineration | Alloc Def, U
% --------------------------- Chlorine(pick one) ------------------------
UP_OP(5) = INV_OP(6); % Chlorine, gaseous {RER}| chlor-alkali electrolysis, membrane cell | Alloc Def, U
% UP_OP(16) = 0; % Chlorine, gaseous {RER}| chlor-alkali electrolysis, diaphragm cell | Alloc Def, U
UP_OP(6) = INV_OP(6); % Chlorine, gaseous {RER}| chlor-alkali electrolysis, mercury cell | Alloc Def, U

UP_OP(7) = INV_OP(7); % Methyl methacrylate {RER}| production | Alloc Def, U
UP_OP(8) = (INV_OP(8)/1000)*16.09344; % Transport, freight, lorry 16-32 metric ton, EURO5 {RER}| transport, freight, lorry 16-32 metric ton, EURO5 | Alloc Def, U; Assume 10 mi transport (16.09 km), 1000 converts km to metric ton
UP_OP(9) = INV_OP(9); % Citric acid {RER}| production | Alloc Def, U

for i = 1:length(UP_OP)
    traci_op{i,:} = UP_OP(i).*traci_op{i,:} ./ V_treated;
    %     CML_op{i,:} = UP_OP(i).*CML_op{i,:} ./ V_treated;
    CED_op{i,:} = UP_OP(i).*CED_op{i,:} ./ V_treated;
end
% IMPACT_OP = UP_OP_matrix .* CF_OP ./ V_treated;


%% Impacts from Direct Emissions
% Characterization factors (CF) for direct emissions (DE)
% DE = [COD; NH3; NH4; Org-N; P; PO4; CH4; CO2];

for i = 1:length(DEmis)
    IMPACT_DE{i,:} = DEmis(i).*IMPACT_DE{i,:} ./ V_treated;
end
% IMPACT_DE = DEmis .* CF_DE ./ V_treated; % Construction impact from each unit process

%% Avoided Impacts
UP_avoided = zeros(7,1); % Ecoinvent unit process matrix

% ------------------------ Electricity ----------------------------------
% Electricity source percentages are based on US average in 2013
UP_avoided(1) = Electricity{2,1} * E_offset_kWh; % Electricity, natural gas, at power plant/US U
UP_avoided(2) = Electricity{1,1} * E_offset_kWh; % Electricity, hard coal, at power plant/US U
UP_avoided(3) = Electricity{4,1} * E_offset_kWh; % Electricity, hydropower, at pumped storage power plant/US U
UP_avoided(4) = Electricity{7,1} * E_offset_kWh; % Electricity, high voltage {GB}| electricity production, oil | Alloc Def, U
UP_avoided(5) = Electricity{3,1} * E_offset_kWh; % Electricity, nuclear, at power plant/US U
UP_avoided(6) = Electricity{6,1} * E_offset_kWh; % Electricity, at wind power plant 800kW/RER U
UP_avoided(7) = Electricity{5,1} * E_offset_kWh; % Electricity, high voltage {CH}| treatment of municipal solid waste, incineration | Alloc Def, U

for i = 1:length(UP_avoided)
    traci_avoid{i,:} = -UP_avoided(i).*traci_avoid{i,:} ./ V_treated;
    %     CML_avoid{i,:} = UP_avoided(i).*CML_avoid{i,:} ./ V_treated;
    CED_avoid{i,:} = -UP_avoided(i).*CED_avoid{i,:} ./ V_treated;
end
% IMPACT_avoided = - UP_avoided_matrix .* CF_OP ./ V_treated;
end