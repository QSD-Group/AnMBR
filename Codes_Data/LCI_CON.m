function INV_CON = LCI_CON (VC_m3, VEX_m3, M_Steel_kg, M_Membrane_kg, M_LDPE_kg, M_HDPE_kg, M_GAC_kg, id)
% Inputs:
% Volume of concrete, VC_m3_m3 [m^3]
% Volume of excavation, VEX_m3 [m^3]
% Mass of steel, M_Steel_kg [kg]
% Mass of membrane, M_Membrane_kg [kg]
% Mass of GAC, M_GAC_kg [kg]
% step_D (membrane material): input 'PET' or 'PTFE'

% Output: Construction Inventories, INV_CON


%% Construction inventory multiplication factors are from Foley, J.
% Comprehensive life cycle inventories of alternative wastewater treatment
% systems. Water Research 2010, 44, (5), 1654-1666.
% =========================================================================
% (1)excavation, hydraulic digger [m^3]
% does not use multiplication factor, use VEX_m^3 instead

% (2)transport, lorry 28t [t-km] (Note: subscript MF represents the results calculated using the mulplication factors)
Lorry_MF = 49.28380109 * VC_m3;

% (3)transport, freight, rail [t-km]
Rail_MF = 49.28380109 * VC_m3;

% (4)electricity, medium voltage, at grid [kWh]
Electricity_MF = 0.03796545 * VC_m3;

% (5)reinforcing steel, at plant [kg]
Rein_Steel_MF = 77.70945427 * VC_m3;

% (7)tap water, at user [kg]
Water_MF = 122.03433256 * VC_m3;

% (8)aluminum, production mix, cast alloy, at plant [kg]
Aluminum_MF = 0.86856264 * VC_m3;

% (9)limestone, crushed, washed [kg]
Limestone_MF = 21.45272867 * VC_m3;

% (10)chromium steel 18/8, at plant [kg]
% does not use multiplication factor, use M_Steel_kg instead

% (11)glass fibre, at plant (flat glass??) [kg]
Glass_MF = 1.95874968 * VC_m3;

% (12)copper, at plant [kg]
Copper_MF = 0.92313605 * VC_m3;

% (13)synthetic rubber, at plant [kg]
Rubber_MF = 0.88444786 * VC_m3;

% (14)rock wool mat, packed, at plant [kg]
Wool_MF = 0.87445555 * VC_m3;

% (15)chemicals organic, at plant [kg]
Chem_Org_MF =  4.06097873 * VC_m3;

% (16)bitumen, at refinery [kg]
Bitumen_MF = 0.50190985 * VC_m3;

% (17)chemicals inorganic, at plant [kg]
Chem_Inorg_MF = 0.49807840 * VC_m3;

% (18)polyethylene, LDPE, granulate, at plant [kg] (factor: 0.01608445)
% does not use multiplication factor, use M_LDPE_kg instead

% (19)polyethylene, HDPE, granulate, at plant [kg] (facotr: 2.4476044)
% does not use multiplication factor, use M_HDPE_kg instead
% =========================================================================

%% Additional construction inventories include:
% Membrane material (PET, PTFE, ceramic, sintered steel, PES, or PVDF): M_Membrane_kg
% Assuming a 10-year life span for membranes, so multiply mass of membrane
% by 3
M_Membrane_kg1 = M_Membrane_kg*3;
% Granular activated carbon (GAC): M_GAC_kg
% Combined heat and power plant (CHP): quantity=1


%% Generate the construction inventory matrix:
% INV_CON =
% [1-Concrete; 2-Reinforcing steel; 3-Tap water; 4-Aluminum; 5-Limestone;
% 6-Chromium Steel; 7-Flat glass; 8-Copper; 9-Sythetic rubber; 10-Rock wool;
% 11-Bitumen; 12-LDPE; 13-HDPE; 14-Excavation; 15-Operation; 16-Transport;
% 17-Extrusion; 18-Transport; 19-Electricity; 20-Organic chemicals;
% 21-Inorganic chemicals; 22-PET; 23-PTFE; 24-GAC; 25-CHP; 26-Ceramic;
% 27-Sintered Steel; 28-PES; 29-PVDF];

if id==1
    INV_CON = [VC_m3; Rein_Steel_MF; Water_MF; Aluminum_MF; Limestone_MF; ...
        M_Steel_kg; Glass_MF; Copper_MF; Rubber_MF; Wool_MF; ...
        Bitumen_MF; M_LDPE_kg; M_HDPE_kg; VEX_m3; 0; Rail_MF; ...
        0; Lorry_MF; Electricity_MF; Chem_Org_MF; ...
        Chem_Inorg_MF; M_Membrane_kg1; 0; M_GAC_kg; 1; 0; 0; 0; 0];
elseif id==2
    INV_CON = [VC_m3; Rein_Steel_MF; Water_MF; Aluminum_MF; Limestone_MF; ...
        M_Steel_kg; Glass_MF; Copper_MF; Rubber_MF; Wool_MF; ...
        Bitumen_MF; M_LDPE_kg; M_HDPE_kg; VEX_m3; 0; Rail_MF; ...
        0; Lorry_MF; Electricity_MF; Chem_Org_MF; ...
        Chem_Inorg_MF; 0; M_Membrane_kg1; M_GAC_kg; 1; 0; 0; 0; 0];
elseif id==3
    INV_CON = [VC_m3; Rein_Steel_MF; Water_MF; Aluminum_MF; Limestone_MF; ...
        M_Steel_kg; Glass_MF; Copper_MF; Rubber_MF; Wool_MF; ...
        Bitumen_MF; M_LDPE_kg; M_HDPE_kg; VEX_m3; 0; Rail_MF; ...
        0; Lorry_MF; Electricity_MF; Chem_Org_MF; ...
        Chem_Inorg_MF; 0; 0; M_GAC_kg; 1; M_Membrane_kg1; 0; 0; 0];
elseif id==4
    INV_CON = [VC_m3; Rein_Steel_MF; Water_MF; Aluminum_MF; Limestone_MF; ...
        M_Steel_kg; Glass_MF; Copper_MF; Rubber_MF; Wool_MF; ...
        Bitumen_MF; M_LDPE_kg; M_HDPE_kg; VEX_m3; 0; Rail_MF; ...
        0; Lorry_MF; Electricity_MF; Chem_Org_MF; ...
        Chem_Inorg_MF; 0; 0; M_GAC_kg; 1; 0; M_Membrane_kg1; 0; 0];
elseif id==5
    INV_CON = [VC_m3; Rein_Steel_MF; Water_MF; Aluminum_MF; Limestone_MF; ...
        M_Steel_kg; Glass_MF; Copper_MF; Rubber_MF; Wool_MF; ...
        Bitumen_MF; M_LDPE_kg; M_HDPE_kg; VEX_m3; 0; Rail_MF; ...
        0; Lorry_MF; Electricity_MF; Chem_Org_MF; ...
        Chem_Inorg_MF; 0; 0; M_GAC_kg; 1; 0; 0; M_Membrane_kg1; 0];
elseif id==6
    INV_CON = [VC_m3; Rein_Steel_MF; Water_MF; Aluminum_MF; Limestone_MF; ...
        M_Steel_kg; Glass_MF; Copper_MF; Rubber_MF; Wool_MF; ...
        Bitumen_MF; M_LDPE_kg; M_HDPE_kg; VEX_m3; 0; Rail_MF; ...
        0; Lorry_MF; Electricity_MF; Chem_Org_MF; ...
        Chem_Inorg_MF; 0; 0; M_GAC_kg; 1; 0; 0; 0; M_Membrane_kg1];
end

end