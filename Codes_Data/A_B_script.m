%% Initialize variables
clear
tic
trials = 1;
Year = 30;
current_system = 'A_B';
LCA_impacts = {'Ozone','Global_Warming','Smog','Acidification','Eutrophication','Carbon_Emis','NonCarbon_Emis','Respiratory_Disease','Ecotoxicity'};
Const_names = {'CON_sum_mean_OzoneDepletion','CON_sum_mean_GlobalWarming','CON_sum_mean_Smog','CON_sum_mean_Acidification','CON_sum_mean_Eutrophication','CON_sum_mean_Carcinogenics','CON_sum_mean_NonCarcinogenics','CON_sum_mean_RespiratoryEffects','CON_sum_mean_Ecotoxicity','CON_sum_mean_FossilFuelDepletion','CON_sum_mean_AbioticDepletion','CON_sum_mean_AbioticDepletion_fossilFuels_','CON_sum_mean_GlobalWarming_GWP100a_','CON_sum_mean_OzoneLayerDepletion_ODP_','CON_sum_mean_HumanToxicity','CON_sum_mean_FreshWaterAquaticEcotox_','CON_sum_mean_MarineAquaticEcotoxicity','CON_sum_mean_TerrestrialEcotoxicity','CON_sum_mean_PhotochemicalOxidation','CON_sum_mean_Acidification_CML','CON_sum_mean_Eutrophication_CML','CON_sum_mean_NonRenewable_Fossil','CON_sum_mean_Non_renewable_Nuclear','CON_sum_mean_Non_renewable_Biomass','CON_sum_mean_Renewable_Biomass','CON_sum_mean_Renewable_Wind_Solar_Geothe','CON_sum_mean_Renewable_Water','CON_sum_var_OzoneDepletion','CON_sum_var_GlobalWarming','CON_sum_var_Smog','CON_sum_var_Acidification','CON_sum_var_Eutrophication','CON_sum_var_Carcinogenics','CON_sum_var_NonCarcinogenics','CON_sum_var_RespiratoryEffects','CON_sum_var_Ecotoxicity','CON_sum_var_FossilFuelDepletion','CON_sum_var_AbioticDepletion','CON_sum_var_AbioticDepletion_fossilFuels_','CON_sum_var_GlobalWarming_GWP100a_','CON_sum_var_OzoneLayerDepletion_ODP_','CON_sum_var_HumanToxicity','CON_sum_var_FreshWaterAquaticEcotox_','CON_sum_var_MarineAquaticEcotoxicity','CON_sum_var_TerrestrialEcotoxicity','CON_sum_var_PhotochemicalOxidation','CON_sum_var_Acidification_CML','CON_sum_var_Eutrophication_CML','CON_sum_var_NonRenewable_Fossil','CON_sum_var_Non_renewable_Nuclear','CON_sum_var_Non_renewable_Biomass','CON_sum_var_Renewable_Biomass','CON_sum_var_Renewable_Wind_Solar_Geothe','CON_sum_var_Renewable_Water','CON_sum_median_OzoneDepletion','CON_sum_median_GlobalWarming','CON_sum_median_Smog','CON_sum_median_Acidification','CON_sum_median_Eutrophication','CON_sum_median_Carcinogenics','CON_sum_median_NonCarcinogenics','CON_sum_median_RespiratoryEffects','CON_sum_median_Ecotoxicity','CON_sum_median_FossilFuelDepletion','CON_sum_median_AbioticDepletion','CON_sum_median_AbioticDepletion_fossilFuels_','CON_sum_median_GlobalWarming_GWP100a_','CON_sum_median_OzoneLayerDepletion_ODP_','CON_sum_median_HumanToxicity','CON_sum_median_FreshWaterAquaticEcotox_','CON_sum_median_MarineAquaticEcotoxicity','CON_sum_median_TerrestrialEcotoxicity','CON_sum_median_PhotochemicalOxidation','CON_sum_median_Acidification_CML','CON_sum_median_Eutrophication_CML','CON_sum_median_NonRenewable_Fossil','CON_sum_median_Non_renewable_Nuclear','CON_sum_median_Non_renewable_Biomass','CON_sum_median_Renewable_Biomass','CON_sum_median_Renewable_Wind_Solar_Geothe','CON_sum_median_Renewable_Water'};
OM_names = {'OP_sum_mean_OzoneDepletion','OP_sum_mean_GlobalWarming','OP_sum_mean_Smog','OP_sum_mean_Acidification','OP_sum_mean_Eutrophication','OP_sum_mean_Carcinogenics','OP_sum_mean_NonCarcinogenics','OP_sum_mean_RespiratoryEffects','OP_sum_mean_Ecotoxicity','OP_sum_mean_FossilFuelDepletion','OP_sum_mean_AbioticDepletion','OP_sum_mean_AbioticDepletion_fossilFuels_','OP_sum_mean_GlobalWarming_GWP100a_','OP_sum_mean_OzoneLayerDepletion_ODP_','OP_sum_mean_HumanToxicity','OP_sum_mean_FreshWaterAquaticEcotox_','OP_sum_mean_MarineAquaticEcotoxicity','OP_sum_mean_TerrestrialEcotoxicity','OP_sum_mean_PhotochemicalOxidation','OP_sum_mean_Acidification_CML','OP_sum_mean_Eutrophication_CML','OP_sum_mean_NonRenewable_Fossil','OP_sum_mean_Non_renewable_Nuclear','OP_sum_mean_Non_renewable_Biomass','OP_sum_mean_Renewable_Biomass','OP_sum_mean_Renewable_Wind_Solar_Geothe','OP_sum_mean_Renewable_Water','OP_sum_var_OzoneDepletion','OP_sum_var_GlobalWarming','OP_sum_var_Smog','OP_sum_var_Acidification','OP_sum_var_Eutrophication','OP_sum_var_Carcinogenics','OP_sum_var_NonCarcinogenics','OP_sum_var_RespiratoryEffects','OP_sum_var_Ecotoxicity','OP_sum_var_FossilFuelDepletion','OP_sum_var_AbioticDepletion','OP_sum_var_AbioticDepletion_fossilFuels_','OP_sum_var_GlobalWarming_GWP100a_','OP_sum_var_OzoneLayerDepletion_ODP_','OP_sum_var_HumanToxicity','OP_sum_var_FreshWaterAquaticEcotox_','OP_sum_var_MarineAquaticEcotoxicity','OP_sum_var_TerrestrialEcotoxicity','OP_sum_var_PhotochemicalOxidation','OP_sum_var_Acidification_CML','OP_sum_var_Eutrophication_CML','OP_sum_var_NonRenewable_Fossil','OP_sum_var_Non_renewable_Nuclear','OP_sum_var_Non_renewable_Biomass','OP_sum_var_Renewable_Biomass','OP_sum_var_Renewable_Wind_Solar_Geothe','OP_sum_var_Renewable_Water','OP_sum_median_OzoneDepletion','OP_sum_median_GlobalWarming','OP_sum_median_Smog','OP_sum_median_Acidification','OP_sum_median_Eutrophication','OP_sum_median_Carcinogenics','OP_sum_median_NonCarcinogenics','OP_sum_median_RespiratoryEffects','OP_sum_median_Ecotoxicity','OP_sum_median_FossilFuelDepletion','OP_sum_median_AbioticDepletion','OP_sum_median_AbioticDepletion_fossilFuels_','OP_sum_median_GlobalWarming_GWP100a_','OP_sum_median_OzoneLayerDepletion_ODP_','OP_sum_median_HumanToxicity','OP_sum_median_FreshWaterAquaticEcotox_','OP_sum_median_MarineAquaticEcotoxicity','OP_sum_median_TerrestrialEcotoxicity','OP_sum_median_PhotochemicalOxidation','OP_sum_median_Acidification_CML','OP_sum_median_Eutrophication_CML','OP_sum_median_NonRenewable_Fossil','OP_sum_median_Non_renewable_Nuclear','OP_sum_median_Non_renewable_Biomass','OP_sum_median_Renewable_Biomass','OP_sum_median_Renewable_Wind_Solar_Geothe','OP_sum_median_Renewable_Water'};
interest = 0.08;

%Initial conditions
Q_mgd = latin_hs(20,5,trials,1); % [mgd] (Normal dist., 25% std. dev)
MLSS = latin_hs(1210,120,trials,1); % [mg/L] (normal dist, 10% std. dev)
COD = 300; % [mg/L]
TSS = 18.3; % [mg/L]
X_V_e = 15; % [mg VSS/L]
X_V_w = 10000; % [mg VSS/L]
X_inert = 100;
S_SO = COD; % [mg-COD/L]
SLR = 20; % [lb/ft2-d], assumed
SF = 20; % [safety factor, unitless]
J = lhs_triangle(5,12,17,trials); % [L/m^2-hr] (Triang. dist.)


%Kinetic and stoichiometric parameters
q_hat = 12; % [mg BOD/mg VSS-d]
K = 20; % [mg COD/L] (Smith, 2014)
Y = 0.5; % [mg TSS/mg BOD] (Smith, 2014)
b = 0.396; % [1/d] (Smith, 2014)
f_d = 0.2; % [unitless] (Smith, 2014)
q_UAP = 1.8; % [mg COD/mg VSS-d], following values from Rittmann and McCarty, pg. 350
q_BAP = 0.1; % [mg COD/mg VSS-d]
K_UAP = 100; % [mg COD/L]
K_BAP = 85; % [mg COD/L]
k1 = 0.12; % [mg COD/mg BOD]
k2 = 0.09; % [mg COD/mg VSS-d]

% Effluent Values (not calculated)
NH3_eff = 0; % [mg-N/L]
NH4_eff = 25; % [mg-N/L] % Hospido,2012
Org_N_eff = 0; % [mg-N/L]
P_eff = 0; % [mg-N/L]
PO4_eff = 3; % [mg-N/L] % Hospido,2012
CO2_eff = 127 * 10^-3; % [mg-CO2/L] % Hospido,2012

%Tank dimensions
w = 21; % [ft], assumed, consistent with rest of code
d = 12; % [ft], assumed

wait=waitbar(0,'');

for i = 1:trials
    V_treated = Q_mgd(i) * 3785.41178 * Year * 365; % [m3]
    %% Calculations
    Q_LD = Q_mgd(i)*3785411.78; %[L/d], conversion is L/10^6 gallons
    % SRT_min = 1/(Y*q_hat-b); % [d]
    % SRT = round(SRT_min*SF); % [d]
    SRT_A = 0.5; %[d], from Wett et al., 2007
    SRT_B = 10; % [d], from Wett et al., 2007
    
    %% A Stage
    S_A = K*((1+b*SRT_A)/(Y*q_hat*SRT_A-(1+b*SRT_A))); % [mg/L]
    % HRT = (SRT_B/MLSS(i))*(X_inert+((1+(1-f_d)*b*SRT_B)/(1+b*SRT_B))*Y*(S_SO-S)); % [d]
    HRT_A = 0.5/24; % [d], from Wett et al., 2007
    V_A = Q_LD*HRT_A; % [L]
%     Q_was_A = (((MLSS(i)*V_A)/SRT_A)-(Q_LD*X_V_e))/(X_V_w-X_V_e); % [L/d]
    Q_was_A = (COD*0.607*Q_LD)/(X_V_w*1.42); %[L/d]
    Q_ras_A = ((X_V_w*Q_was_A)+(X_V_e*(Q_LD-Q_was_A))-(MLSS(i)*Q_LD))/(MLSS(i)-X_V_w); % [L/d]
    X_a_A = ((SRT_A/HRT_A)*Y*(S_SO-S_A))/(1+b*SRT_A); % [mg VSS/L]
    r_ut_A = (S_SO-S_A)/HRT_A; % [mg/L-d]
    
    %SMP calc
    UAP_A = (-(q_UAP*X_a_A*HRT_A+K_UAP+k1*r_ut_A*HRT_A)+sqrt(((q_UAP*X_a_A*HRT_A+K_UAP+k1*r_ut_A*HRT_A)^2)+4*K_UAP*k1*r_ut_A*HRT_A))/2; % [mg COD/L]
    BAP_A = (-(K_BAP+(q_BAP-k2)*X_a_A*HRT_A)+sqrt(((K_BAP+(q_BAP-k2)*X_a_A*HRT_A)^2)+4*K_BAP*k2*X_a_A*HRT_A))/2; % [mg COD/L]
    SMP_A = UAP_A+BAP_A; % [mg COD/L]
    
    %Clarifier sizing
    A_Clarifier_A = ((Q_LD+Q_ras_A)*MLSS(i)*2.20462e-6)/SLR; % [ft^2], rectangular with floor slope of 1% (Clarifier Design, 2nd Ed., WEF, MOP FD-8), 2.20462e-6 converts mg to lb
    
    %Reactor sizing
    L_tank_A = (V_A*0.0353147)/(w*d); % [ft], 0.0353147 converts L to ft3
    
    %Sizing of other buildings
    W_PB = 38 + 4/12; % [ft]
    W_BB = 22; % [ft]
    L_BB = 76 + 8/12; % [ft]
    t_wall = 1; % [ft]
    t_slab = t_wall + 2/12; % [ft]
    
    %O2 requirement
%     dMLSS_dt = (MLSS(i)*V_A)/SRT_A; % [mg VSS/d]
%     dX_dt_bio = dMLSS_dt-X_inert*Q_LD; % [mg VSS/d]
%     Inf_O2 = Q_LD*S_SO+Q_LD*X_inert*1.42; % [mg O2/d], 1.42 converts VSS to O2 equiv.
%     Eff_O2 = Q_LD*(S_A+SMP_A)+dMLSS_dt*1.42; % [mg O2/d]
%     O2_uptake = Inf_O2-Eff_O2; % [mg O2/d]
    O2_uptake = 0.29*COD*Q_LD; %[mg O2/d] Based on data from WWTP in Strass, Austria. 5% of COD is degraded in A stage from aerobic degradation
    O2_cfm = O2_uptake*22.4*0.0353147/(1000*32*1440);% [cfm O2]
    [~, M_SS_blower, OD_gh, OD_gsm] = Gas_Sparging_Submerged(O2_cfm, L_tank_A, W_PB, L_BB);
    Power_reqd = (O2_uptake*(10^(-6)))/24; % [kW], assumes 1 kWh/1 kg O2, Rittmann and McCarty, 2001, pg. 352
    P_input_blower_A = Power_reqd;
    
    % Calculate number of trains (which we're assuming is the number of clarifiers)
    N_train_A=2;
    Module_SA = 370; %Membrane surface area per Zenon Membrane Module [ft2]
    Mod_per_cas = 44; %Zenon 500D Modules per Cassette (Max = 48)
    Contact_SA = Module_SA*Mod_per_cas; %Contact Surface Area per Zenon Membrane Cassette [ft2]
    Cas_per_tank = 16; %Cassettes per membrane tank
    Spare_cassettes = 2; %Spare cassettes per membrane tank
    Flux_oneoff = (Q_mgd*1000000)/((N_train_A-1)*Cas_per_tank*Contact_SA);
    while Flux_oneoff>J*0.5754 %J converted to gal/ft^2-d by multiplying by 24 hr/d and dividing by (10.7639 ft2/m2 and 3.785 L/gal)
        Mod_per_cas=Mod_per_cas+1;
        if Mod_per_cas==49
            if Cas_per_tank==23
                N_train_A=N_train_A+1;
                Cas_per_tank=16;
                Mod_per_cas = 44;
            else
                Cas_per_tank = Cas_per_tank+1;
                Mod_per_cas = 44;
            end
        end
        Contact_SA = Module_SA*Mod_per_cas;
        Flux_oneoff = (Q_mgd*1000000)/((N_train_A-1)*Cas_per_tank*Contact_SA);
    end
    
    N_cassettes = Cas_per_tank*N_train_A;
    if Cas_per_tank<11
        L_membrane_unit = ceil(Cas_per_tank*8.5+Spare_cassettes*8.5); %[ft]
        W_membrane_unit = 10; %[ft]
    else
        L_membrane_unit = ceil(Cas_per_tank*8.5/2+Spare_cassettes*8.5/2); %[ft]
        W_membrane_unit = 21; %[ft]
    end
    
    Tank_footprint_A = N_train_A*L_membrane_unit*W_membrane_unit;
    
    %Pumping
    P_input_LIFT_A = 0; % Eliminated lift pumping 11/20/15
    M_SS_LIFT_A = 0;
    Q_cfh = Q_mgd(i) * 133681 / 24; % [cfh]
    L_MT_A = A_Clarifier_A/(w*N_train_A); % [ft]
    %% B Stage
%     S_B = K*((1+b*SRT_B)/(Y*q_hat*SRT_B-(1+b*SRT_B))); % [mg/L]
%     HRT_B = (SRT_B/MLSS(i))*(X_inert+((1+(1-f_d)*b*SRT_B)/(1+b*SRT_B))*Y*(S_A-S_B)); % [d]
%     % HRT = 2/24; % [d]
%     V_B = Q_LD*HRT_B; % [L]
% %     Q_was_B = (((MLSS(i)*V_B)/SRT_B)-(Q_LD*X_V_e))/(X_V_w-X_V_e); % [L/d]
%     Q_was_B = (COD*0.141*Q_LD)/(X_V_w*1.42); %[L/d]
%     Q_ras_B = ((X_V_w*Q_was_B)+(X_V_e*(Q_LD-Q_was_B))-(MLSS(i)*Q_LD))/(MLSS(i)-X_V_w); % [L/d]
%     X_a_B = ((SRT_B/HRT_B)*Y*(S_A-S_B))/(1+b*SRT_B); % [mg VSS/L]
%     r_ut_B = (S_A-S_B)/HRT_B; % [mg/L-d]
%     
%     %SMP calc
%     UAP_B = (-(q_UAP*X_a_B*HRT_B+K_UAP+k1*r_ut_B*HRT_B)+sqrt(((q_UAP*X_a_B*HRT_B+K_UAP+k1*r_ut_B*HRT_B)^2)+4*K_UAP*k1*r_ut_B*HRT_B))/2; % [mg COD/L]
%     BAP_B = (-(K_BAP+(q_BAP-k2)*X_a_B*HRT_B)+sqrt(((K_BAP+(q_BAP-k2)*X_a_B*HRT_B)^2)+4*K_BAP*k2*X_a_B*HRT_B))/2; % [mg COD/L]
%     SMP_B = UAP_B+BAP_B; % [mg COD/L]
%     
%     %Clarifier sizing
%     A_Clarifier_B = ((Q_LD+Q_ras_B)*MLSS(i)*2.20462e-6)/SLR; % [ft^2], rectangular with floor slope of 1% (Clarifier Design, 2nd Ed., WEF, MOP FD-8), 2.20462e-6 converts mg to lb
%     
%     %Reactor sizing
%     L_tank_B = (V_B*0.0353147)/(w*d); % [ft], 0.0353147 converts L to ft3
%     
%     %Sizing of other buildings
% %     W_PB = 38 + 4/12; % [ft]
% %     W_BB = 22; % [ft]
% %     L_BB = 76 + 8/12; % [ft]
% %     t_wall = 1; % [ft]
% %     t_slab = t_wall + 2/12; % [ft]
%     
%     %O2 requirement
%     % dMLSS_dt = (MLSS(i)*V)/SRT_B; % [mg VSS/d]
%     % dX_dt_bio = dMLSS_dt-X_inert*Q_LD; % [mg VSS/d]
%     % Inf_O2 = Q_LD*S_SO+Q_LD*X_inert*1.42; % [mg O2/d], 1.42 converts VSS to O2 equiv.
%     % Eff_O2 = Q_LD*(S+SMP)+dMLSS_dt*1.42; % [mg O2/d]
%     % O2_uptake = Inf_O2-Eff_O2; % [mg O2/d]
%     O2_uptake = 0.175*COD*Q_LD; %[mg O2/d] Based on data from WWTP in Strass, Austria. 17.5% of COD is degraded in B stage from aerobic degradation
%     O2_cfm = O2_uptake*22.4*0.0353147/(1000*32*1440);% [cfm O2]
%     [P_input_blower, M_SS_blower, OD_gh, OD_gsm] = Gas_Sparging_Submerged(O2_cfm, L_tank_B, W_PB, L_BB);
%     Power_reqd = (O2_uptake*(10^(-6)))/24; % [kW], assumes 1 kWh/1 kg O2, Rittmann and McCarty, 2001, pg. 352
%     P_input_blower_B = Power_reqd;
%     
%     % Calculate number of trains (which we're assuming is the number of clarifiers)
%     N_train_B=2;
% %     Module_SA = 370; %Membrane surface area per Zenon Membrane Module [ft2]
% %     Mod_per_cas = 44; %Zenon 500D Modules per Cassette (Max = 48)
% %     Contact_SA = Module_SA*Mod_per_cas; %Contact Surface Area per Zenon Membrane Cassette [ft2]
%     Cas_per_tank = 16; %Cassettes per membrane tank
%     Spare_cassettes = 2; %Spare cassettes per membrane tank
%     Flux_oneoff = (Q_mgd*1000000)/((N_train_B-1)*Cas_per_tank*Contact_SA);
%     while Flux_oneoff>J*0.5754 %J converted to gal/ft^2-d by multiplying by 24 hr/d and dividing by (10.7639 ft2/m2 and 3.785 L/gal)
%         Mod_per_cas=Mod_per_cas+1;
%         if Mod_per_cas==49
%             if Cas_per_tank==23
%                 N_train_B=N_train_B+1;
%                 Cas_per_tank=16;
%                 Mod_per_cas = 44;
%             else
%                 Cas_per_tank = Cas_per_tank+1;
%                 Mod_per_cas = 44;
%             end
%         end
%         Contact_SA = Module_SA*Mod_per_cas;
%         Flux_oneoff = (Q_mgd*1000000)/((N_train_B-1)*Cas_per_tank*Contact_SA);
%     end
%     
%     N_cassettes = Cas_per_tank*N_train_B;
%     if Cas_per_tank<11
%         L_membrane_unit = ceil(Cas_per_tank*8.5+Spare_cassettes*8.5); %[ft]
%         W_membrane_unit = 10; %[ft]
%     else
%         L_membrane_unit = ceil(Cas_per_tank*8.5/2+Spare_cassettes*8.5/2); %[ft]
%         W_membrane_unit = 21; %[ft]
%     end
%     
%     Tank_footprint_B = N_train_B*L_membrane_unit*W_membrane_unit;
%     
%     %Pumping
%     P_input_LIFT_B = 0;
%     M_SS_LIFT_B = 0;
    
Q_was_B = 0;
S_B =0;
A_Clarifier_B = 0;
N_train_B=0;
L_tank_B = 0;
M_SS_LIFT_B = 0;
P_input_LIFT_B = 0;
P_input_blower_B = 0;
Q_ras_B = 0;
    %% AD
    %Assume Qin for AD is Qwas, SRT = HRT
    %6 m side depth (Smith, 2014)
    
    AD_side = 10; % [m]
    AD_SRT = 30; % [d], approximate number based on Strass Plant, Wett et al., 2007, Figure 8
    AD_HRT = AD_SRT; % [d]
    X_thickened_WAS = 25000; %[mg/L]
    Q_was_thickened_LD = (Q_was_A+Q_was_B)*(X_V_w/X_thickened_WAS);
    AD_Q_in = Q_was_thickened_LD/1000; % [m3/d]
    AD_V = AD_Q_in*AD_HRT; % [m3]
    AD_SA = (AD_V/AD_side); % [m2]
    AD_diameter = sqrt(AD_SA*4/pi); % [m]
    
    % Heating
    %Heat transfer coefficients
    H_wall = 0.7; % Walls with insulation, [W/m2-deg C]
    H_floor = 1.7; % Floor in contact with dry earth, [W/m2-deg C]
    H_ceiling = 0.95; % Floating cover with insulation, [W/m2-deg C]
    SH_sludge = 4200; % Specific heat of sludge, [J/kg-deg C]
    
    %Temps
    T_air = 17; % deg C
    T_earth = 10; % deg C
    T_sludge = 25; % deg C, assuming sludge is at ambient air temp
    T_AD = 35; % deg C
    
    %Sludge calcs
    S_mass = AD_Q_in*10^3; % [kg/d]
%     S_vol = S_mass/(1.02*10^3*0.05); % [m3/d], Assumes sludge is 95% moisture, has a specific gravity of 1.02, 10^3 = specific weight of water
%     COD_loading = 0.743*COD*Q_LD/10^6; % [kg/d],  Wett et al., 2007
%     Vol_loading = COD_loading/AD_V; % [kg/m3-d]
    COD_degraded_AD = 0.354*COD*Q_LD/10^6; %[kg/d], 35.4% influent COD removal, Wett et al., 2007
    VSS_produced = (0.08*COD_degraded_AD)/(1+0.03*AD_SRT); % [kg/d] 0.08 = Y, 0.03 = b, 0.7 = efficiency of waste conversion (Metcalf, pg 1508)
    Q_CH4 = 0.4*COD_degraded_AD-1.42*VSS_produced; % [m3/d]
    Q_CH4_burned = 0;
    increment = Q_CH4*0.01;
    Q_CH4_1 = Q_CH4;
    
    H_req = S_mass*(T_AD-T_sludge)*SH_sludge; % [J/d]
    A_wall = pi*AD_diameter*AD_side; % [m2]
    A_floor = pi*AD_diameter^2/4; % [m2]
    A_roof = pi*AD_diameter^2/4; % [m2]
    Heat_loss_wall = H_wall*A_wall*(T_AD-T_air)*86400; % [J/d], 86,400 convert seconds to days
    Heat_loss_floor = H_floor*A_floor*(T_AD-T_earth)*86400;
    Heat_loss_roof = H_ceiling*A_roof*(T_AD-T_air)*86400;
    Heat_loss_building = Heat_loss_wall+Heat_loss_floor+Heat_loss_roof;
    Heating_total = H_req+Heat_loss_building;
    [P_offset, heat] = CHP(Q_CH4);
    
    while Heating_total > heat
        if Q_CH4 <= 0
            Q_CH4 = 0;
            [P_offset, ~] = CHP(Q_CH4);
            heat = Q_CH4_burned*35846*1000;
            break
        elseif Q_CH4_burned == 0
            Q_CH4_burned = increment;
            Q_CH4 = Q_CH4-increment;
        else
            Q_CH4_burned = Q_CH4_burned+increment;
            Q_CH4 = Q_CH4-increment;
        end
        [P_offset, heat] = CHP(Q_CH4);
        heat = heat + Q_CH4_burned*35846*1000; % 35846 kJ/m3 CH4
    end
    
    % --> Concrete
    % Concrete of anaerobic filter
    FB_AD = 3; % [ft] freeboard
    t_wall_AD = 6/12; % [ft] wall thickness
    t_slab_AD = 8/12; % [ft] slab thickness
    
    % External wall (wall concrete), VWC_E_AF [ft^3]
    VWC_AD = t_wall_AD * pi * AD_diameter * (AD_side + FB_AD);
    
    % Floor and cover (slab concrete), VSC_F_AF [ft^3]
    VSC_AD =2* t_slab_AD * (pi/4) * AD_diameter^2;
    
    % --> Excavation
    SL = 1.5; % Slope = horizontal/vertical
    CA = 3; % [ft] Construction Access
    
    %   Excavation of Pump Building
    AD_PBL = 50; % [ft] Pump Building Length
    AD_PBW = 30; % [ft] Pump Building Width
    AD_PBD = 10; % [ft] Pump Building Depth
    AD_Area_B_P = (AD_PBL + 2 * CA) * (AD_PBW + 2 * CA); % [ft^2] Bottom Area of frustum
    AD_Area_T_P = (AD_PBL + 2 * CA + AD_PBW * SL) * (AD_PBW + 2 * CA + AD_PBD * SL); % [ft^2] Top Area of frustum
    AD_VEX_PB = 0.5 * (AD_Area_B_P + AD_Area_T_P) * AD_PBD; % [ft^2] Volume of excavaion of Pump Building
    
    %% Sludge handling
    Q_WAS_mgd = (Q_was_A+Q_was_B)/3785411.78;
    Q_cmd = Q_mgd(i) * 3785.41178; % [m^3/day] Unit conversion
    M_COD_removed = (S_SO - S_B)/10^3 * Q_cmd; % [kg-COD/day] Daily COD removal
    [P_input_GBT, M_SS_GBT, M_COD_WAS] = Sludge_Handling (Q_CH4_1, Q_WAS_mgd, M_COD_removed);
    Q_ras_A_mgd = Q_ras_A/3785411.78;
    [P_input_RAS, ~, ~, ~] = Pumping(Q_ras_A_mgd, 1, 0, 50, 50, 0);
    P_input_GBT = P_input_GBT+P_input_RAS;
    
    %Air (H&S)
    TCFM=O2_cfm;
    if (TCFM<=30000)
        N=1;
        CFMB = TCFM/N;
        while (CFMB>7500)
            N=N+1;
            CFMB = TCFM/N;
        end
        NB = N +1;
    elseif (TCFM>30000)&&(TCFM<=72000)
        N=1;
        CFMB = TCFM/N;
        while (CFMB>18000)
            N=N+1;
            CFMB = TCFM/N;
        end
        NB = N +1;
    elseif (TCFM>72000)
        N=1;
        CFMB = TCFM/N;
        while (CFMB>100000)
            N=N+1;
            CFMB = TCFM/N;
        end
        NB = N +1;
    end
    
    L_MT_B = A_Clarifier_B/(w*N_train_B); % [ft]
    
    %% Concrete
    % Concrete [ft^3] - Distribution Channel
    W_dist = 4.5; % [ft] Width of distribution channel
    W_N_trains_A = (w + 2 * t_wall) * N_train_A - t_wall * (N_train_A - 1);
    W_N_trains_B = (w + 2 * t_wall) * N_train_B - t_wall * (N_train_B - 1);
    VWC_dist = d * t_wall * (2 * W_N_trains_A + 2 * W_N_trains_B + 2 * W_dist);
    VSC_dist = W_N_trains_A * (W_dist + 2 * t_wall)+ W_N_trains_A * (W_dist + 2 * t_wall) * t_slab + W_N_trains_B * (W_dist + 2 * t_wall)+ W_N_trains_B * (W_dist + 2 * t_wall) * t_slab;
    
    % Concrete [ft^3] - CSTR Trains
    VWC_CSTR = d * t_wall * (N_train_A + 1) * L_tank_A+d * t_wall * (N_train_B + 1) * L_tank_B;
    VSC_CSTR = N_train_A * (d + 2.4 + 4.81 + 2.4 + 7.26) * L_tank_A +N_train_A * (d + 2.4 + 4.81 + 2.4 + 7.26) * t_slab * L_tank_A +N_train_B * (d + 2.4 + 4.81 + 2.4 + 7.26) * L_tank_B +N_train_B * (d + 2.4 + 4.81 + 2.4 + 7.26) * t_slab * L_tank_B;
    
    % Concrete [ft^3] - Membrane Trains
    VWC_MT = d * t_wall * (N_train_A + 1) * L_MT_A + d * t_wall * (N_train_B + 1) * L_MT_B;
    VSC_MT = N_train_A * (d + 2.4 + 4.81 + 2.4 + 7.26) * L_MT_A+N_train_A * (d + 2.4 + 4.81 + 2.4 + 7.26) * t_slab * L_MT_A+ N_train_B * (d + 2.4 + 4.81 + 2.4 + 7.26) * L_MT_B+N_train_B * (d + 2.4 + 4.81 + 2.4 + 7.26) * t_slab * L_MT_B;
    
    % Concrete [ft^3] - Effluent Channel
    W_eff = 5; % [ft] Width of effluent channel
    VWC_eff = d * t_wall * (2 * W_N_trains_A + 2*W_N_trains_B +  2 * W_eff);
    VSC_eff = W_N_trains_A * (W_eff + 2 * t_wall)+ W_N_trains_A * (W_eff + 2 * t_wall) * t_slab +W_N_trains_B * (W_eff + 2 * t_wall)+ W_N_trains_B * (W_eff + 2 * t_wall) * t_slab;
    
    % Pump/Blower Building
    VWC_PBB = d * t_wall * (2 * W_N_trains_A + 2 * W_N_trains_B + 2 * W_PB + 2 * W_BB);
    VSC_PBB = W_N_trains_A * (W_PB + t_wall + W_BB)+ W_N_trains_A * (W_PB + t_wall + W_BB) * t_slab + W_N_trains_B * (W_PB + t_wall + W_BB)+W_N_trains_B * (W_PB + t_wall + W_BB) * t_slab;
    
    % Wet Well (for mix liquor storage)
    D_WW = 12; % [ft] Depth of wet well
    W_WW = 8; % [ft] Width of wet well
    L_WW = 8; % [ft] Length of wet well
    VWC_WW = D_WW * (L_WW * t_wall + W_WW * t_wall + 4 * t_wall);
    VSC_WW = t_slab * (L_WW + 2 * t_wall) * (W_WW + 2 * t_wall)+(L_WW + 2 * t_wall) * (W_WW + 2 * t_wall);
    
    % Total Volume of Wall Concrete [f^3]
    VWC = 2 * VWC_dist + VWC_CSTR + VWC_MT + 2 * VWC_eff + VWC_PBB + VWC_WW+VWC_AD;
    
    % Total Volume of Slab Concrete [f^3]
    VSC = 2 * VSC_dist + VSC_CSTR + VSC_MT + 2 * VSC_eff + VSC_PBB + VSC_WW+VSC_AD;
    
    % Volume of excavation of membrane trains [ft^3]
    Area_B_train = (W_dist + L_MT_A + L_MT_B + W_eff + 2 * CA) * (W_N_trains_A + W_N_trains_B + 2 * CA); % top area
    Area_T_train = (W_dist + L_MT_A + L_MT_B + W_eff + 2 * CA + d * SL) * (W_N_trains_A + W_N_trains_B + 2 * CA + d * SL); % bottom area
    VEX_train = 0.5 * (Area_B_train + Area_T_train) * d;
    
    % Volume of excavation of pump/blower building [ft^3]
    Area_B_PB = (W_PB + W_BB + 2 * CA) * (W_N_trains_A + W_N_trains_B + 2 * CA);
    Area_T_P = (W_PB + W_BB + 2 * CA + d * SL) * (W_N_trains_A + W_N_trains_B + 2 * CA + d * SL);
    VEX_PB = 0.5 * (Area_B_PB + Area_T_P) * d;
    
    % Total Volume of Excavation [ft^3]
    VEX = VEX_train + VEX_PB+AD_VEX_PB;
    
    %% Construction inventory
    %Concrete
    VC_m3 = (VSC + VWC) * 0.0283168; % [m^3] Unit conversion from [ft^3] to [m^3]
    
    % Excavation
    VEX_m3 = VEX * 0.0283168; % [m^3] Unit conversion from [ft^3] to [m^3]
    
    % Stainless
    M_Steel_kg = (M_SS_LIFT_A + M_SS_LIFT_B + M_SS_blower + M_SS_GBT);
    
    %Constrution Inventory
    INV_CON = LCI_CON (VC_m3, VEX_m3, M_Steel_kg, 0, 0, 0, 0, 1);
    
    %% O&M Inventory
    
    P_input_LIFT = P_input_LIFT_A + P_input_LIFT_B;
    P_input_blower = P_input_blower_A + P_input_blower_B;
    P_input_tot = P_input_LIFT + P_input_blower + P_input_GBT; % [kW] Total power input
    Power_pct = [P_input_LIFT/P_input_tot, P_input_blower/P_input_tot, P_input_GBT/P_input_tot,];
    
    E_input_kWh = P_input_tot * Year * 365 * 24; % [kWh] Total electricity consumption over N years
    E_input_per_m3 = E_input_kWh/(Q_cmd*Year*365);
    E_offset_kWh = P_offset * Year * 365 * 24; % [kWh] Total electricity offset over N years
    E_offset_per_m3 = E_offset_kWh/(Q_cmd*Year*365);
    
    INV_OP = [0; 0; 0; 0; E_input_kWh; 0; 0; M_COD_WAS; 0]; % Operational inventory matrix
    
    %% Direct Emissions
    % Total emissions over N years
    COD_water = S_B /10^3 * V_treated; % [kg-COD]
    NH3_water = NH3_eff /10^3 * V_treated; % [kg-N]
    NH4_water = NH4_eff /10^3 * V_treated; % [kg-N]
    Org_N_water = Org_N_eff /10^3 * V_treated; % [kg-N]
    P_water = P_eff /10^3 * V_treated; % [kg-N]
    PO4_water = PO4_eff /10^3 * V_treated; % [kg-N]
    CO2_air = CO2_eff /10^3 * V_treated; % [kg-CO2]
    
    DEmis = [COD_water; NH3_water; NH4_water; Org_N_water; P_water; PO4_water; 0; CO2_air]; %The 0 is for methane to the atmosphere, but we're assuming it's either burned by the CHP or flared
    
    %% output1 for Cost Estimation
    %NOTE: Equations are either from H&S (indicated in parentheses) or from
    %Capdet (all others). Some equations from Capdet have been modified to be
    %smooth.
    output1_LCC=struct();
    total_OM=0;
    total_const=0;
    total=0;
    
    output1_LCC.Air_Piping=0;
    output1_LCC.Blowers=0;
    
    %Air piping (CapdetWorks, unchanged and H&S)
    CFMD = O2_cfm; % Design Capacity of Blowers, scfm
    if (CFMD<=1000)
        COSTAP = 617.2 * 3.33*CFMD^0.2553; % 3.33: Air Flow Fraction, calculated as STE/6, and STE stands for Standard Oxygen Transfer Efficiency _%
    elseif (CFMD>1000)&&(CFMD<=10000)       % The default STE value is 20. If users have their own STE value, then plug the value into STE/6.
        COSTAP = 1.43 * 3.33* CFMD^1.1337;  % If STE/6 > 1, take this value and replace 3.33, if the value is smaller than 1, then use 1 to replace 3.33.
    elseif (CFMD>10000)
        COSTAP = 28.59 * 3.33*CFMD^0.8085; % In calculating COSTAP, the CEPCIP (Current CE Plant Cost Index for Pipe, Valves, etc) value is set to be 241
    end
    output1_LCC.Air_Piping=COSTAP;
    total_const=total_const+COSTAP;
    
    %Blowers (CapdetWorks, changed and H&S)
    %     NB=2;
    %     BLW_size =ceil(Air_demand_10_10_instantaneous/(NB-1)); %Blower size to handle 10:10
    %     Tot_cfm=BLW_size*(NB-1);
    %     while Tot_cfm<Air_demand_10_10_instantaneous
    %         NB=NB+1;
    %         Tot_cfm=BLW_size*(NB-1);
    %     end
    if (TCFM<=30000)
        COSTRO = 0.7*CFMB^0.6169; %These equations are changed to smooth out the curve
    elseif (TCFM>30000)&&(TCFM<=72000)
        COSTRO = 0.377* CFMB^0.5928;
    elseif (TCFM>72000)
        COSTRO = 0.964* CFMB^0.4286;
    end
    
    if (TCFM<=30000)
        COSTBE  = 58000*(COSTRO) / 100; % Purchase Cost of Blower of CFMB Capacity, $, 50000:Purchase cost of 3,000 scfm at 8 psig
    elseif (TCFM>30000)&&(TCFM<=72000)
        COSTBE  = 218000*(COSTRO) / 100; % Purchase Cost of 12,000 scfm at 8 psig
    elseif (TCFM>72000)
        COSTBE  = 480000*(COSTRO) / 100; % Purchase Cost of 50,000 scfm at 8 psig
    end
    
    IBC = 2 * COSTBE; % IBC: Installed Blower Costs, $
    BBA = 128* TCFM^0.256; % Blower Building Area, ft^2
    BBC = BBA * 90; % 90: Unit Price Input for Building Costs, $/ft^2
    TBCC_BLW = (IBC * NB*2 + BBC) * 1.11; % 1.11: Correction Factor for Other Minor Construction Costs includes piping, concrete, steel, electrical, paint and installation labor. Multiply NB by 2 to represent replacement after 15 years (Daigger, 2011).
    output1_LCC.Blowers=TBCC_BLW;
    total_const=total_const+TBCC_BLW;
    
    %Pumping (CapdetWorks, changed)
    RSR = (Q_ras_A+Q_ras_B)/Q_LD; %input('Enter the Return Sludge Ratio to Average Wastewater Flow: ');
    Qavg=Q_mgd(i);
    GPMI = (2*Qavg*10^6) / 1440; % GPMI: Design Capacity of Intermediate Pumps in gpm. 2 = excess capacity factor to handle peak flows
    GPMR = (Qavg*RSR*10^6) / 1440; % GPMR: Design Capacity of Return Sludge Pumps in gpm
    FPCI = 1440*GPMI / (10^6); % FPC: Firm Pumping Capacity in mgd
    FPCR = 1440*GPMR / (10^6);
    if (FPCI == 0)
        OLCI = 0;
    elseif (FPCI > 0) && (FPCI <= 7)
        OLCI = 440*25* FPCI^0.1285;
    elseif (FPCI > 7) && (FPCI <= 41) %Changed the upper limit to 41 (from 30)
        OLCI = 294.4*25* FPCI^0.3335;
    elseif (FPCI > 41) && (FPCI <= 80)
        OLCI = 40.5 *25*FPCI^0.8661;
    elseif (FPCI > 80)
        OLCI = 21.3*25*FPCI^1.012;
    end
    if (FPCI == 0)
        MLCI = 0;
    elseif (FPCI > 0) && (FPCI <= 7)
        MLCI = 360 *25*FPCI^0.1478;
    elseif (FPCI > 7) && (FPCI <= 41) %Changed the upper limit to 41 (from 30)
        MLCI = 255.2* 25*FPCI^0.3247;
    elseif (FPCI > 41) && (FPCI <= 80)
        MLCI = 85.7* 25*FPCI^0.6456;
    elseif (FPCI > 80)
        MLCI = 30.6*25* FPCI^0.8806;
    end
    if (FPCR == 0)
        OLCR = 0;
    elseif (FPCR > 0) && (FPCR <= 7)
        OLCR = 440* 25*FPCR^0.1285;
    elseif (FPCR > 7) && (FPCR <= 41) %Changed the upper limit to 41 (from 30)
        OLCR = 294.4* 25*FPCR^0.3335;
    elseif (FPCR > 41) && (FPCR <= 80)
        OLCR = 40.5 *25*FPCR^0.8661;
    elseif (FPCR > 80)
        OLCR = 21.3*25*FPCR^1.012;
    end
    if (FPCR == 0)
        MLCR = 0;
    elseif (FPCR > 0) && (FPCR <= 7)
        MLCR = 360 *25*FPCR^0.1478;
    elseif (FPCR > 7) && (FPCR <= 41) %Changed the upper limit to 41 (from 30)
        MLCR = 255.2*25* FPCR^0.3247;
    elseif (FPCR > 41) && (FPCR <= 80)
        MLCR = 85.7*25* FPCR^0.6456;
    elseif (FPCR > 80)
        MLCR = 30.6*25* FPCR^0.8806;
    end
    TOLC = OLCI + OLCR;
    TMLC = MLCI + MLCR;
    a = 1;
    GPMBI = GPMI/a;
    while (GPMBI > 80000)
        a = a+1;
        GPMBI = GPMI/a;
    end
    if (GPMR == 0)
        b = 0;
    elseif (GPMR>0)
        b = 1;
        GPMBR = GPMR/b;
        while (GPMBR > 80000)
            b = b+1;
            GPMBR = GPMR/b;
        end
    end
    c=2;
    GPMPI = GPMBI / c;
    while (GPMPI >20000)
        c=c+1;
        GPMPI = GPMBI / c;
    end
    c=c+1;
    if (GPMR == 0)
        d = 0;
    elseif (GPMR>0)
        d = 2;
        GPMPR = GPMBR/d;
        while (GPMPR > 20000)
            d = d+1;
            GPMPR = GPMR/d;
        end
        d = d+1;
    end
    PBAI = (0.0284*GPMBI+640)*a;
    PBAR = (0.0284*GPMBR+640)*b;
    COSTPB=(PBAI+PBAR)*90; % 90: Unit price for pump building cost. Information from CapdetWorks 2007.
    
    if (GPMPI>0)&&(GPMPI<5000)
        COSTROI=2.93*GPMPI^0.4404;
    elseif (GPMPI>5000)
        COSTROI=0.0064*GPMPI^1.16;
    end
    if (GPMPR>0)&&(GPMPR<5000)
        COSTROR=2.93*GPMPR^0.4404;
    elseif (GPMPR>5000)
        COSTROR=0.0064*GPMPR^1.16;
    end
    IPC = 2.065*10^5+7.721*10^4*Qavg; % Fit a curve to the data instead of using the equation below
    %IPC = 2*(COSTROI/100) *48549.837*a*c + 2*(COSTROR/100) *48549.837*b*d; % IPC: Installed equipment cost. 48549.837:Standard Cost of a 3000 gpm pump and driver.
    COSTE = VEX/27*0.3; %Cost of earthwork
    TBCC_PUP = (COSTE + COSTPB + IPC*2)*1.18; %Multiply cost of pumps by 2 to represent replacement cost after 15-year useful life is up (Diagger, 2011).
    MSC = TBCC_PUP*(0.007/100);
    TBCC_PUP=TBCC_PUP+MSC;
    output1_LCC.Pumping_Op_Labor=TOLC;
    output1_LCC.Pumping_Main_Labor=TMLC;
    output1_LCC.Pumps=TBCC_PUP;
    total_const=total_const+TBCC_PUP+MSC;
    total_OM=total_OM+TOLC+TMLC;
    
    %Concrete/earthwork (CapdetWorks, unchanged)
    COSTCW = VWC/27 * 650;  % Cost of wall concrete. 18.52: Unit price for wall concrete in dollar, information obtained from the 2007 vendor data in CapdetWorks
    % 27: Conversion Factor from ft^3 to yd^3
    COSTCS = VSC/27*350;    % Cost of slab concrete. 12.96: Unit price for slab concrete in dollar, information obtained from the 2007 vendor data in CapdetWorks
    % 27: Conversion Factor from ft^3 to yd^3
    COSTRC = COSTCW + COSTCS; % Total concrete cost
    COSTE = VEX/27*8.00; %Cost of earthwork
    output1_LCC.Concrete=COSTRC;
    output1_LCC.Earthwork=COSTE;
    total_const=total_const+COSTRC+COSTE;
    
    Power_cost = 0.10; %$/kWh
    Cost_power = P_input_tot * 365 * 24*Power_cost; %[$]
    output1_LCC.Electricity(1,1)=Cost_power;
    total_OM=total_OM+Cost_power;
    
    %Screenings disposal (additional screenings for MBR, i.e., in addition to 6
    %mm screens) (Hazen and Sawyer)
    Screen_rate_max = 4; %[CF/hr/MGD]
    Screen_rate = 0.5;
    Screen_per_day = Screen_rate_max*Q_mgd(i)*24*Screen_rate; %[CF/day]
    Compaction = 0.75;
    Daily_screen = Screen_per_day*(1-Compaction);
    Daily_screen_CY = Daily_screen/27;
    Cost_disposal = 225; %$/20 CY container
    Cost_screening_and_sludge = Daily_screen_CY*365*Cost_disposal/20+(M_COD_WAS/907.18474)*125*365; %$125/ton dewatered sludge disposal, 907.18474 conversion of kg to ton (US)
    output1_LCC.Screening=Cost_screening_and_sludge;
    total_OM=total_OM+Cost_screening_and_sludge;
    
    
    Cost_CHP_const = P_offset * 1225;
    Cost_CHP_OM = P_offset * 365 * 24*.0185;
    output1_LCC.CHP_Const=Cost_CHP_const;
    output1_LCC.CHP_OM=Cost_CHP_OM;
    total_const = total_const+Cost_CHP_const;
    total_OM = total_OM+Cost_CHP_OM;
    
    % Net present value of O&M costs
    Present_OM = total_OM*(((1+interest)^Year-1)/(interest*(1+interest)^Year)); %From Engineering Economic Analysis, 11th edition, Newnan et al.
    
    %Total Cost
    total=Present_OM+total_const;
    output1_LCC.Total_Const=total_const;
    output1_LCC.Total_Present_OM=Present_OM;
    output1_LCC.Total_Yearly_OM=total_OM;
    output1_LCC.Total=total;
    output1_LCC.Percent_Const=total_const./total;
    output1_LCC.Percent_OM=Present_OM./total;
    
    %% Inventory Analysis
    if i==1
        output_A_B.(current_system).LCC=output1_LCC;
    else
        output_A_B.(current_system).LCC=[output_A_B.(current_system).LCC; output1_LCC];
    end
    
    %% Impact Assessment
    [traci_con, CML_con,CED_con, traci_op, CML_op,CED_op, traci_avoid, CML_avoid,CED_avoid, IMPACT_DE] = Impact_Assessment (INV_CON, INV_OP, DEmis, E_offset_kWh, V_treated);
    %Impacts are per m3
    
    %% Impact by Sources
    GWP_total=0;
    all_con = [traci_con CML_con CED_con];
    all_op = [traci_op CML_op CED_op];
    all_avoid = [traci_avoid CML_avoid CED_avoid];
    all_con.Properties.RowNames = {};
    all_op.Properties.RowNames = {};
    all_avoid.Properties.RowNames = {};
    all_VarNames = all_con.Properties.VariableNames;
    
    if i==1
        output_A_B.(current_system).Concrete=all_con(2,:);
        output_A_B.(current_system).Steel=array2table(all_con{3,:}+all_con{9,:},'VariableNames',all_VarNames);
        output_A_B.(current_system).Excavation=all_con(19,:);
        output_A_B.(current_system).Construction_Misc=array2table(sum(all_con{4:8,:})+sum(all_con{10:18,:})+sum(all_con{20:38,:}),'VariableNames',all_VarNames);
        output_A_B.(current_system).Lift_Pumping=array2table(Power_pct(1)*sum(all_op{10:16,:}),'VariableNames',all_VarNames);
        output_A_B.(current_system).Gas_Sparging=array2table(Power_pct(2)*sum(all_op{10:16,:}),'VariableNames',all_VarNames);
        output_A_B.(current_system).Gravity_Belt_Thickeners=array2table(Power_pct(3)*sum(all_op{10:16,:}),'VariableNames',all_VarNames);
        output_A_B.(current_system).Operation_Misc=array2table(sum(all_op{1:9,:}),'VariableNames',all_VarNames);
    else
        output_A_B.(current_system).Concrete=[output_A_B.(current_system).Concrete; all_con(2,:)];
        output_A_B.(current_system).Steel=[output_A_B.(current_system).Steel;array2table(all_con{3,:}+all_con{9,:},'VariableNames',all_VarNames)];
        output_A_B.(current_system).Excavation=[output_A_B.(current_system).Excavation;all_con(19,:)];
        output_A_B.(current_system).Construction_Misc=[output_A_B.(current_system).Construction_Misc;array2table(sum(all_con{4:8,:})+sum(all_con{10:18,:})+sum(all_con{20:38,:}),'VariableNames',all_VarNames)];
        output_A_B.(current_system).Lift_Pumping=[output_A_B.(current_system).Lift_Pumping;array2table(Power_pct(1)*sum(all_op{10:16,:}),'VariableNames',all_VarNames)];
        output_A_B.(current_system).Gas_Sparging=[output_A_B.(current_system).Gas_Sparging;array2table(Power_pct(2)*sum(all_op{10:16,:}),'VariableNames',all_VarNames)];
        output_A_B.(current_system).Gravity_Belt_Thickeners=[output_A_B.(current_system).Gravity_Belt_Thickeners;array2table(Power_pct(3)*sum(all_op{10:16,:}),'VariableNames',all_VarNames)];
        output_A_B.(current_system).Operation_Misc=[output_A_B.(current_system).Operation_Misc;array2table(sum(all_op{1:9,:}),'VariableNames',all_VarNames)];
    end
    
    % Direct emission
    if i==1
        output_A_B.(current_system).Direct_Emissions=array2table(sum(IMPACT_DE{1:8,:}),'VariableNames',all_VarNames);
    else
        output_A_B.(current_system).Direct_Emissions=[output_A_B.(current_system).Direct_Emissions;array2table(sum(IMPACT_DE{1:8,:}),'VariableNames',all_VarNames)];
    end
    
    % Avoided Energy
    if i==1
        output_A_B.(current_system).Avoided_Energy=array2table(sum(all_avoid{1:7,:}),'VariableNames',all_VarNames);
    else
        output_A_B.(current_system).Avoided_Energy=[output_A_B.(current_system).Avoided_Energy;array2table(sum(all_avoid{1:7,:}),'VariableNames',all_VarNames)];
    end
    
    Energy_IO(i,1) = E_input_kWh;
    Energy_IO(i,2) = E_offset_kWh;
    Energy_IO(i,3) = E_input_kWh-E_offset_kWh;
    Energy_IO(i,4) = E_offset_kWh/E_input_kWh;
    
    % %For Comparison of LCC to GWP
    % GWP_total=GWP_total+Concrete_pct(2)+Steel_pct(2)+Excavation_pct(2)+CON_MISC_pct(2)+Pumping_PERM_pct(2)+Pumping_IR_pct(2)+...
    %     Pumping_CHEM_pct(2)+Pumping_GBT_pct(2)+OP_MISC_pct(2)+DE_pct(2)+Avoided_pct(2);
    % LCC_and_GWP=[output1_cost(1,length(output1_cost)),GWP_total];
    % if i==1
    %     output1.(current_system).LCC_and_GWP=LCC_and_GWP;
    % else
    %     output1.(current_system).LCC_and_GWP=[output1.(current_system).LCC_and_GWP LCC_and_GWP];
    % end
    
    %% Clear so we don't have any mistakes in the data
    clear {'INV_CON', 'INV_OP', 'DEmis', 'Power_pct', 'E_input_kWh', 'E_offset_kWh', 'V_treated','output1_cost','IMPACT_CON', 'IMPACT_OP', 'IMPACT_DE', 'IMPACT_avoided','CON_TOT','OP_TOT','IMPACT_TOT','output1_LCC'};
    if i==trials
        names=fieldnames(output_A_B.(current_system));
        for j=1:length(names)
            if isstruct(output_A_B.(current_system).(names{j}))==1
                output_A_B.(current_system).(names{j})=struct2table(output_A_B.(current_system).(names{j}));
                %         else
                %             output1.(current_system).(names{j}).Properties.VariableNames=LCA_impacts;
            end
        end
    end
    waitbar(i/trials)
end

%% Calculations
designs=fieldnames(output_A_B);
output_A_B.LCC_LCA=table();
output_A_B.LCA=table();
for i=1:length(designs)
    categories=fieldnames(output_A_B.(designs{i}));
    for j=1:length(categories)
        if j==1
            output_A_B.LCC_LCA=[output_A_B.LCC_LCA; varfun(@mean,output_A_B.(designs{i}).(categories{j})) varfun(@std,output_A_B.(designs{i}).(categories{j})) varfun(@median,output_A_B.(designs{i}).(categories{j}))]; %Calculates average and standard dev. and median from all trials
            output_A_B.(designs{i}).LCA_Avg=table();
        else
            output_A_B.(designs{i}).LCA_Avg=[output_A_B.(designs{i}).LCA_Avg; varfun(@mean,output_A_B.(designs{i}).(categories{j})) varfun(@var,output_A_B.(designs{i}).(categories{j})) varfun(@median,output_A_B.(designs{i}).(categories{j}))]; %Calculates average, variance, and median of impacts from all trials
        end
    end
    output_A_B.(designs{i}).LCA_Avg.Properties.RowNames=categories(2:length(categories));
    output_A_B.(designs{i}).LCA_sum=table();
    output_A_B.(designs{i}).LCA_sum=varfun(@sum,output_A_B.(designs{i}).LCA_Avg); %Finds the sum of means and variances for each impact category for a given design
    output_A_B.(designs{i}).LCA_sum{1,10:18}=output_A_B.(designs{i}).LCA_sum{1,10:18}.^2; %Squares variance to calculate standard dev.
    output_A_B.(designs{i}).LCA_Const=varfun(@sum,output_A_B.(designs{i}).LCA_Avg(1:4,:));
    output_A_B.(designs{i}).LCA_Const.Properties.VariableNames=Const_names;
    output_A_B.(designs{i}).LCA_OM=varfun(@sum,output_A_B.(designs{i}).LCA_Avg(5:9,:));
    output_A_B.(designs{i}).LCA_OM.Properties.VariableNames=OM_names;
    output_A_B.LCA=[output_A_B.LCA; output_A_B.(designs{i}).LCA_sum output_A_B.(designs{i}).LCA_Const output_A_B.(designs{i}).LCA_OM];
end
output_A_B.LCC_LCA.Properties.RowNames=designs;
output_A_B.LCC_LCA=[output_A_B.LCC_LCA output_A_B.LCA];
close(wait)
toc