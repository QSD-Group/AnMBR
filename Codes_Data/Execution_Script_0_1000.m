%clc; 
clear
load results_8-18-15
tic
output=struct();
LCA_impacts = {'Ozone','Global_Warming','Smog','Acidification','Eutrophication','Carbon_Emis','NonCarbon_Emis','Respiratory_Disease','Ecotoxicity'};
Const_names = {'CON_sum_mean_OzoneDepletion','CON_sum_mean_GlobalWarming','CON_sum_mean_Smog','CON_sum_mean_Acidification','CON_sum_mean_Eutrophication','CON_sum_mean_Carcinogenics','CON_sum_mean_NonCarcinogenics','CON_sum_mean_RespiratoryEffects','CON_sum_mean_Ecotoxicity','CON_sum_mean_FossilFuelDepletion','CON_sum_mean_AbioticDepletion','CON_sum_mean_AbioticDepletion_fossilFuels_','CON_sum_mean_GlobalWarming_GWP100a_','CON_sum_mean_OzoneLayerDepletion_ODP_','CON_sum_mean_HumanToxicity','CON_sum_mean_FreshWaterAquaticEcotox_','CON_sum_mean_MarineAquaticEcotoxicity','CON_sum_mean_TerrestrialEcotoxicity','CON_sum_mean_PhotochemicalOxidation','CON_sum_mean_Acidification_CML','CON_sum_mean_Eutrophication_CML','CON_sum_mean_NonRenewable_Fossil','CON_sum_mean_Non_renewable_Nuclear','CON_sum_mean_Non_renewable_Biomass','CON_sum_mean_Renewable_Biomass','CON_sum_mean_Renewable_Wind_Solar_Geothe','CON_sum_mean_Renewable_Water','CON_sum_var_OzoneDepletion','CON_sum_var_GlobalWarming','CON_sum_var_Smog','CON_sum_var_Acidification','CON_sum_var_Eutrophication','CON_sum_var_Carcinogenics','CON_sum_var_NonCarcinogenics','CON_sum_var_RespiratoryEffects','CON_sum_var_Ecotoxicity','CON_sum_var_FossilFuelDepletion','CON_sum_var_AbioticDepletion','CON_sum_var_AbioticDepletion_fossilFuels_','CON_sum_var_GlobalWarming_GWP100a_','CON_sum_var_OzoneLayerDepletion_ODP_','CON_sum_var_HumanToxicity','CON_sum_var_FreshWaterAquaticEcotox_','CON_sum_var_MarineAquaticEcotoxicity','CON_sum_var_TerrestrialEcotoxicity','CON_sum_var_PhotochemicalOxidation','CON_sum_var_Acidification_CML','CON_sum_var_Eutrophication_CML','CON_sum_var_NonRenewable_Fossil','CON_sum_var_Non_renewable_Nuclear','CON_sum_var_Non_renewable_Biomass','CON_sum_var_Renewable_Biomass','CON_sum_var_Renewable_Wind_Solar_Geothe','CON_sum_var_Renewable_Water','CON_sum_median_OzoneDepletion','CON_sum_median_GlobalWarming','CON_sum_median_Smog','CON_sum_median_Acidification','CON_sum_median_Eutrophication','CON_sum_median_Carcinogenics','CON_sum_median_NonCarcinogenics','CON_sum_median_RespiratoryEffects','CON_sum_median_Ecotoxicity','CON_sum_median_FossilFuelDepletion','CON_sum_median_AbioticDepletion','CON_sum_median_AbioticDepletion_fossilFuels_','CON_sum_median_GlobalWarming_GWP100a_','CON_sum_median_OzoneLayerDepletion_ODP_','CON_sum_median_HumanToxicity','CON_sum_median_FreshWaterAquaticEcotox_','CON_sum_median_MarineAquaticEcotoxicity','CON_sum_median_TerrestrialEcotoxicity','CON_sum_median_PhotochemicalOxidation','CON_sum_median_Acidification_CML','CON_sum_median_Eutrophication_CML','CON_sum_median_NonRenewable_Fossil','CON_sum_median_Non_renewable_Nuclear','CON_sum_median_Non_renewable_Biomass','CON_sum_median_Renewable_Biomass','CON_sum_median_Renewable_Wind_Solar_Geothe','CON_sum_median_Renewable_Water'};
OM_names = {'OP_sum_mean_OzoneDepletion','OP_sum_mean_GlobalWarming','OP_sum_mean_Smog','OP_sum_mean_Acidification','OP_sum_mean_Eutrophication','OP_sum_mean_Carcinogenics','OP_sum_mean_NonCarcinogenics','OP_sum_mean_RespiratoryEffects','OP_sum_mean_Ecotoxicity','OP_sum_mean_FossilFuelDepletion','OP_sum_mean_AbioticDepletion','OP_sum_mean_AbioticDepletion_fossilFuels_','OP_sum_mean_GlobalWarming_GWP100a_','OP_sum_mean_OzoneLayerDepletion_ODP_','OP_sum_mean_HumanToxicity','OP_sum_mean_FreshWaterAquaticEcotox_','OP_sum_mean_MarineAquaticEcotoxicity','OP_sum_mean_TerrestrialEcotoxicity','OP_sum_mean_PhotochemicalOxidation','OP_sum_mean_Acidification_CML','OP_sum_mean_Eutrophication_CML','OP_sum_mean_NonRenewable_Fossil','OP_sum_mean_Non_renewable_Nuclear','OP_sum_mean_Non_renewable_Biomass','OP_sum_mean_Renewable_Biomass','OP_sum_mean_Renewable_Wind_Solar_Geothe','OP_sum_mean_Renewable_Water','OP_sum_var_OzoneDepletion','OP_sum_var_GlobalWarming','OP_sum_var_Smog','OP_sum_var_Acidification','OP_sum_var_Eutrophication','OP_sum_var_Carcinogenics','OP_sum_var_NonCarcinogenics','OP_sum_var_RespiratoryEffects','OP_sum_var_Ecotoxicity','OP_sum_var_FossilFuelDepletion','OP_sum_var_AbioticDepletion','OP_sum_var_AbioticDepletion_fossilFuels_','OP_sum_var_GlobalWarming_GWP100a_','OP_sum_var_OzoneLayerDepletion_ODP_','OP_sum_var_HumanToxicity','OP_sum_var_FreshWaterAquaticEcotox_','OP_sum_var_MarineAquaticEcotoxicity','OP_sum_var_TerrestrialEcotoxicity','OP_sum_var_PhotochemicalOxidation','OP_sum_var_Acidification_CML','OP_sum_var_Eutrophication_CML','OP_sum_var_NonRenewable_Fossil','OP_sum_var_Non_renewable_Nuclear','OP_sum_var_Non_renewable_Biomass','OP_sum_var_Renewable_Biomass','OP_sum_var_Renewable_Wind_Solar_Geothe','OP_sum_var_Renewable_Water','OP_sum_median_OzoneDepletion','OP_sum_median_GlobalWarming','OP_sum_median_Smog','OP_sum_median_Acidification','OP_sum_median_Eutrophication','OP_sum_median_Carcinogenics','OP_sum_median_NonCarcinogenics','OP_sum_median_RespiratoryEffects','OP_sum_median_Ecotoxicity','OP_sum_median_FossilFuelDepletion','OP_sum_median_AbioticDepletion','OP_sum_median_AbioticDepletion_fossilFuels_','OP_sum_median_GlobalWarming_GWP100a_','OP_sum_median_OzoneLayerDepletion_ODP_','OP_sum_median_HumanToxicity','OP_sum_median_FreshWaterAquaticEcotox_','OP_sum_median_MarineAquaticEcotoxicity','OP_sum_median_TerrestrialEcotoxicity','OP_sum_median_PhotochemicalOxidation','OP_sum_median_Acidification_CML','OP_sum_median_Eutrophication_CML','OP_sum_median_NonRenewable_Fossil','OP_sum_median_Non_renewable_Nuclear','OP_sum_median_Non_renewable_Biomass','OP_sum_median_Renewable_Biomass','OP_sum_median_Renewable_Wind_Solar_Geothe','OP_sum_median_Renewable_Water'};
CML_names = {'AbioticDepletion','AbioticDepletion_fossilFuels_','GlobalWarming_GWP100a_','OzoneLayerDepletion_ODP_','HumanToxicity','FreshWaterAquaticEcotox_','MarineAquaticEcotoxicity','TerrestrialEcotoxicity','PhotochemicalOxidation','Acidification_CML','Eutrophication_CML'};
CED_names = {'NonRenewable_Fossil','Non_renewable_Nuclear','Non_renewable_Biomass','Renewable_Biomass','Renewable_Wind_Solar_Geothe','Renewable_Water'};
CED_con=array2table(zeros(38,6), 'VariableNames',CED_names);
CED_op=array2table(zeros(16,6), 'VariableNames',CED_names);
CED_avoid=array2table(zeros(7,6), 'VariableNames',CED_names);
CML_con=array2table(zeros(38,11), 'VariableNames',CML_names);
CML_op=array2table(zeros(16,11), 'VariableNames',CML_names);
CML_avoid=array2table(zeros(7,11), 'VariableNames',CML_names);
Designs_matrix;
%% Variables/Parameters
%Decision variables
% trials=3000;
% 
% % J = 8.5;
% % TMP_high = 25;
% % TMP_low = 21;
% % OL_AF = 5;
% % HL_AF = 10;
% % OL_AER = 2;
% % HL_AER = 2;
% % HRT = 10;
% % IRR = 4;
% % interest = 0.08;
% % GAC = 100;
% % CFV = 3;
% % SGD = 1.7;
% % Sparg_freq = .75;
% % SRT = 35;
% % Biogas_CH4 = 0.6;
% % 
% % Q_mgd = 30;
% % S_SO = 300;
% % X_SO = 100;
% % dummy = 100;
% 
% J = lhs_triangle(5,12,17,trials); % [L/m^2-hr] (Triang. dist.)
% TMP_high = lhs_triangle(0.5656,2.5,5.37,trials); % [psi] (Triang. dist.)
% TMP_low = lhs_triangle(0.5656,2.5,5.37,trials); % [psi] (Triang. dist.)
% membrane_life = lhs_triangle(5,10,15,trials); % [years] (Triang. dist)
% membrane_cost = lhsu(6,10,trials); % [$/sf], based on Ashland memo
% OL_AF = lhsu(0.2,8,trials); % [g-COD/L-d] or [kg/m^3-day] (Uniform dist.)
% HL_AF = lhsu(2,6,trials); % [m/hr] (Uniform dist.)
% OL_AER = lhsu(0.5,4,trials); % [g-COD/L-d] or [kg/m^3-day] (Uniform dist.)
% HL_AER = lhsu(0.11,0.44,trials); % [m/hr] (Uniform dist.)
% HRT = lhsu(8,16,trials); % [hr] (Uniform dist.)
% HRT_for_GAC = lhsu(2.2,3.3,trials); % [hr] (Uniform dist.)
% IRR = lhsu(0.5,4,trials); %(Uniform dist.)
% interest=lhsu(0.06,0.10,trials); % [%] (Uniform dist.)
% GAC = lhsu(187.8,225,trials); %[g/L] Concentration used in Kim et al, 2011
% CFV = lhsu(0.4,2,trials); %[m/s](Uniform dist.)
% SGD = lhsu(0.05,1.2,trials); %[m3/m2-hr](Uniform dist.)
% Sparg_freq = lhsu(0.5,1,trials); %[%] (Uniform dist.)
% yield = lhsu(0.02,0.08,trials); %[d] (Uniform dist.)
% Biogas_CH4 = lhsu(0.5,0.7,trials); %Dissolved CH4 (Unifrom dist.) [From Smith et al., 2011 and Kim et al., 2011, saying that 50% and 70% of methane (respectively) is recovered in gas]
% Upflow_vel_for_GAC = lhs_triangle(6,8,10,trials); %[m/hr] Upflow velocity for GAC bed expansion, based on Aslan, 2014
% 
% %Uncertainty parameters
% Q_mgd = latin_hs(20,5,trials,1); % [mgd] (Normal dist., 25% std. dev)
% S_SO = latin_hs(300,30,trials,1); % [mg-COD/L] or [g/m^3] (Normal dist., 10% std. dev)
% X_SO = latin_hs(100,10,trials,1); % [mg-COD/L] or [g/m^3] (Normal dist., 10% std. dev)
% dummy = latin_hs(100,10,trials,1); % Dummy parameter for sensitivity analysis
% 
% %Table Initialization
% vars=table(J,TMP_high,TMP_low,OL_AF,HL_AF,OL_AER,HL_AER,HRT,HRT_for_GAC, Upflow_vel_for_GAC,IRR,interest,GAC,CFV,SGD,Sparg_freq,yield,membrane_life,membrane_cost,Biogas_CH4,Q_mgd,S_SO,X_SO,dummy);
% varnames=vars.Properties.VariableNames;
% vars=varfun(@abs,vars);
% vars.Properties.VariableNames=varnames;

%% Initiation of decision loops
wait=waitbar(0,'');
% for k = 1:height(all_designs)
%     step_A = all_designs.step_A{k};
%     step_A1 = all_designs.step_A1{k};
%     step_B = all_designs.step_B{k};
%     step_C = all_designs.step_C{k};
%     step_D = all_designs.step_D{k};
%     step_H = all_designs.step_H{k};
%     %step_I = all_designs.step_I{k};
for i = 1:1000
    for ia = 1:2
        for ia1 = 1:3
            for ib = 1:2
                for ic = 1:3
                    for id=1:6
                        for ih=1:2
                            %for ii = 1:4
if ib==1 && ic==3 %Prevents the 'Submerged MT' phenotype
    break
elseif ib==2 && ic==1 %Prevents the 'Cross-flow HF' phenotype
    break
elseif ia==1 && ia1==3 %Prevents the CSTR AER phenotype
    break 
elseif ic==1 && id==3 || ic==1 && id==4 %Makes sure hollow fiber is only made of plastic
    break
elseif ic==2 && id==3 %Prevents flat sheet made of ceramic phenotype
    break
elseif ic==3 && id==4 %Prevents multi-tube made of steel phenotype
    break
elseif ia1==2 && ib==2 %Prevents GAC in cross-flow membranes
    break
end

if ia==1 %Reactor type
    step_A='CSTR';
else
    step_A='AF';
end

if ia1==1 %Additional options
    step_A1='none';
elseif ia1==2
    step_A1='GAC';
else
    step_A1='AER';
end

if ib==1 %Membrane configuration
    step_B='Submerged';
elseif ib==2
    step_B='Crossflow';
end

if ic==1 %Membrane type
    step_C='HF'; %Hollow fiber
elseif ic==2
    step_C='FS'; %Flat sheet
else
    step_C='MT'; %Multi-tube
end

if id==1 %Membrane material
    step_D='PET';
elseif id==2
    step_D='PTFE';
elseif id==3
    step_D='ceramic';
elseif id==4
    step_D='sinter_steel';
elseif id==5
    step_D='PES';
else
    step_D='PVDF';
end

if ih==1 %use DM or not
    step_H='DM';
else
    step_H='no_DM';
end

% if ii==1 %CHP type
%     step_I='IC'; %Internal Combustion
% elseif ii==2
%     step_I='CG'; %Combustion Gas
% elseif ii==3
%     step_I='micro'; %Microturbine
% else
%     step_I='FC'; %Fuel cell
% end
current_system=horzcat(step_A,'_',step_A1,'_',step_B,'_',step_C,'_',step_D,'_',step_H); %,'_',step_I);
%% Assign values
UParameter = [vars.Q_mgd(i), vars.S_SO(i), vars.X_SO(i),vars.CFV(i),vars.SGD(i),vars.Sparg_freq(i),vars.yield(i),vars.Biogas_CH4(i), vars.Upflow_vel_for_GAC(i),vars.membrane_life(i),vars.membrane_cost(i)];

if ia==1 && ib==1 && ia1==1 %CSTR_none_Submerged
    DVariable = [vars.J(i), vars.TMP_high(i), vars.IRR(i), vars.HRT(i), vars.HRT_for_GAC(i)];
elseif ia==1 && ib==2 && ia1==1 %CSTR_none_Crossflow
    DVariable = [vars.J(i), vars.TMP_low(i), vars.IRR(i), vars.HRT(i), vars.HRT_for_GAC(i)];
elseif ia==1 && ib==1 && ia1==2 %CSTR_GAC_Submerged
    DVariable = [vars.J(i), vars.TMP_high(i), vars.IRR(i), vars.HRT(i), vars.HRT_for_GAC(i),vars.GAC(i)];
% elseif ia==1 && ib==2 && ia1==2 %CSTR_GAC_Crossflow
%     DVariable = [vars.J(i), vars.TMP_low(i), vars.IRR(i), vars.HRT(i), vars.GAC(i),vars.HRT_for_GAC(i)];
elseif ia==2 && ib==2 && ia1==1 %AF_none_Crossflow
    DVariable = [vars.J(i), vars.TMP_low(i), vars.OL_AF(i), vars.HL_AF(i), vars.IRR(i)];
% elseif ia==2 && ia1==2 %AF_GAC_Crossflow
%     DVariable = [vars.J(i), vars.TMP_low(i), vars.OL_AF(i), vars.HL_AF(i), vars.IRR(i), vars.GAC(i)];
elseif ia==2 && ib==2 && ia1==3 %AF_AER_Crossflow
    DVariable = [vars.J(i), vars.TMP_low(i), vars.OL_AF(i), vars.HL_AF(i), vars.IRR(i), vars.OL_AER(i), vars.HL_AER(i)];
elseif ia==2 && ib==1 && ia1==1 % AF_none_Submerged
    DVariable = [vars.J(i), vars.TMP_high(i), vars.OL_AF(i), vars.HL_AF(i), vars.IRR(i)];
elseif ia==2 && ib==1 && ia1==2 % AF_GAC_Submerged
    DVariable = [vars.J(i), vars.TMP_high(i), vars.OL_AF(i), vars.HL_AF(i), vars.IRR(i),vars.HRT_for_GAC(i),vars.GAC(i)];
elseif ia==2 && ib==1 && ia1==3 % AF_AER_Submerged
    DVariable = [vars.J(i), vars.TMP_high(i), vars.OL_AF(i), vars.HL_AF(i), vars.IRR(i),vars.OL_AER(i), vars.HL_AER(i)];
else 
    break
end
%% Inventory Analysis

[INV_CON, INV_OP, DEmis, Power_pct, E_input_kWh, E_offset_kWh, V_treated, Output_cost] = LCI_anMBR_HandS_Best(ia,ia1,ib,ic,id,ih, DVariable, UParameter,interest(i));

if i==1
    output.(current_system).LCC=Output_cost;
else
    output.(current_system).LCC=[output.(current_system).LCC; Output_cost];
end

%% Impact Assessment
[traci_con, CED_con, traci_op, CED_op, traci_avoid, CED_avoid, IMPACT_DE] = Impact_Assessment (INV_CON, INV_OP, DEmis, E_offset_kWh, V_treated);
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
    output.(current_system).Concrete=all_con(2,:);
    output.(current_system).Steel=array2table(all_con{3,:}+all_con{9,:},'VariableNames',all_VarNames);
    output.(current_system).Excavation=all_con(19,:);
    output.(current_system).Construction_Misc=array2table(sum(all_con{4:8,:})+sum(all_con{10:18,:})+sum(all_con{20:38,:}),'VariableNames',all_VarNames);
    output.(current_system).Permeate_Pumping=array2table(Power_pct(1)*sum(all_op{10:16,:}),'VariableNames',all_VarNames);
    output.(current_system).Recirculation_Pumping=array2table(Power_pct(2)*sum(all_op{10:16,:}),'VariableNames',all_VarNames);
    output.(current_system).Chemical_Pumping=array2table(Power_pct(3)*sum(all_op{10:16,:}),'VariableNames',all_VarNames);
    output.(current_system).Gravity_Belt_Thickeners=array2table(Power_pct(4)*sum(all_op{10:16,:}),'VariableNames',all_VarNames);
    output.(current_system).Operation_Misc=array2table(sum(all_op{1:9,:}),'VariableNames',all_VarNames);
else
    output.(current_system).Concrete=[output.(current_system).Concrete; all_con(2,:)];
    output.(current_system).Steel=[output.(current_system).Steel;array2table(all_con{3,:}+all_con{9,:},'VariableNames',all_VarNames)];
    output.(current_system).Excavation=[output.(current_system).Excavation;all_con(19,:)];
    output.(current_system).Construction_Misc=[output.(current_system).Construction_Misc;array2table(sum(all_con{4:8,:})+sum(all_con{10:18,:})+sum(all_con{20:38,:}),'VariableNames',all_VarNames)];
    output.(current_system).Permeate_Pumping=[output.(current_system).Permeate_Pumping;array2table(Power_pct(1)*sum(all_op{10:16,:}),'VariableNames',all_VarNames)];
    output.(current_system).Recirculation_Pumping=[output.(current_system).Recirculation_Pumping;array2table(Power_pct(2)*sum(all_op{10:16,:}),'VariableNames',all_VarNames)];
    output.(current_system).Chemical_Pumping=[output.(current_system).Chemical_Pumping;array2table(Power_pct(3)*sum(all_op{10:16,:}),'VariableNames',all_VarNames)];
    output.(current_system).Gravity_Belt_Thickeners=[output.(current_system).Gravity_Belt_Thickeners;array2table(Power_pct(4)*sum(all_op{10:16,:}),'VariableNames',all_VarNames)];
    output.(current_system).Operation_Misc=[output.(current_system).Operation_Misc;array2table(sum(all_op{1:9,:}),'VariableNames',all_VarNames)];
end

%CSTR_none_Submerged
if ia==1 && ia1==1 && ib==1 %&& ic==1||ia==1 && ia1==1 && ib==1 && ic==2 
    if i==1
        output.(current_system).Gas_Sparging=array2table(Power_pct(5)*sum(all_op{10:16,:}),'VariableNames',all_VarNames);
        output.(current_system).Lift_Pumping=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);%table([0],[0],[0],[0],[0],[0],[0],[0],[0]);
        output.(current_system).Retentate_Pumping=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
    else
        output.(current_system).Gas_Sparging=[output.(current_system).Gas_Sparging;array2table(Power_pct(5)*sum(all_op{10:16,:}),'VariableNames',all_VarNames)];
    end
%CSTR_GAC_Submerged
elseif ia==1 && ia1==2 && ib==1 %&& ic==1||ia==1 && ia1==2 && ib==1 && ic==2 
    if i==1
        output.(current_system).Gas_Sparging=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
        output.(current_system).Lift_Pumping=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
        output.(current_system).Retentate_Pumping=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
    end
%CSTR_none_Cross-flow
elseif ia==1 && ia1==1 && ib==2 %&& ic==3||ia==1 && ia1==1 && ib==2 && ic==2
    if i==1
        output.(current_system).Gas_Sparging=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
        output.(current_system).Lift_Pumping=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
        output.(current_system).Retentate_Pumping=array2table(Power_pct(5)*sum(all_op{10:16,:}),'VariableNames',all_VarNames);
    else
        output.(current_system).Retentate_Pumping=[output.(current_system).Retentate_Pumping;array2table(Power_pct(5)*sum(all_op{10:16,:}),'VariableNames',all_VarNames)];
    end
%AF_none_Crossflow OR AF_AeF_Crossflow
elseif ia==2 && ib==2 %ia==2 && ia1==1 && ib==2 && ic==3||ia==2 && ia1==1 && ib==2 && ic==2||ia==2 && ia1==3 && ib==2 && ic==3||ia==2 && ia1==3 && ib==2 && ic==2
    if i==1
        output.(current_system).Gas_Sparging=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
        output.(current_system).Lift_Pumping=array2table(Power_pct(5)*sum(all_op{10:16,:}),'VariableNames',all_VarNames);
        output.(current_system).Retentate_Pumping=array2table(Power_pct(6)*sum(all_op{10:16,:}),'VariableNames',all_VarNames);
    else
        output.(current_system).Lift_Pumping=[output.(current_system).Lift_Pumping;array2table(Power_pct(5)*sum(all_op{10:16,:}),'VariableNames',all_VarNames)];
        output.(current_system).Retentate_Pumping=[output.(current_system).Retentate_Pumping;array2table(Power_pct(6)*sum(all_op{10:16,:}),'VariableNames',all_VarNames)];
    end
%AF_none_Submerged OR AF_AeF_Submerged
elseif ia==2 && ia1~=2 && ib==1 %&& ic==1||ia==2 && ia1==1 && ib==1 && ic==2
    if i==1
        output.(current_system).Gas_Sparging=array2table(Power_pct(5)*sum(all_op{10:16,:}),'VariableNames',all_VarNames);
        output.(current_system).Lift_Pumping=array2table(Power_pct(6)*sum(all_op{10:16,:}),'VariableNames',all_VarNames);
        output.(current_system).Retentate_Pumping=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
    else
        output.(current_system).Gas_Sparging=[output.(current_system).Gas_Sparging;array2table(Power_pct(5)*sum(all_op{10:16,:}),'VariableNames',all_VarNames)];
        output.(current_system).Lift_Pumping=[output.(current_system).Lift_Pumping;array2table(Power_pct(6)*sum(all_op{10:16,:}),'VariableNames',all_VarNames)];
    end
%AF_GAC_Submerged
elseif ia==2 && ia1==2 && ib==1 %&& ic==1||ia==2 && ia1==2 && ib==1 && ic==2||ia==2 && ia1==3 && ib==1 && ic==1||ia==2 && ia1==3 && ib==1 && ic==2
    if i==1
        output.(current_system).Gas_Sparging=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
        output.(current_system).Lift_Pumping=array2table(Power_pct(5)*sum(all_op{10:16,:}),'VariableNames',all_VarNames);
        output.(current_system).Retentate_Pumping=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
    else
        output.(current_system).Lift_Pumping=[output.(current_system).Lift_Pumping;array2table(Power_pct(5)*sum(all_op{10:16,:}),'VariableNames',all_VarNames)];
    end
end

if ih==1
    if i==1
        output.(current_system).Degassing_Memb_Op=array2table(Power_pct(length(Power_pct))*sum(all_op{10:16,:}),'VariableNames',all_VarNames);
    else
        output.(current_system).Degassing_Memb_Op=[output.(current_system).Degassing_Memb_Op;array2table(Power_pct(length(Power_pct))*sum(all_op{10:16,:}),'VariableNames',all_VarNames)];
    end
else
    if i==1
        output.(current_system).Degassing_Memb_Op=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
    end
end

% Direct emission
if i==1
    output.(current_system).Direct_Emissions=array2table(sum(IMPACT_DE{1:8,:}),'VariableNames',all_VarNames);
else
    output.(current_system).Direct_Emissions=[output.(current_system).Direct_Emissions;array2table(sum(IMPACT_DE{1:8,:}),'VariableNames',all_VarNames)];
end

% Avoided Energy
if i==1
    output.(current_system).Avoided_Energy=array2table(sum(all_avoid{1:7,:}),'VariableNames',all_VarNames);
else
    output.(current_system).Avoided_Energy=[output.(current_system).Avoided_Energy;array2table(sum(all_avoid{1:7,:}),'VariableNames',all_VarNames)];
end

% %For Comparison of LCC to GWP
% GWP_total=GWP_total+Concrete_pct(2)+Steel_pct(2)+Excavation_pct(2)+CON_MISC_pct(2)+Pumping_PERM_pct(2)+Pumping_IR_pct(2)+...
%     Pumping_CHEM_pct(2)+Pumping_GBT_pct(2)+OP_MISC_pct(2)+DE_pct(2)+Avoided_pct(2);
% LCC_and_GWP=[Output_cost(1,length(Output_cost)),GWP_total];
% if i==1
%     output.(current_system).LCC_and_GWP=LCC_and_GWP;
% else
%     output.(current_system).LCC_and_GWP=[output.(current_system).LCC_and_GWP LCC_and_GWP];
% end

%% Clear so we don't have any mistakes in the data
clear {'INV_CON', 'INV_OP', 'DEmis', 'Power_pct', 'E_input_kWh', 'E_offset_kWh', 'V_treated','Output_cost','IMPACT_CON', 'IMPACT_OP', 'IMPACT_DE', 'IMPACT_avoided','CON_TOT','OP_TOT','IMPACT_TOT','Output_LCC'};
if i==1000
    names=fieldnames(output.(current_system));
    for j=1:length(names)
        if isstruct(output.(current_system).(names{j}))==1
            output.(current_system).(names{j})=struct2table(output.(current_system).(names{j}));
%         else
%             output.(current_system).(names{j}).Properties.VariableNames=LCA_impacts;
        end
    end
end
                        end
                    end
                end
            end
        end
    end
    waitbar(i/1000)
end

%% Calculations
designs=fieldnames(output);
output.LCC_LCA=table();
output.LCA=table();
for i=1:length(designs)
    categories=fieldnames(output.(designs{i}));
    for j=1:length(categories)
        if j==1
            output.LCC_LCA=[output.LCC_LCA; varfun(@mean,output.(designs{i}).(categories{j})) varfun(@std,output.(designs{i}).(categories{j})) varfun(@median,output.(designs{i}).(categories{j}))]; %Calculates average and standard dev. and median from all trials
            output.(designs{i}).LCA_Avg=table();
        else
            output.(designs{i}).LCA_Avg=[output.(designs{i}).LCA_Avg; varfun(@mean,output.(designs{i}).(categories{j})) varfun(@var,output.(designs{i}).(categories{j})) varfun(@median,output.(designs{i}).(categories{j}))]; %Calculates average, variance, and median of impacts from all trials
        end
    end
    output.(designs{i}).LCA_Avg.Properties.RowNames=categories(2:length(categories));
    output.(designs{i}).LCA_sum=table();
    output.(designs{i}).LCA_sum=varfun(@sum,output.(designs{i}).LCA_Avg); %Finds the sum of means and variances for each impact category for a given design
    output.(designs{i}).LCA_sum{1,10:18}=output.(designs{i}).LCA_sum{1,10:18}.^2; %Squares variance to calculate standard dev.
    output.(designs{i}).LCA_Const=varfun(@sum,output.(designs{i}).LCA_Avg(1:4,:));
    output.(designs{i}).LCA_Const.Properties.VariableNames=Const_names;
    output.(designs{i}).LCA_OM=varfun(@sum,output.(designs{i}).LCA_Avg(5:14,:));
    output.(designs{i}).LCA_OM.Properties.VariableNames=OM_names;
    output.LCA=[output.LCA; output.(designs{i}).LCA_sum output.(designs{i}).LCA_Const output.(designs{i}).LCA_OM];
end
output.LCC_LCA.Properties.RowNames=designs;
output.LCC_LCA=[output.LCC_LCA output.LCA];
%one_more = output.LCC_LCA{:,[64 121:129 142:147]};
close(wait)
toc
