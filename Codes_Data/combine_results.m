output=struct();
for ia = 1:2
    for ia1 = 1:3
        for ib = 1:2
            for ic = 1:3
                for id=1:6
                    for ih=1:2
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
                        
                        current_system=horzcat(step_A,'_',step_A1,'_',step_B,'_',step_C,'_',step_D,'_',step_H);
                        output.(current_system).LCC=[output_1000.(current_system).LCC; output_2000.(current_system).LCC; output_3000.(current_system).LCC];
                        output.(current_system).Concrete=[output_1000.(current_system).Concrete; output_2000.(current_system).Concrete; output_3000.(current_system).Concrete];
                        output.(current_system).Steel=[output_1000.(current_system).Steel; output_2000.(current_system).Steel; output_3000.(current_system).Steel];
                        output.(current_system).Excavation=[output_1000.(current_system).Excavation; output_2000.(current_system).Excavation; output_3000.(current_system).Excavation];
                        output.(current_system).Construction_Misc=[output_1000.(current_system).Construction_Misc; output_2000.(current_system).Construction_Misc; output_3000.(current_system).Construction_Misc];
                        output.(current_system).Permeate_Pumping=[output_1000.(current_system).Permeate_Pumping; output_2000.(current_system).Permeate_Pumping; output_3000.(current_system).Permeate_Pumping];
                        output.(current_system).Recirculation_Pumping=[output_1000.(current_system).Recirculation_Pumping; output_2000.(current_system).Recirculation_Pumping; output_3000.(current_system).Recirculation_Pumping];
                        output.(current_system).Chemical_Pumping=[output_1000.(current_system).Chemical_Pumping; output_2000.(current_system).Chemical_Pumping; output_3000.(current_system).Chemical_Pumping];
                        output.(current_system).Gravity_Belt_Thickeners=[output_1000.(current_system).Gravity_Belt_Thickeners; output_2000.(current_system).Gravity_Belt_Thickeners; output_3000.(current_system).Gravity_Belt_Thickeners];
                        output.(current_system).Operation_Misc=[output_1000.(current_system).Operation_Misc; output_2000.(current_system).Operation_Misc; output_3000.(current_system).Operation_Misc];
                        %CSTR_none_Submerged
                        if ia==1 && ia1==1 && ib==1
                            output.(current_system).Gas_Sparging=[output_1000.(current_system).Gas_Sparging; output_2000.(current_system).Gas_Sparging; output_3000.(current_system).Gas_Sparging];
                            output.(current_system).Lift_Pumping=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
                            output.(current_system).Retentate_Pumping=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
                        elseif ia==1 && ia1==2 && ib==1 
                            %CSTR_GAC_Submerged
                            output.(current_system).Gas_Sparging=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
                            output.(current_system).Lift_Pumping=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
                            output.(current_system).Retentate_Pumping=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
                            %CSTR_none_Cross-flow
                        elseif ia==1 && ia1==1 && ib==2
                            output.(current_system).Gas_Sparging=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
                            output.(current_system).Lift_Pumping=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
                            output.(current_system).Retentate_Pumping=[output_1000.(current_system).Retentate_Pumping; output_2000.(current_system).Retentate_Pumping; output_3000.(current_system).Retentate_Pumping];
                            %AF_none_Crossflow OR AF_AeF_Crossflow
                        elseif ia==2 && ib==2
                            output.(current_system).Gas_Sparging=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
                            output.(current_system).Lift_Pumping=[output_1000.(current_system).Lift_Pumping; output_2000.(current_system).Lift_Pumping; output_3000.(current_system).Lift_Pumping];
                            output.(current_system).Retentate_Pumping=[output_1000.(current_system).Retentate_Pumping; output_2000.(current_system).Retentate_Pumping; output_3000.(current_system).Retentate_Pumping];
                            %AF_none_Submerged OR AF_AeF_Submerged
                        elseif ia==2 && ia1~=2 && ib==1
                            output.(current_system).Gas_Sparging=[output_1000.(current_system).Gas_Sparging; output_2000.(current_system).Gas_Sparging; output_3000.(current_system).Gas_Sparging];
                            output.(current_system).Lift_Pumping=[output_1000.(current_system).Lift_Pumping; output_2000.(current_system).Lift_Pumping; output_3000.(current_system).Lift_Pumping];
                            output.(current_system).Retentate_Pumping=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
                            %AF_GAC_Submerged
                        elseif ia==2 && ia1==2 && ib==1
                            output.(current_system).Gas_Sparging=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
                            output.(current_system).Lift_Pumping=[output_1000.(current_system).Lift_Pumping; output_2000.(current_system).Lift_Pumping; output_3000.(current_system).Lift_Pumping];
                            output.(current_system).Retentate_Pumping=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
                        end
                        
                        if ih==1
                            output.(current_system).Degassing_Memb_Op=[output_1000.(current_system).Degassing_Memb_Op; output_2000.(current_system).Degassing_Memb_Op; output_3000.(current_system).Degassing_Memb_Op];
                        else
                            output.(current_system).Degassing_Memb_Op=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',all_VarNames);
                        end
                        
                        % Direct emission
                        output.(current_system).Direct_Emissions=[output_1000.(current_system).Direct_Emissions; output_2000.(current_system).Direct_Emissions; output_3000.(current_system).Direct_Emissions];
                        
                        % Avoided Energy
                        output.(current_system).Avoided_Energy=[output_1000.(current_system).Avoided_Energy; output_2000.(current_system).Avoided_Energy; output_3000.(current_system).Avoided_Energy];
                    end
                end
            end
        end
    end
end

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
