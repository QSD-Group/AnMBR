current_system=cell(150,1);
i=1;
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
                            step_A='A'; %CSTR
                        else
                            step_A='B'; %AF
                        end
                        
                        if ia1==1 %Additional options
                            step_A1='A'; %none
                        elseif ia1==2
                            step_A1='B'; %GAC
                        else
                            step_A1='C'; %AER
                        end
                        
                        if ib==1 %Membrane configuration
                            step_B='A'; %Submerged
                        elseif ib==2
                            step_B='B'; %Crossflow
                        end
                        
                        if ic==1 %Membrane type
                            step_C='A'; %Hollow fiber
                        elseif ic==2
                            step_C='B'; %Flat sheet
                        else
                            step_C='C'; %Multi-tube
                        end
                        
                        if id==1 %Membrane material
                            step_D='A'; %'PET';
                        elseif id==2
                            step_D='B'; %'PTFE';
                        elseif id==3
                            step_D='C'; %'ceramic';
                        elseif id==4
                            step_D='D'; %'sinter_steel';
                        elseif id==5
                            step_D='E'; %'PES';
                        else
                            step_D='F'; %'PVDF';
                        end
                        
                        if ih==1 %use DM or not
                            step_H='A'; %DM
                        else
                            step_H='B'; %No DM
                        end
                        current_system(i,1)={horzcat(step_A,step_A1,step_B,step_C,step_D,step_H)};
                        i=i+1;
                    end
                end
            end
        end
    end
end