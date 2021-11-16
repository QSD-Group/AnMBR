function [summary] = RankingDataExtraction(output)

x=fieldnames(output);
summary=zeros(3000,112);
addition=zeros(3000,16);
count=1;
for i=1:150
    summary(:,count)=output.(x{i}).LCC{:,20};
    addition(:,1)=output.(x{i}).Concrete{:,2};
    addition(:,2)=output.(x{i}).Steel{:,2};
    addition(:,3)=output.(x{i}).Excavation{:,2};
    addition(:,4)=output.(x{i}).Construction_Misc{:,2};
    addition(:,5)=output.(x{i}).Permeate_Pumping{:,2};
    addition(:,6)=output.(x{i}).Recirculation_Pumping{:,2};
    addition(:,7)=output.(x{i}).Chemical_Pumping{:,2};
    addition(:,8)=output.(x{i}).Gravity_Belt_Thickeners{:,2};
    addition(:,9)=output.(x{i}).Operation_Misc{:,2};
    addition(:,10)=output.(x{i}).Gas_Sparging{:,2};
    addition(:,11)=output.(x{i}).Lift_Pumping{:,2};
    addition(:,12)=output.(x{i}).Retentate_Pumping{:,2};
    addition(:,13)=output.(x{i}).Degassing_Memb_Op{:,2};
    addition(:,14)=output.(x{i}).Direct_Emissions{:,2};
    addition(:,15)=output.(x{i}).Avoided_Energy{:,2};
    addition(:,16)=addition(:,1)+addition(:,2)+addition(:,3)+addition(:,4)+addition(:,5)+addition(:,6)+addition(:,7)+addition(:,8)+addition(:,9)+addition(:,10)+addition(:,11)+addition(:,12)+addition(:,13)+addition(:,14)+addition(:,15);
    summary(:,count+150)=addition(:,16);
    count=count+1;
end
end