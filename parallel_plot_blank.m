A = [0,0,0,0,0,0,0,0,0,0]; %SI from Smith et al., 2014
CAS = [113074732.181326,1.47635867609083e-08,0.250567573769560,0.0137041405942390,0.00265179614659425,0.0274762386909206,9.14297560436428e-08,9.99965792595472e-08,0.000360749388307178,3.04570336451803]; %SI from Smith et al., 2014
A_B = [113858978.137449,1.11686100121969e-08,0.144001626919317,0.00896265672173792,0.000267505429023883,0.0273975150379856,7.84281292272563e-08,8.60527671580387e-08,0.000252116442027889,2.65893321450320];
A = vertcat(A,CAS,A_B);
labels = {'NPV','Ozone Depletion','Global Warming','Smog','Acidification','Eutrophication','HH Cancer','HH Non-Cancer','HH Criteria','Ecotoxicity'};
m = min(A);
M = max(A);
B = (A-repmat(m,size(A,1),1))./(repmat(M,size(A,1),1)-repmat(m,size(A,1),1));

b1 = B(3,:);

figure
px=parallelcoords(b1,'Labels',labels,'Color','k');
ax=gca;
ax.XTickLabelRotation = 45;
box on
ax=gca;
ax.XTickLabelRotation = 45;
ax.YLabel.String = 'Normalized Value [0-1]';
ax.YLabel.Position = [0.3 0.5];
ax.YGrid = 'on';
set(gcf,'units','points','position',[0,0,800,250])