A = output.LCC_LCA{:,[64 121:129]};
CAS = [113074732.181326,1.47635867609083e-08,0.250567573769560,0.0137041405942390,0.00265179614659425,0.0274762386909206,9.14297560436428e-08,9.99965792595472e-08,0.000360749388307178,3.04570336451803]; %SI from Smith et al., 2014
A_B = [113858978.137449,1.11686100121969e-08,0.144001626919317,0.00896265672173792,0.000267505429023883,0.0273975150379856,7.84281292272563e-08,8.60527671580387e-08,0.000252116442027889,2.65893321450320];
A = vertcat(A,CAS,A_B);
labels = {'NPV','Ozone Depletion','Global Warming','Smog','Acidification','Eutrophication','HH Cancer','HH Non-Cancer','HH Criteria','Ecotoxicity'};
m = min(A);
M = max(A);
B = (A-repmat(m,size(A,1),1))./(repmat(M,size(A,1),1)-repmat(m,size(A,1),1));
C = B(1:150,:);

b1 = B(86,:);
b2 = B(49,:);
b3 = B(151,:);
b4 = B(132,:);
b5 = B(152,:);
B1 = B([29:38 85:94 141:150],:);
B2 = B([39:56 95:112],:);

figure
px=parallelcoords(C,'Labels',labels,'Color',0.7*ones(1,3));
ax=gca;
ax.XTickLabelRotation = 45;
box on
hold all
p1=parallelcoords(b1,'Color',[0.329 0.604 0.773]); % MT Cross-flow
p2=parallelcoords(b2,'Color',[0.722 0.294 0.306]); % Submerged, GAC
p4=parallelcoords(b4,'Color',0.7*ones(1,3)); % All other solutions
p6=parallelcoords(B2,'Color',[0.722 0.294 0.306]); % Submerged, GAC
p5=parallelcoords(B1,'Color',[0.329 0.604 0.773]); % MT Cross-flow
p3=parallelcoords(b3,'Color','k','LineStyle','--','LineWidth',2'); % CAS
p7=parallelcoords(b5,'Color','k','LineStyle',':','LineWidth',2'); %A/B

hold off
legend([p1,p2,p3,p7,p4],'Cross-flow Multi-tube','Submerged with GAC','Conventional Activated Sludge','A/B','All Other Solutions','Location','northeastoutside')
ax=gca;
ax.XTickLabelRotation = 45;
ax.YLabel.String = 'Normalized Value [0-1]';
ax.YLabel.Position = [0.3 0.5];
ax.YGrid = 'on';
set(gcf,'units','points','position',[0,0,800,250])

%MT, X-flow: 29:38 85:94 141:150
%Sub, GAC: 39:56 95:112