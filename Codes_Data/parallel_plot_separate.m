A = output.LCC_LCA{:,[64 121:129]};
m = min(A);
M = max(A);
B = (A-repmat(m,size(A,1),1))./(repmat(M,size(A,1),1)-repmat(m,size(A,1),1));
% B1 = B([87:96 128 130 132 134 136],:);
% B2 = B([1 3 5 7 9:18],:);

B1 = B([75:84 131:140],:); %AF Cross-flow FS
B2 = B([85:94 141:150],:); %AF Cross-flow MT
B3 = B([1:18],:); %Submerged no GAC
B4 = B([39:56],:); %Submerged with GAC
B5 = B([1:8 39:46 58:64 95:102 113:120],:); %HF
B6 = B([9:28 47:56 65:84 103:112 121:140],:); %FS
B7 = B([29:38 85:94 141:150],:); %MT
B8 = B([1 2 5:10 13:20 23:30 33:40 43:48 51:58 61:66 69:76 79:86 89:96 99:104 107:114 117:122 125:132 135:142 145:150],:); %Other materials
B9 = B([3 4 11 12 21 22 31 32 41 42 49 50 59 60 67 68 77 78 87 88 97 98 105 106 115 116 123 124 133 134 143 144],:); %PTFE

%%
figure %AF Cross-flow
hold on
P2 = parallelcoords(B1,'Color',[0.639 0.831 0.961]); %Lighter Blue
P1=parallelcoords(B2,'Color',[0.329 0.604 0.773]); %Blue
hold off
box on
set(gca,'XtickLabel',[],'YtickLabel',[],'XLabel',[],'YLabel',[])
set(gcf,'units','points','position',[0,0,650,200])
axis([1,10,0,1])

%%
figure %CSTR Submerged
hold on
P4 = parallelcoords(B3,'Color',[0.890 0.525 0.533]); %Lighter Red
P3=parallelcoords(B4,'Color',[0.722 0.294 0.306]); %Red
hold off
box on
set(gca,'XtickLabel',[],'YtickLabel',[],'XLabel',[],'YLabel',[])
set(gcf,'units','points','position',[0,0,650,200])
axis([1,10,0,1])

%%
figure %Comparing membrane types
hold on
P7 = parallelcoords(B5,'Color',[0.557 0.765 0.612]); %Lighter Green
P6 = parallelcoords(B6,'Color',[0.322 0.549 0.392]); %Green
P5=parallelcoords(B7,'Color','k'); %Black
hold off
box on
set(gca,'XtickLabel',[],'YtickLabel',[],'XLabel',[],'YLabel',[])
set(gcf,'units','points','position',[0,0,650,200])
axis([1,10,0,1])

%%
figure %Comparing membrane materials
hold on
P9 = parallelcoords(B9,'Color',[0.698 0.612 0.918]); %Lighter Purple
P8=parallelcoords(B8,'Color',[0.463 0.357 0.659]); %Purple
hold off
box on
set(gca,'XtickLabel',[],'YtickLabel',[],'XLabel',[],'YLabel',[])
set(gcf,'units','points','position',[0,0,650,200])
axis([1,10,0,1])

% Green: [0.322 0.549 0.392]
% Purple: [0.463 0.357 0.659]
% Orange: [0.863 0.561 0.376]