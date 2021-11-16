%% Variables
% AFF		Air flow fraction
% AFR		Air flow required
% AL		Filter aluminum dosage [lb/day]
% ALUM 		Alum dosage as Al2(SO4)3·14 H2O [lb/day]
% BBA  		Blower building area [ft2]
% BBC  		Blower building costs [$]
% CALFED 	Capital cost of alum storage and feed system [$]
% CCLIME  	Capital cost of lime feed system [$]
% CEC 		Construction and equipment cost [$]
% CEPCIP	Current CE plant cost index for pipe, valves, etc.
% CF 		Correction factor for other minor construction costs
% CFEFED  	Capital cost of iron salt feed system [$]
% CFM		Minimum air flow (blower) capacity [scfm/1,000 ft3 tank volume]
% CFMB	 	Capacity of individual blower [scfm]
% CFMD		Design capacity of blowers [scfm]
% COSTAP	Cost of air piping [$]
% COSTBE 	Purchase cost of blower of CFMB capacity [$]
% COSTCS	Cost of reinforced concrete slab in place [$]
% COSTCW	Cost of reinforced concrete wall in place [$]
% COSTE  	Cost of earthwork [$]
% COSTHR	Cost of handrails in place [$]
% COSTPB  	Cost of pump building [$]
% COSTRC	Total cost of reinforced concrete
% COSTRO  	Purchase cost of blower of CFMB capacity as percent of cost of standard blower [%]
% COSTSB  	Purchase cost of standard blower [$]
% CPMFED  	Capital cost of polymer feed system [$]
% CTFED  	Total capital cost of chemical feed system [$]
% FBMAF     Fine bubble system minimum air flow
% FE		Iron salt dosage rate expressed as equivalent Fe molecules [lb/day]
% FPC		Firm pumping capacity [MGD]
% GPM		Design capacity of pumps [gpm]
% IBC  		Installed blower costs [$]
% IPC  		Installed pumping equipment cost [$]
% IRON		Ferric chloride dosage rate [lb/day]
% LCAT		EPA cost index for larger city advanced treatment
% LCV  		Liquid chemical solution feed per day [gpd]
% LHR		Length of handrails [ft]
% LIME		Lime dosage rate as Ca(OH)2 [lb/day]
% LRR  		Labor required ratio
% MCFMD     Modified CFMD
% MLC  		Maintenance labor cost [$]
% MLCA  	Maintenance labor cost for alum
% MLCF  	Maintenance labor cost for iron salt
% MLCI  	Maintenance labor cost for intermediate pumps
% MLCL 		Maintenance labor cost for lime
% MLCP 		Maintenance labor cost for polymer
% MLCR 		Maintenance labor cost for sludge return pumps
% MLH  		Maintenance labor required [person-hr/year]
% MLH_PUP	Maintenance labor required for pumping
% MSC  		Material and supply cost [$]
% N 		Number of blowers needed to provide maximum air requirements
% N_module  Number of modules
% NB  		Total number of blowers required
% NP_IR     Total number of pumps needed for internal recirculation pumping
% NP_PERM   Total number of pumps needed for permeate pumping
% NP_R      Number of pumps for cross-flow pumping
% OLC  		Operational labor cost [$]
% OLCA  	Operation labor cost for alum 
% OLCF  	Operation labor cost for iron salt
% OLCI		Operation labor cost for intermediate pumps
% OLCL  	Operation labor cost for lime
% OLCP  	Operation labor cost for polymers
% OLCR  	Operation labor cost for sludge return pumps
% OLH  		Operating labor required [person-hr/year]
% OLH_PUP	Operating labor required for pumping
% OMCC_PUP  Operation and maintenance material and supply cost [$]
% OMMP      Operation and maintenance material and supply costs, as percent of the total bare construction cost [%]
% OMMP_CMF  Operation and maintenance material supply cost as fraction of capital cost
% OMMP_PUP  Operation and maintenance material and supply costs, as percent of the total bare construction cost [%]
% PBA		Pump building area [ft2]
% Q_gas_cfm Design capacity for blower system [scfm]
% Q_IR      Flow rate for internal recirculation pumping [MGD]
% Q_PERM    Flow rate for permeate pumping [MGD]
% Q_IR      Flow rate for internal recirculation pumping [MGD]
% Q_R       Firm pumping capacity for retentate
% PLMER     Polymer dosage rate [lb/day]
% RSR		Return sludge ratio to average wastewater flow (see table on p. 1197)
% SA		Surface area [m2]
% STE		Standard oxygen transfer efficiency [%]
% TBCC_BLW  Total bare construction cost of blowers [$]
% TBCC_PUP  Total bare construction cost of pumping [$]
% TCFM  	Design capacity for blower system [scfm]
% TMLC  	Total maintenance labor cost [$/year]
% TOLC  	Total operation labor cost [$/year]
% UPIBC 	Unit price input for building costs [$/ft2]
% UPICS		Unit price input of reinforced concrete slab in place [$/yd3]
% UPICW		Unit price input of reinforced concrete wall in place [$/yd3]
% UPIEX		Unit price input for earthwork [$/yd3]
% UPIHR		Unit price input for handrails in place [$/ft]
% V_ANA     Volume of anaerobic filter [m3]
% V_CSTR    Volume of CSTR [m3]
% VEX       Volume earthwork [ft3]
% V_m_AER   Volume of filter media for aerobic polishing filter [ft3]
% V_m_AAN   Volume of filter media for anaerobic filter [ft3]
% VSC       Volume slab concrete [ft3]
% VWC       Volume wall concrete [ft3]
%% Crossflow or Submerged?
LCC_Output=cell(40,2);
decide = input('Crossflow or Submerged? ');
if strncmp(decide,'Submerged',9)==1
    VWC=21476;
    VSC=33872.74722;
    VEX=379004;
    Q_PERM_mgd = 20;
    NP_PERM = 7;
    Q_IR_mgd = 20;
    NP_IR = 1;
    Q_gas_cfm = 37263.82333;
    NB=4;
    Q_NaOCl_mgd = 0.000120879;
    NP_NaOCl = 1;
    Q_CA_mgd = 0.000032967032967033;
    NP_CA = 1;
    N_module = 8008;
elseif strncmp(decide,'Crossflow',9)==1
    VWC=89682.94443;
    VSC=45273.151;
    VEX=35835;
    Q_PERM_mgd = 20;
    NP_PERM = 386.5820854;
    Q_IR_mgd = 20;
    NP_IR = 21;
    Q_NaOCl_mgd = 0.000120879;
    NP_NaOCl = 1;
    Q_CA_mgd = 0.000032967032967033;
    NP_CA = 1;
    N_module = 11598;
else
    disp('Enter Crossflow or Submerged (with single quotes).')
    return
end

LCC_Output{1,1}='Flow rate [MGD]';
LCC_Output{2,1}='Permeate flow rate [MGD]';
LCC_Output{2,2}=Q_PERM_mgd;
LCC_Output{3,1}='Internal recirculation flow rate [MGD]';
LCC_Output{3,2}=Q_IR_mgd;
if strncmp(decide,'Submerged',9)==1
    LCC_Output{4,1}='Gas flow rate [cfm]';
    LCC_Output{4,2}=Q_gas_cfm;
end
LCC_Output{5,1}='Volume of tank [m3]';
LCC_Output{6,1}='Volume of wall concrete [m3]';
LCC_Output{6,2}=VWC;
LCC_Output{7,1}='Volume of slab concrete [m3]';
LCC_Output{7,2}=VSC;
LCC_Output{8,1}='Volume of earthwork [m3]';
LCC_Output{8,2}=VEX;
LCC_Output{9,1}='Number of permeate pumps';
LCC_Output{9,2}=NP_PERM;
LCC_Output{10,1}='Number of internal recirculation pumps';
LCC_Output{10,2}=NP_IR;
if strncmp(decide,'Submerged',9)==1
    LCC_Output{11,1}='Number of blowers';
    LCC_Output{11,2}=NB;
end
LCC_Output{12,1}='NaOCl requirement';
LCC_Output{12,2}=Q_NaOCl_mgd;
LCC_Output{13,1}='NaOCl pumps';
LCC_Output{13,2}=NP_NaOCl;
LCC_Output{14,1}='Citric acid requirment';
LCC_Output{14,2}=Q_CA_mgd;
LCC_Output{15,1}='Citric acid pumps';
LCC_Output{15,2}=NP_CA;
LCC_Output{16,1}='Number of modules';
LCC_Output{16,2}=N_module;

%% AP
% disp('Instruction: Enter the requested values below, and type in COSTAP to get the cost of Air Piping')
% CFM = input('Enter the value of the Minimum Air Flow (blower) Capacity, scfm/1000 ft^3 tank volume: ');
% V = input('Enter the value of the Volume of the tank: ');
% CFMD = CFM * V/1000; % Design Capacity of Blowers, scfm
% if (CFMD<=1000)
%     COSTAP = 617.2 * 3.33*CFMD^0.2553; % 3.33: Air Flow Fraction, calculated as STE/6, and STE stands for Standard Oxygen Transfer Efficiency _%
% elseif (CFMD>1000)&&(CFMD<=10000)       % The default STE value is 20. If users have their own STE value, then plug the value into STE/6.
%     COSTAP = 1.43 * 3.33* CFMD^1.1337;  % If STE/6 > 1, take this value and replace 3.33, if the value is smaller than 1, then use 1 to replace 3.33.
% elseif (CFMD>10000)
%     COSTAP = 28.59 * 3.33*CFMD^0.8085; % In calculating COSTAP, the CEPCIP (Current CE Plant Cost Index for Pipe, Valves, etc) value is set to be 241
% end                                    % If users have updated CEPCIP value, then multiply the COSTAP result with CEPCIP/241.
% LCC_Output{17,1}='Cost of air piping [$]';

%% BLW
disp('Instruction: Enter the requested value below, and enter the abbreviations of the cost information you want to know ')
disp('Here is a list of abbreviations you will need')
disp('BBC: Blower Building Cost, $')
disp('COSTBE: Purchase Cost of Blower of CFMB Capacity, $')
disp('TBCC: Total Bare Construction Cost, $')
disp('IBC: Installed Blower Cost, $')
disp('If you want to enter a new value for the Design Capacity for Blower System, type in BLW to re-enter the program')
disp('For blower system, the Material and Supply Cost comes mainly from the Purchase cost of Blowers.')
disp('And the Operation Cost comes mainly from Electricity Cost. Please type in CEHE for Electricity Cost information')
TCFM = input('Enter the Design Capacity for Blower System in scfm: ');

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
if (TCFM<=30000)
    COSTRO = 0.7*CFMB^0.6169;
elseif (TCFM>30000)&&(TCFM<=72000)
    COSTRO = 0.377* CFMB^0.5928;
elseif (TCFM>72000)
    COSTRO = 0.964* CFMB^0.4286;
end

if (TCFM<=30000)
    COSTBE  = 50000*(COSTRO) / 100; % Purchase Cost of Blower of CFMB Capacity, $, 50000:Purchase cost of 3,000 scfm at 8 psig
elseif (TCFM>30000)&&(TCFM<=72000)
    COSTBE  = 190000*(COSTRO) / 100; % Purchase Cost of 12,000 scfm at 8 psig
elseif (TCFM>72000)
    COSTBE  = 4140000*(COSTRO) / 100; % Purchase Cost of 50,000 scfm at 8 psig
end

IBC = 2 * COSTBE; % IBC: Installed Blower Costs, $
BBA = 128* TCFM^0.256; % Blower Building Area, ft^2
BBC = BBA * 90; % 90: Unit Price Input for Building Costs, $/ft^2
TBCC_BLW = (IBC * NB + BBC) * 1.11; % 1.11: Correction Factor for Other Minor Construction Costs includes piping, concrete, steel, electrical, paint and installation labor

LCC_Output{18,1}='Design capacity for blower system [scfm]';
LCC_Output{18,2}=TCFM;
LCC_Output{19,1}='Number of blowers';
LCC_Output{19,2}=NB;
LCC_Output{20,1}='Installed blower costs [$]';
LCC_Output{20,2}=IBC;
LCC_Output{21,1}='Total bare construction costs for blowers [$]';
LCC_Output{21,2}=TBCC_BLW;

%% CEHE
disp('Instruction: Enter the values as instructed below')
disp('Here is a table of all the abbreviations you will need')
disp('COSTCW - Cost of Wall Reinforced Concrete, $')
disp('COSTCS - Cost of Slab Reinforced Concrete, $')
disp('COSTRC - Total Cost of Reinforced Concrete, $')
disp('COSTE - Cost of Earthwork, $')
disp('COSTHR - Cost of Handrail, $')
disp('COSTEC - Cost of Electricity, $')
disp('Enter 0 for items that are irrelevant to the cost information you need')
%VWC = input('Enter a value for Volume of wall reinforced concrete in ft^3: ');
COSTCW = VWC/27 * 18.52;  % 18.52: Unit price for wall concrete in dollar, information obtained from the 2007 vendor data in CapdetWorks  
                          % 27: Conversion Factor from ft^3 to yd^3
%VSC = input('Enter a value for Volume of slab reinforced concrete in ft^3: ');
COSTCS = VSC/27*12.96;    % 12.96: Unit price for slab concrete in dollar, information obtained from the 2007 vendor data in CapdetWorks
                          % 27: Conversion Factor from ft^3 to yd^3
COSTRC = COSTCW + COSTCS;

%VEX = input('Enter a value for the volume of excavation of earthwork in ft^3: ');
COSTE = VEX/27*0.3; % Assuptions: The basin will be constructed using equal cut and fill
                    % The side slopes will be 3 to 1
                    % The basin will have 2 ft. of freeboard
                    % 27: Conversion Factor from ft^3 to yd^3
%LHR = input('Enter a value for Length of Handrails in place in ft: ');
%COSTHR = LHR * 75; % 75: Unit price for handrail in place in dollar, information obtained from the 2007 vendor data in CapdetWorks
%UPEC = input('Enter a value for Electricity Consumption in kWh: ');
%COSTEC = 0.1 * UPEC; % 0.1: Unit price for Electricity per kWh in dollar, information obtained from the 2007 vendor data in CapdetWorks
LCC_Output{22,1}='Cost of concrete [$]';
LCC_Output{22,2}=COSTRC;
LCC_Output{23,1}='Cost of earthwork [$]';
LCC_Output{23,2}=COSTE;

%% CMF
% disp ('Instruction: Enter the values below as instructed, if you do not add the chemical(s) being asked to your treatment system, enter the value 0.');
% disp ('After you finish entering the values, you can type in the abbreviations of the cost information you want and get the results you need')
% disp ('If you want to enter a new combination of values of each chemicals added to the system, just type in CMF to re-enter the program')
% disp ('Here is a table of all the abbreviations:')
% disp ('OLCA - Operation Labor Cost for Alum, $/yr')
% disp ('MLCA - Maintenance Labor Cost for Alum, $/yr')
% disp ('OLCF - Operation Labor Cost for Iron Salt , $/yr')
% disp ('MLCF - Maintenance Labor Cost for Iron Salt, $/yr')
% disp ('OLCL - Operation Labor Cost for Lime, $/yr')
% disp ('MLCL - Maintenance Labor Cost for Lime, $/yr')
% disp ('OLCP - Operation Labor Cost for Polymers, $/yr')
% disp ('MLCP - Maintenance Labor Cost for Polymers , $/yr')
% disp ('CALFED -Capital Cost of Alum Storage and Feed System, $ ')
% disp ('CFEFED -Capital Cost of Iron Salt Feed System, $') 
% disp ('CCLIME -Capital Cost of Lime Feed System, $')
% disp ('CPMFED -Capital Cost of Polymer Feed System, $')
% disp ('TOLC - Total Operation Labor Cost, $/yr')
% disp ('TMLC - Total Maintenance Labor Cost, $/yr')
% disp ('CTFED - Total Capital Cost of Chemical Feed System, $')
% disp ('MSC - Material and Supply Cost. $/yr') 
% % Output = zeros(500,20);
% % for i = 1:500
% %     ALUM = i*200;
% %     IRON = i*250;
% %     LIME = i*25;
% %     PLMER = i;
% 
% ALUM = input('Enter a value for Alum Dosage in lb/day as Al2(SO4)3·14H2O: ');
% LCVa = (ALUM / 0.4902) * (27/603);    % LCVa: Liquid Chemical Solution Feed Per Day, gpd
% if (LCVa == 0)                        % It is assumed that liquid alum contains 0.4902 lb of aluminum per gallon
%     OLCA = 0;                         % 27/603 = Conversion from alum to filter aluminum
% elseif (LCVa <= 90)
%     OLCA = 600*0.54*25;               % 0.54: LRR (Labor Required Ratio) = OLH_PUP / ( OLH_PUP + MLH_PUP) = 0.54
% elseif (LCVa > 90) && (LCVa <= 350)   % The same value 0.54 is being used for dividing OLC and MLC cost for all chemicals
%     OLCA = (189.2 * LCVa ^0.2565)*0.54*25;
% elseif (LCVa > 350) && (LCVa <= 1050)
%     OLCA =(33.4 * LCVa^0.5527)*0.54*25;
% elseif (LCVa > 1050) && (LCVa <= 10000)
%     OLCA = (51.8* LCVa^0.4894)*0.54*25;
% elseif (LCVa>10000);
%     OLCA =(12.2 * LCVa^0.647)*0.54*25;
% end
% 
% if (LCVa == 0)
%     MLCA = 0;
% elseif (LCVa <= 90)
%     MLCA = 600 * (1-0.54)*25;
% elseif (LCVa > 90) && (LCVa <= 350)
%     MLCA = (189.2 * LCVa ^0.2565)*(1-0.54)*25;
% elseif (LCVa > 350) && (LCVa <= 1050)
%     MLCA =(33.4 * LCVa^0.5527) *(1-0.54)*25;
% elseif (LCVa > 1050) && (LCVa <= 10000)
%     MLCA =(51.8* LCVa^0.4894)*(1-0.54)*25;
% elseif (LCVa > 10000)
%     MLCA =(12.2 * LCVa^0.647)*(1-0.54)*25;
% end
% AL = ALUM * (27/603);
% if (LCVa == 0)
%     CALFED = 0;
% elseif (AL<=425)
%     CALFED = 59000;
% elseif (AL>425) && (AL<=945)
%     CALFED = 10260*AL^0.289; % Here we assumed that the LCAT (EPA Cost Index for Larger City Advanced Treatment) value 
% elseif (AL>945) && (AL<=4027) % to be 132, which is the value at the first quarter of 1977.   
%     CALFED = 661.7*AL^0.689;   % If more recent LCAT value is obtained, time the results by LCATnew/132
% elseif (AL>4027)               % The same rule apply to all the Capital Cost calculations.
%     CALFED = 2.419*AL^1.365;
% end
% 
% IRON = input('Enter a value for Ferric Chloride Dosage Rate, lb/day: ');
% FE = IRON *(55.8/162.2); % 55.8/162.2= Conversion from iron to FeCl3
% LCVf = FE/4.11;   %LCVf: Liquid Chemical Solution Feed Per Day, gpd
% if (LCVf == 0)
%     OLCF = 0;
% elseif (LCVf <= 90)
%     OLCF = 600*0.54*25;
% elseif (LCVf > 90) && (LCVf <= 350)
%     OLCF = (189.2 * LCVf ^0.2565)*0.54*25;
% elseif (LCVf > 350) && (LCVf <= 1050)
%     OLCF =(33.4 * LCVf^0.5527)*0.54*25;
% elseif (LCVf > 1050) && (LCVf <= 10000)
%     OLCF = (51.8* LCVf^0.4894)*0.54*25;
% elseif (LCVf>10000);
%     OLCF =(12.2 * LCVf^0.647)*0.54*25;
% end
% 
% if (LCVf == 0)
%     MLCF = 0;
% elseif (LCVf <= 90)
%     MLCF = 600 * (1-0.54)*25;
% elseif (LCVf > 90) && (LCVf <= 350)
%     MLCF = (189.2 * LCVf ^0.2565)*(1-0.54)*25;
% elseif (LCVf > 350) && (LCVf <= 1050)
%     MLCF =(33.4 * LCVf^0.5527) *(1-0.54)*25;
% elseif (LCVf > 1050) && (LCVf <= 10000)
%     MLCF =(51.8* LCVf^0.4894)*(1-0.54)*25;
% elseif (LCVf > 10000)
%     MLCF =(12.2 * LCVf^0.647)*(1-0.54)*25;
% end
% 
% if (FE == 0)
%     CFEFED = 0;
% elseif (FE<=1000)
%     CFEFED = 59000;
% elseif (FE>1000) && (FE<=2350)
%     CFEFED =3352 * FE^0.4152 ;
% elseif (FE>2350) && (FE<=16767)
%     CFEFED = 86.92 * FE^0.8857;
% elseif (FE>16767)
%     CFEFED =0.458* FE^1.425;
% end
% 
% LIME = input('Enter a value for Lime Feed Rate as lb/day of Ca(OH)2: ');
% LCVl = LIME / 0.5;   %LCVl: Liquid Chemical Solution Feed Per Day, gpd
% if (LCVl == 0)
%     OLCL = 0;
% elseif (LCVl <= 90)
%     OLCL = (600+92.5 *LCVl^0.2827)*0.54*25;
% elseif (LCVl > 90) && (LCVl <= 350)
%     OLCL = (189.2 * LCVl ^0.2565+92.5* LCVl^0.2827)*0.54*25;
% elseif (LCVl > 350) && (LCVl <= 1050)
%     OLCL =(33.4 * LCVl^0.5527+92.5 *LCVl^0.2827)*0.54*25;
% elseif (LCVl > 1050) && (LCVl <= 10000)
%     OLCL = (51.8* LCVl^0.4894+92.5* LCVl^0.2827)*0.54*25;
% elseif (LCVl >10000);
%     OLCL =(12.2 * LCVl^0.647+92.5 *LCVl^0.2827)*0.54*25;
% end
% 
% if (LCVl == 0)
%     MLCL = 0;
% elseif (LCVl <= 90)
%     MLCL = (600+92.5* LCVl^0.2827) * (1-0.54)*25;
% elseif (LCVl > 90) && (LCVl <= 350)
%     MLCL = (189.2 * LCVl ^0.2565+92.5 *LCVl^0.2827)*(1-0.54)*25;
% elseif (LCVl > 350) && (LCVl <= 1050)
%     MLCL =(33.4 * LCVl^0.5527+92.5 *LCVl^0.2827) *(1-0.54)*25;
% elseif (LCVl > 1050) && (LCVl <= 10000)
%     MLCL =(51.8* LCVl^0.4894+92.5* LCVl^0.2827)*(1-0.54)*25;
% elseif (LCVl > 10000)
%     MLCL =(12.2 * LCVl^0.647+92.5 *LCVl^0.2827)*(1-0.54)*25;
% end
% 
% if (LIME == 0)
%     CCLIME = 0;
% elseif (LIME<=750)
%     CCLIME = 26100;
% elseif (LIME>750)
%     CCLIME =326.7 * LIME^0.6614;
% end
% 
% 
% PLMER = input('Enter a value for Polymer Dosage Rate, lb/day: ');
% LCVp = (100 * PLMER) / (0.25 * 8.34);   %LCVp: Liquid Chemical Solution Feed Per Day, gpd
% if (LCVp == 0)
%     OLCP = 0;
% elseif (LCVp <= 90)
%     OLCP = (600+92.5 *LCVp^0.2827)*0.54*25;
% elseif (LCVp > 90) && (LCVp <= 350)
%     OLCP = (189.2 * LCVp ^0.2565+92.5* LCVp^0.2827)*0.54*25;
% elseif (LCVp > 350) && (LCVp <= 1050)
%     OLCP =(33.4 * LCVp^0.5527+92.5 *LCVp^0.2827)*0.54*25;
% elseif (LCVp > 1050) && (LCVp <= 10000)
%     OLCP = (51.8* LCVp^0.4894+92.5* LCVp^0.2827)*0.54*25;
% elseif (LCVp >10000);
%     OLCP =(12.2 * LCVp^0.647+92.5 *LCVp^0.2827)*0.54*25;
% end
% 
% if (LCVp == 0)
%     MLCP = 0;
% elseif (LCVp <= 90)
%     MLCP = (600+92.5* LCVp^0.2827) * (1-0.54)*25;
% elseif (LCVp > 90) && (LCVp <= 350)
%     MLCP = (189.2 * LCVp ^0.2565+92.5 *LCVp^0.2827)*(1-0.54)*25;
% elseif (LCVp > 350) && (LCVp <= 1050)
%     MLCP =(33.4 * LCVp^0.5527+92.5 *LCVp^0.2827) *(1-0.54)*25;
% elseif (LCVp > 1050) && (LCVp <= 10000)
%     MLCP =(51.8* LCVp^0.4894+92.5* LCVp^0.2827)*(1-0.54)*25;
% elseif (LCVp > 10000)
%     MLCP =(12.2 * LCVp^0.647+92.5 *LCVp^0.2827)*(1-0.54)*25;
% end
% 
% if (PLMER == 0)
%     CPMFED = 0;
% elseif (PLMER<=28)
%     CPMFED = 4524 * PLMER ^0.4075;
% elseif (PLMER>28)
%     CPMFED =1018 * PLMER ^0.8562;
% end
% 
% 
% CTFED = CALFED + CFEFED + CCLIME + CPMFED;
% MSC = 0.02 * CTFED; % 0.02: Material and Supply Cost as Fraction of Capital Cost, Fraction = 0.02
% TOLC = OLCA+OLCF+OLCL+OLCP;
% TMLC = MLCA+MLCF+MLCL+MLCP;
% % 
% % Output(i,1)=ALUM;
% % Output(i,2)=OLCA;
% % Output(i,3)=MLCA;
% % Output(i,4)=CALFED;
% % Output(i,5)=FE;
% % Output(i,6)=OLCF;
% % Output(i,7)=MLCF;
% % Output(i,8)=CFEFED;
% % Output(i,9)=LIME;
% % Output(i,10)=OLCL;
% % Output(i,11)=MLCL;
% % Output(i,12)=CCLIME;
% % Output(i,13)=PLMER;
% % Output(i,14)=OLCP;
% % Output(i,15)=MLCP;
% % Output(i,16)=CPMFED;
% % Output(i,17)=i;
% % Output(i,18)=TOLC;
% % Output(i,19)=TMLC;
% % Output(i,20)=CTFED;
% % 
% % 
% % end

% LCC_Output{24,1}='Alum dosage [lb/d]';
% LCC_Output{24,2}=ALUM;
% LCC_Output{25,1}='Alum feed solution [gpd]';
% LCC_Output{25,2}=LCVa;
% LCC_Output{26,1}='Operation labor cost for alum [$/yr]';
% LCC_Output{26,2}=OLCA;
% LCC_Output{27,1}='Maintenance labor cost for alum [$/yr]';
% LCC_Output{27,2}=MLCA;
% LCC_Output{28,1}='Capital cost of alum storage and feed system [$]';
% LCC_Output{28,2}=CALFED;
% LCC_Output{29,1}='Ferric chloride dosage [lb/d]';
% LCC_Output{29,2}=IRON;
% LCC_Output{30,1}='Ferric chloride feed solution [gpd]';
% LCC_Output{30,2}=LCVf;
% LCC_Output{31,1}='Operation labor cost for ferric chloride [$/yr]';
% LCC_Output{31,2}=OLCF;
% LCC_Output{32,1}='Maintenance labor cost for ferric chloride [$/yr]';
% LCC_Output{32,2}=MLCF;
% LCC_Output{33,1}='Capital cost of ferric chloride feed system [$]';
% LCC_Output{33,2}=CFEFED;
% LCC_Output{34,1}='Lime dosage [lb/d]';
% LCC_Output{34,2}=LIME;
% LCC_Output{35,1}='Lime feed solution [gpd]';
% LCC_Output{35,2}=LCVl;
% LCC_Output{36,1}='Operation labor cost for lime [$/yr]';
% LCC_Output{36,2}=OLCL;
% LCC_Output{37,1}='Maintenance labor cost for lime [$/yr]';
% LCC_Output{37,2}=MLCL;
% LCC_Output{38,1}='Capital cost of lime feed system [$]';
% LCC_Output{38,2}=CCLIME;
% LCC_Output{39,1}='Polymer dosage [lb/d]';
% LCC_Output{39,2}=PLMER;
% LCC_Output{40,1}='Polymer feed solution [gpd]';
% LCC_Output{40,2}=LCVp;
% LCC_Output{41,1}='Operation labor cost for polymer [$/yr]';
% LCC_Output{41,2}=OLCP;
% LCC_Output{42,1}='Maintenance labor cost for polymer [$/yr]';
% LCC_Output{42,2}=MLCP;
% LCC_Output{43,1}='Capital cost of polymer feed system [$]';
% LCC_Output{43,2}=CPMFED;
% LCC_Output{44,1}='Total capital cost of chemical feed system [$]';
% LCC_Output{44,2}=CTFED;
% LCC_Output{45,1}='Material supply cost';
% LCC_Output{45,2}=MSC;
% LCC_Output{46,1}='Total operational labor cost [$/yr]';
% LCC_Output{46,2}=TOLC;
% LCC_Output{47,1}='Total maintenance labor cost [$/yr]';
% LCC_Output{47,2}=TMLC;

%% PUP
% disp('Instruction: Enter the values requested. Enter 0 for items that are irrelevant to the cost information you need')
% disp('Here is a table of all the abbreivations you will need:')
% disp('OLCI: Operation Labor Cost for Intermediate Pumps, $/yr')
% disp('OLCR: Operation Labor Cost for Return Sludge Pumps, $/yr')
% disp('MLCI: Maintenance Labor Cost for Intermediate Pumps, $/yr')
% disp('MLCR: Maintenance Labor Cost for Return Sludge Pumps, $/yr')
% disp('TOLC: Total Operation Labor Cost, $/yr')
% disp('TMLC: Total Maintenance Labor Cost, $/yr')
% disp('COSTPB: Cost of Pump Building, $')
% disp('TBCC: Total Bare Construction Cost, $')
% disp('MSC: Material and Supply Cost, $')
% disp('After entering the values requested, enter the abbreviations of the cost information you want')
% disp('Once you want to enter new values for Average Daily Wastewater Flow or Return Sludge Ratio to Average Wastewater Flow')
% disp('Just enter PUP to re-enter the program')
% disp('RSR: Return Sludge Ratio to Average Wastewater Flow')
% disp('Below is a table of different RSR values for different Activated Sludge Processes:')
% disp('Activated Sludge Processes       RSR')
% disp('Conventional                     1.0')
% disp('Complete-Mix                     1.0')
% disp('Step-Aeration                    1.0')
% disp('Modified-Aeration                1.0')
% disp('Contact-Stabilization            1.0')
% disp('Kraus-Process                    1.0')
% disp('Pure-Oxygen System               1.0')
% disp('Extended-Aeration                1.5')
% disp('High-Rate Aeration               5.0')
% disp('')
%Qavg = input('Enter the Average Daily Wastewater Flow in mgd: ');
 Output=zeros(100,16);
RSR = 1;%input('Enter the Return Sludge Ratio to Average Wastewater Flow: ');
for i = 1:100
    Qavg=i;
GPMI = (2*Qavg*10^6) / 1440; % GPMI: Design Capacity of Intermediate Pumps in gpm. 2 = excess capacity factor to handle peak flows
GPMR = (Qavg*RSR*10^6) / 1440; % GPMR: Design Capacity of Return Sludge Pumps in gpm
FPCI = 1440*GPMI / (10^6); % FPC: Firm Pumping Capacity in mgd
FPCR = 1440*GPMR / (10^6);
if (FPCI == 0)
    OLCI = 0;
elseif (FPCI > 0) && (FPCI <= 7)
    OLCI = 440*25* FPCI^0.1285;
elseif (FPCI > 7) && (FPCI <= 41)
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
elseif (FPCI > 7) && (FPCI <= 41)
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
elseif (FPCR > 7) && (FPCR <= 41)
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
elseif (FPCR > 7) && (FPCR <= 41)
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
IPC = 2.065*10^5+7.721*10^4*Qavg;
%IPC = 2*(COSTROI/100) *48549.837*a*c + 2*(COSTROR/100) *48549.837*b*d; % IPC: Installed equipment cost. 48549.837:Standard Cost of a 3000 gpm pump and driver.
VEX = 1;%input('Enter a value for the volume of excavation of earthwork for Pump Building in ft^3: ');
COSTE = VEX/27*0.3;
TBCC_PUP = (COSTE + COSTPB + IPC)*1.18;
MSC = TBCC_PUP*(0.007/100);

 Output(i,1)=i;
% Output(i,2)=FPCI;
% Output(i,3)=OLCI;
% Output(i,4)=MLCI;
% Output(i,5)=FPCR;
% Output(i,6)=OLCR;
% Output(i,7)=MLCR;
% Output(i,8)=GPMI;
% Output(i,9)=a;
% Output(i,10)=GPMBI;
% Output(i,11)=COSTROI;
% Output(i,12)=COSTROR;
% Output(i,13)=IPC;
 Output(i,14)=TBCC_PUP;
% Output(i,15)=MSC;
% 
 end

% LCC_Output{48,1}='Return sludge ratio';
% LCC_Output{48,2}=RSR;
% LCC_Output{49,1}='Design capacity of intermediate pumps [gpm]';
% LCC_Output{49,2}=GPMI;
% LCC_Output{50,1}='Design capacity of return sludge pumps [gpm]';
% LCC_Output{50,2}=GPMR;
% LCC_Output{51,1}='Firm pumping capacity of intermediate pumps [gpm]';
% LCC_Output{51,2}=FPCI;
% LCC_Output{52,1}='Firm pumping capacity of return sludge pumps [gpm]';
% LCC_Output{52,2}=FPCR;
% LCC_Output{53,1}='Operation labor cost for intermediate pumps [$]';
% LCC_Output{53,2}=OLCI;
% LCC_Output{54,1}='Maintenance labor cost for intermediate pumps [$]';
% LCC_Output{54,2}=MLCI;
% LCC_Output{55,1}='Operation labor cost for return sludge pumps [$]';
% LCC_Output{55,2}=OLCR;
% LCC_Output{56,1}='Maintenance labor cost for return sludge pumps [$]';
% LCC_Output{56,2}=MLCR;
% LCC_Output{57,1}='Installed pump cost [$]';
% LCC_Output{57,2}=IPC;
% LCC_Output{58,1}='Total bare construction cost for pumps [$]';
% LCC_Output{58,2}=TBCC_PUP;
% LCC_Output{59,1}='Material supply costs [$]';
% LCC_Output{59,2}=MSC;

%% H&S Equations

%Flow and load conditions
Q_mgd=input('Enter influent flow rate [MGD]: ');

%Membrane tank design criteria
Module_SA = 370; %Membrane surface area per Zenon Membrane Module [ft2]
Mod_per_cas = 44; %Zenon 500D Modules per Cassette (Max = 48)
Contact_SA = Module_SA*Mod_per_cas; %Contact Surface Area per Zenon Membrane Cassette [ft2]
Cas_per_tank = 19; %Cassettes per membrane tank
Spare_cassettes = 3; %Spare cassettes per membrane tank
if Cas_per_tank<11
    L_membrane_unit = Cas_per_tank*8.5+Spare_cassettes*8.5; %[ft]
    W_membrane_unit = 10; %[ft]
else
    L_membrane_unit = Cas_per_tank*8.5/2+Spare_cassettes*8.5/2; %[ft]
    W_membrane_unit = 21; %[ft]
end
N_tanks = 7; %Total number of membrane tanks
N_cassettes = Cas_per_tank*N_tanks;
Flux_all = (Q_mgd*1000000)/(N_tanks*Cas_per_tank*Contact_SA);
Flux_oneoff = (Q_mgd*1000000)/((N_tanks-1)*Cas_per_tank*Contact_SA);
Design_criteria = input('Enter design criterion [gal/ft2/d]: ');
Tank_footprint = N_tanks*L_membrane_unit*W_membrane_unit;

%Membrane Air
Air_demand_10_10 = 6; %[cfm/module]
Air_demand_10_30 = 3; %[cfm/module]
Air_demand_10_10_instantaneous = N_tanks*Air_demand_10_10*(Mod_per_cas*Cas_per_tank); %[cfm/module]
Air_demand_10_30_instantaneous = N_tanks*Air_demand_10_30*(Mod_per_cas*Cas_per_tank); %[cfm/module]
NB = 5; %For this size plant, assume one blower per two trains, re-evaluate during detailed design.
BLW_size = Air_demand_10_10_instantaneous/(NB-1); %Blower size to handle 10:10
Tot_cfm = BLW_size*(NB-1);

%MLSS recycle pump criteria
Q_R_necess_MGD = 4*Q_mgd; %Necessary recycle pump flow rate [MGD]
Q_R_necess_gpm = (Q_R_necess_MGD*1000000)/(24*60);%Necessary recycle pump flow rate [gpm]

%Footprint for permeate pumps, membrane chem feed systems/cleaning, membrane blowers, 
%MBR controls, valves, elec room, HVAC, etc:
SF_per_MGD = 1200;
FT_req = SF_per_MGD*Q_mgd;

%Screenings disposal (additional screenings for MBR, i.e., in addition to 6 mm screens)
Screen_rate_max = 4; %[CF/hr/MGD]
Screen_rate = 0.5;
Screen_per_day = Screen_rate_max*Q_mgd*24*Screen_rate; %[CF/day]
Compaction = 0.75;
Daily_screen = Screen_per_day*(1-Compaction);
Daily_screen_CY = Daily_screen/27;
Cost_disposal = 225; %$/20 CY container
Cost_screening = Daily_screen_CY*365*Cost_disposal/20;

%Power costs
Fine_screens = 1*Q_mgd; %at 1 hp/MGD
Est_TDH_perm = 40; %[ft]
Est_eff_perm = 0.8;
Avg_bhp_perm_pumps = (Q_mgd*1000000/1440)*Est_TDH_perm/(Est_eff_perm*3960);
Est_TDH_blower = 6; %[psig]
Est_eff_blower = 0.7;
Avg_bhp_mem_blowers = (Air_demand_10_30_instantaneous*0.23*(((14.7+Est_TDH_blower)/14.7)^0.283-1))/Est_eff_blower;
HP_subtotal = Fine_screens+Avg_bhp_perm_pumps+Avg_bhp_mem_blowers;
Misc_power = HP_subtotal*0.015;
HP_total = HP_subtotal+Misc_power;
Total_avg_power = HP_total*0.746*24;
Power_cost = 0.07; %$/kWh
Daily_power_cost = Total_avg_power*Power_cost;
Cost_power = Daily_power_cost*365;

%Membrane cleaning
Usage_rate_NaOCl = 2200; %[gal/yr/MGD]
Chem_cost_NaOCl = 0.54; %[$/gal], 12.5% solution
Annual_cost_NaOCl = Chem_cost_NaOCl*Q_mgd*Usage_rate_NaOCl*(12.5/15); %[$/yr]
Usage_rate_citric = 600; %[gal/yr/MGD]
Chem_cost_citric = 0.85; %[$/gal], 100% solution, 13.8 lb/gal
Annual_cost_citric = Chem_cost_citric*Q_mgd*Usage_rate_citric*(13.8); %[$/yr]
Usage_rate_Bisulfite = 350; %[gal/yr/MGD]
Chem_cost_Bisulfite = 0.3; %[$/gal], 38% solution, 3.5 lb/gal
Annual_cost_Bisulfite = Chem_cost_Bisulfite*Q_mgd*Usage_rate_Bisulfite*3.5*0.38; %[$/yr]
Cost_cleaning = Annual_cost_NaOCl+Annual_cost_citric+Annual_cost_Bisulfite;

%Membrane replacement cost
Cost_per_module = 1100; %[$/module]
Labor_replacement = Cost_per_module*0.15;
Total_cost_rep = Cost_per_module+Labor_replacement;
Expected_membrane_life = 7; %[yrs]
N_modules = Mod_per_cas*N_cassettes;
Replacement_cost = N_modules*Total_cost_rep;
Cost_replacement_annual = Replacement_cost/Expected_membrane_life;

%O&M Costs
Total_annual_cost = Cost_screening+Cost_power+Cost_cleaning+Cost_replacement_annual;
Total_daily_cost = Total_annual_cost/365;

LCC_Output{60,1}='Influent flow rate [MGD]';
LCC_Output{60,2}=Q_mgd;
LCC_Output{61,1}='Membrane surface area per module [ft2]';
LCC_Output{61,2}=Module_SA;
LCC_Output{62,1}='Modules per cassette';
LCC_Output{62,2}=Mod_per_cas;
LCC_Output{63,1}='Contact surface area [ft2]';
LCC_Output{63,2}=Contact_SA;
LCC_Output{64,1}='Cassettes per membrane tank';
LCC_Output{64,2}=Cas_per_tank;
LCC_Output{65,1}='Spare cassettes';
LCC_Output{65,2}=Spare_cassettes;
LCC_Output{66,1}='Length of membrane unit [ft]';
LCC_Output{66,2}=L_membrane_unit;
LCC_Output{67,1}='Width of membrane unit [ft]';
LCC_Output{67,2}=W_membrane_unit;
LCC_Output{68,1}='Number of membrane tanks';
LCC_Output{68,2}=N_tanks;
LCC_Output{69,1}='Number of cassettes';
LCC_Output{69,2}=N_cassettes;
LCC_Output{70,1}='Flux when all tanks are online [gal/ft2/d]';
LCC_Output{70,2}=Flux_all;
LCC_Output{71,1}='Flux with one tank offline [gal/ft2/d]';
LCC_Output{71,2}=Flux_oneoff;
LCC_Output{72,1}='Design flux [gal/ft2/d]';
LCC_Output{72,2}=Design_criteria;
LCC_Output{73,1}='Tank footprint [ft2]';
LCC_Output{73,2}=Tank_footprint;
LCC_Output{74,1}='Air demand at 10:10 air cycling [cfm/module]';
LCC_Output{74,2}=Air_demand_10_10;
LCC_Output{75,1}='Instantaneous air demand at 10:10';
LCC_Output{75,2}=Air_demand_10_10_instantaneous;
LCC_Output{76,1}='Air demand at 10:30 air cycling [cfm/module]';
LCC_Output{76,2}=Air_demand_10_30;
LCC_Output{77,1}='Instantaneous air demand at 10:30';
LCC_Output{77,2}=Air_demand_10_30_instantaneous;
LCC_Output{78,1}='Number of blowers';
LCC_Output{78,2}=NB;
LCC_Output{79,1}='Blower size to handle 10:10 [cfm]';
LCC_Output{79,2}=BLW_size;
LCC_Output{80,1}='Total cfm';
LCC_Output{80,2}=Tot_cfm;
LCC_Output{81,1}='Necessary recycle pump flow rate [MGD]';
LCC_Output{81,2}=Q_R_necess_MGD;
LCC_Output{82,1}='Max screening rate [CF/hr/MGD]';
LCC_Output{82,2}=Screen_rate_max;
LCC_Output{83,1}='Screenings rate [CF/hr/MGD]';
LCC_Output{83,2}=Screen_rate;
LCC_Output{84,1}='Compaction';
LCC_Output{84,2}=Compaction;
LCC_Output{85,1}='Daily screenings quantity [CY/day]';
LCC_Output{85,2}=Daily_screen_CY;
LCC_Output{86,1}='Disposal cost [$/20 CY]';
LCC_Output{86,2}=Cost_disposal;
LCC_Output{87,1}='Total annual screening cost [$/yr]';
LCC_Output{87,2}=Cost_screening;
LCC_Output{88,1}='Estimated permeate pump TDH [ft]';
LCC_Output{88,2}=Est_TDH_perm;
LCC_Output{89,1}='Estimated permeate pump efficiency';
LCC_Output{89,2}=Est_eff_perm;
LCC_Output{90,1}='Total average permeate pump bhp';
LCC_Output{90,2}=Avg_bhp_perm_pumps;
LCC_Output{91,1}='Estimated blower TDH [psig]';
LCC_Output{91,2}=Est_TDH_blower;
LCC_Output{92,1}='Estimated blower efficiency';
LCC_Output{92,2}=Est_eff_blower;
LCC_Output{93,1}='Total average blower bhp';
LCC_Output{93,2}=Avg_bhp_mem_blowers;
LCC_Output{94,1}='Miscellaneous power requirements';
LCC_Output{94,2}=Misc_power;
LCC_Output{95,1}='Total horsepower';
LCC_Output{95,2}=HP_total;
LCC_Output{96,1}='Price of electricity [$/kWh]';
LCC_Output{96,2}=Power_cost;
LCC_Output{97,1}='Total annual power cost [$/yr]';
LCC_Output{97,2}=Cost_power;
LCC_Output{98,1}='NaOCl usage rate [gal/yr/MGD]';
LCC_Output{98,2}=Usage_rate_NaOCl;
LCC_Output{99,1}='Cost of NaOCl [$/gal]';
LCC_Output{99,2}=Chem_cost_NaOCl;
LCC_Output{100,1}='Annual cost of NaOCl [$/yr]';
LCC_Output{100,2}=Annual_cost_NaOCl;
LCC_Output{101,1}='Citric acid usage rate [gal/yr/MGD]';
LCC_Output{101,2}=Usage_rate_citric;
LCC_Output{102,1}='Cost of citric acid [$/gal]';
LCC_Output{102,2}=Chem_cost_citric;
LCC_Output{103,1}='Annual cost of citric acid [$/yr]';
LCC_Output{103,2}=Annual_cost_citric;
LCC_Output{104,1}='Sodium bisulfite usage rate [gal/yr/MGD]';
LCC_Output{104,2}=Usage_rate_Bisulfite;
LCC_Output{105,1}='Cost of sodium bisulfite [$/gal]';
LCC_Output{105,2}=Chem_cost_Bisulfite;
LCC_Output{106,1}='Annual cost of sodium bisulfite [$/yr]';
LCC_Output{106,2}=Annual_cost_Bisulfite;
LCC_Output{107,1}='Total membrane cleaning cost [$/yr]';
LCC_Output{107,2}=Cost_cleaning;
LCC_Output{108,1}='Cost per module [$/module]';
LCC_Output{108,2}=Cost_per_module;
LCC_Output{109,1}='Labor for replacement [$/module]';
LCC_Output{109,2}=Labor_replacement;
LCC_Output{110,1}='Total cost/module [$/module]';
LCC_Output{110,2}=Total_cost_rep;
LCC_Output{111,1}='Expected module life [yrs]';
LCC_Output{111,2}=Expected_membrane_life;
LCC_Output{112,1}='Number of modules';
LCC_Output{112,2}=N_modules;
LCC_Output{113,1}='Total replacement cost [$]';
LCC_Output{113,2}=Replacement_cost;
LCC_Output{114,1}='Annual replacement cost [$/yr]';
LCC_Output{114,2}=Cost_replacement_annual;
LCC_Output{115,1}='Total annual cost [$/yr]';
LCC_Output{115,2}=Total_annual_cost;
LCC_Output{116,1}='Total daily cost [$/day]';
LCC_Output{116,2}=Total_daily_cost;
