function [OD,t,ID] = pipe(Q, VEL) % OD:[in], t:[in] VEL: [ft/s] Q:[ft^3/s]
Area = Q / VEL; % [ft^2] X-sectional area 
dia_ft = (4 * Area / pi)^0.5; % [ft] pipe required inner diameter
dia = dia_ft * 12; % [in] unit conversion from [ft] to [in]

if dia < 1/8
    OD = 0.405;
    t = 0.049;
elseif dia > 1/8 && dia < 1/4
    OD = 0.540;
    t = 0.065;
elseif dia > 1/4 && dia < 3/8
    OD = 0.675;
    t = 0.065;
elseif dia > 3/8 && dia < 1/2
    OD = 0.840;
    t = 0.083;
elseif dia > 1/2 && dia < 3/4
    OD = 1.050;
    t = 0.083;
elseif dia > 3/4 && dia < 1
    OD = 1.315;
    t = 0.109;
elseif dia > 1 && dia < 1.25
    OD = 1.660;
    t = 0.109;
elseif dia > 1.25 && dia < 1.5
    OD = 1.900;
    t = 0.109;
elseif dia > 1.5 && dia < 2
    OD = 2.375;
    t = 0.109;
elseif dia > 2 && dia < 2.5
    OD = 2.875;
    t = 0.120;    
elseif dia > 2.5 && dia < 3
    OD = 3.500;
    t = 0.120;
elseif dia > 3 && dia < 3.5
    OD = 4.000;
    t = 0.120; 
elseif dia > 3.5 && dia < 4
    OD = 4.500;
    t = 0.120; 
elseif dia > 4 && dia < 4.5
    OD = 5.000;
    t = 0.120; 
elseif dia > 4.5 && dia < 5
    OD = 5.563;
    t = 0.134; 
elseif dia > 5 && dia < 6
    OD = 6.625;
    t = 0.134; 
elseif dia > 6 && dia < 7
    OD = 7.625;
    t = 0.134; 
elseif dia > 7 && dia < 8
    OD = 8.625;
    t = 0.148;    
elseif dia > 8 && dia < 9
    OD = 9.625;
    t = 0.148;
elseif dia > 9 && dia < 10
    OD = 10.750;
    t = 0.165; 
elseif dia > 10 && dia < 11
    OD = 11.750;
    t = 0.165; 
elseif dia > 11 && dia < 12
    OD = 12.750;
    t = 0.180; 
elseif dia > 12 && dia < 14
    OD = 14.000;
    t = 0.188;        
elseif dia > 14 && dia < 16
    OD = 16.000;
    t = 0.199; 
elseif dia > 16 && dia < 18
    OD = 18.000;
    t = 0.188;     
elseif dia > 18 && dia < 20
    OD = 20.000;
    t = 0.218;   
elseif dia > 20 && dia < 24
    OD = 24.000;
    t = 0.250; 
elseif dia > 24 && dia < 26
    OD = 26.000;
    t = 0.250;
elseif dia > 26 && dia < 28
    OD = 28.000;
    t = 0.250; 
elseif dia > 28 && dia < 30
    OD = 30.000;
    t = 0.312; 
elseif dia > 30 && dia < 32
    OD = 32.000;
    t = 0.312;
elseif dia > 32 && dia < 34
    OD = 34.000;
    t = 0.312;
elseif dia > 34 && dia < 36
    OD = 36.000;
    t = 0.312; 
elseif dia > 36 && dia < 42
    OD = 42.000;
    t = 0.312; 
elseif dia > 42 && dia < 48
    OD = 48.000;
    t = 0.312;
else
    OD = 48.000;
    t = 0.312;
end

ID = OD - 2 * t; % [in] inner diameter 

end