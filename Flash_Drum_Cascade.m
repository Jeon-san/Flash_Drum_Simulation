%% CH2140: Chemical Engineering Unit Operations

%% Assignment 1

% Name:Jeon San
% Mat no.:U1820704G

clear all
clc 

%% a) Specify known values & parameters 
p1 = 300; % pressure of drum 1 (kPa)
p1_bar = 300/100; % pressure of drum 1 (bar)

T_1c = 180; % Temperature of valve inlet 1 (celcius)
T_1k = 180+273.15; % Temperature of valve inlet 1 (K)

F = 500; % Feed molar flow rate (kmol/h)
z = [0.127 0.461 0.244 0.168]; % mol frac of components
Cp_liq = [131.1 168.0 224.7 313.2]; % Liquid heat capacity of pure components
Cp_vap = [111.8 139.8 188.3 266.1]; % Vapour heat capacity of pure components
latent = [22700 26500 35400 51300]; % Latent heat of vapourization

% Save antoine equation paramenters in matrix for all components
% Order: row 1 to 4-butane,pentane,heptane,decane
ant_para = [4.35576 1175.581 -2.071;3.9892 1070.617 -40.454;...
            4.02832 1268.636 -56.199;4.07857 1501.268 -78.670]; 
        
%% Q1) Determining minimum pressure(boiling pressure)

for i = 1:4
            para = ant_para(i,:);
            p_vap_180(i) = 10^(para(1)-para(2)/(T_1k+para(3))); % p_sat in bar
end

z = [0.127 0.461 0.244 0.168]; % mol frac of components
pres_boil = z*p_vap_180';
min_drop = pres_boil-p1_bar; % Min pressure drop
fprintf(' Question 1');
fprintf('\n Boiling pressure =%.3f bar(s), therefore minimum pressure drop is %.3f bar(s) ',pres_boil,min_drop);


%%%%%%%%%%%        TANK 1        %%%%%%%%%%%%
%% Rachford-Rice to get f

%% Input for iteration***

T_guess_c = 100; % Initial guess of T2 (Celcius) 
T_guess_k = T_guess_c + 273.15; % Initial guess of T2 (Kelvin) 

per_error=1;
while per_error >0.001
    
    for i = 1:4
        para = ant_para(i,:);
        p_sat(i) = 10^(para(1)-para(2)/((T_guess_c+273.15)+para(3))); % p_sat in bar using antoine eqn
    end

    K = p_sat/p1_bar; % Calculate Ki values (ideal mixture assumed)

    z = [0.127 0.461 0.244 0.168]; % mol frac of components

    syms f 

        for i= 1:4
            term(i) = ((K(i)-1)*z(i))/(1+(K(i)-1)*f);
        end

    rach_rice = sum(term)==0;

    f_solved = double(solve(rach_rice,f));
    f_01 = f_solved(f_solved>=0 & f_solved<=1); % Select solution between 0 and 1

    %% With f-value, create operating lines and xK=y lines

    % Note: butane=1,pentane=2,heptane=3,decane=4
    % exp - x1,y1 are liquid and vapour mole fraction of butane

    syms x1 x2 x3 x4 y1 y2 y3 y4

    % operating lines
    eqn1 = -((1-f_01)/f_01)*x1 + (1/f_01)*z(1) -y1 == 0;
    eqn2 = -((1-f_01)/f_01)*x2 + (1/f_01)*z(2) -y2 == 0;
    eqn3 = -((1-f_01)/f_01)*x3 + (1/f_01)*z(3) -y3 == 0;
    eqn4 = -((1-f_01)/f_01)*x4 + (1/f_01)*z(4) -y4 == 0;

    % xK=y lines
    eqn5 = x1*K(1)- y1 == 0;
    eqn6 = x2*K(2)- y2 == 0;
    eqn7 = x3*K(3)- y3 == 0;
    eqn8 = x4*K(4)- y4 == 0;

    [A,B] = equationsToMatrix([eqn1, eqn2, eqn3,eqn4, eqn5, eqn6,eqn7, eqn8]...
                              , [x1,x2,x3,x4,y1,y2,y3,y4]);
    XY_solutions = double(linsolve(A,B));
    X = XY_solutions(1:4,:);
    Y = XY_solutions(5:8,:);

    %% Calculate H_vap and H_liq
    V = F*f_01; % Vapour molar flowrate
    L = F-V; % Liquid molar flowrate
    H_liq = Cp_liq*X*T_guess_c; % Enthalpy of liquid (basis = celcius)
    H_vap_latent = latent*(Y);
    H_vap = Cp_liq*Y*T_guess_c +H_vap_latent ; % Enthalpy of vapour (basis = celcius)
   %H_vap = Cp_vap*Y*T_guess_c +H_vap_latent ; % Enthalpy of vapour (basis = celcius)
    RHS_enthalpy = V*H_vap+L*H_liq; % Total enthalpy assuming T_guess 
    LHS_enthalpy = (F*(z*Cp_liq'))*T_1c; % Inlet enthalpy
    per_error = abs((RHS_enthalpy-LHS_enthalpy)/RHS_enthalpy);
    T_guess_c = ((LHS_enthalpy/RHS_enthalpy)*T_guess_c-T_guess_c)*0.1+T_guess_c;
    
end

% export temperature, liquid flowrate and molar fraction to solve balance for tank 2
T2_in = T_guess_c;
F2 = L;
X_1 = XY_solutions(1:4,:);
Y_1 = XY_solutions(5:8,:);

fprintf('\n');
fprintf('\n Question 2');
fprintf('\n The temperature of drum 1 is %.2f celsius', T_guess_c);

%% Balance checks

% Overall mass balance
overall_diff = abs((F-V-L)/F);
% Component mass balance
component_diff = sum(abs(((F*z')-(V*Y_1)-(L*X_1))./(F*z'))); 
% p*/p=K vs yi/xi=K difference
k_diff = sum(abs((K-(Y_1./X_1)')./K)); 

% Rachford-rice
RR_up = (K-1).*z; % Numerator of Rachford-rice
RR_down = 1+(K-1).*(V/F); % Denominator of Rachford-rice
RR = sum(RR_up./RR_down); % Sum of all Rachford-rice terms

% Energy balance
heat_diff = per_error; % Use heat balance error from above directly

% Print balance status
if overall_diff >0.01
    fprintf(2,'\n Overall mass balance does not balance!')
else
    fprintf('\n Overall mass balance OK!')
end

if component_diff >0.01
    fprintf(2,'\n Component mass balance does not balance!')
else
    fprintf('\n Component mass balance OK!')
end

if k_diff >0.01
    fprintf(2,'\n xiKi=yi does not balance!')
else
    fprintf('\n xiKi=yi mass balance OK!')
end

if RR>0.01
    fprintf(2,'\n RR does not balance!')
elseif RR<-0.01
    fprintf(2,'\n RR does not balance!')
else
    fprintf('\n RR balance OK!')
end

if heat_diff >0.01
    fprintf(2,'\n Heat does not balance!')
else
    fprintf('\n Heat balance OK!')
end



%%%%%%%%%%%        TANK 2        %%%%%%%%%%%%

%% Input for iteration***

clearvars -except T2_in F2 X_1 Y_1 p1 p1_bar Cp_liq Cp_vap latent ant_para

T2_guess_c = 100; % Initial guess of T2 (Celcius) 
T2_guess_k = T2_guess_c + 273.15; % Initial guess of T2 (Kelvin) 

p2_kpa = 150; %Pressure of drum 2 (kPa)
p2_bar = 150/100; %Pressure of drum 2 (bar)

per_error=1;
while per_error >0.001
   

        for i = 1:4
            para = ant_para(i,:);
            p2_sat(i) = 10^(para(1)-para(2)/(T2_guess_c+273.15+para(3))); % p_sat in bar

        end

    K2 = p2_sat/p2_bar; % Calculate Ki values (ideal mixture assumed)

    z = X_1; % mol frac of components at bottom stream of drum 1

    syms f2 

        for i= 1:4
            term(i) = ((K2(i)-1)*z(i))./(1+(K2(i)-1)*f2);
        end

    rach_rice = sum(term)==0;

    f_solved = double(solve(rach_rice,f2));
    f_01 = f_solved(f_solved>=0 & f_solved<=1); % Select solution between 0 and 1

    %% With f-value, create operating lines and xK=y lines

    % Note: butane=1,pentane=2,heptane=3,decane=4
    % exp - x1,y1 are liquid and vapour mole fraction of butane

    syms x1 x2 x3 x4 y1 y2 y3 y4

    % operating lines
    eqn1 = -((1-f_01)/f_01)*x1 + (1/f_01)*z(1) -y1 == 0;
    eqn2 = -((1-f_01)/f_01)*x2 + (1/f_01)*z(2) -y2 == 0;
    eqn3 = -((1-f_01)/f_01)*x3 + (1/f_01)*z(3) -y3 == 0;
    eqn4 = -((1-f_01)/f_01)*x4 + (1/f_01)*z(4) -y4 == 0;

    % xK=y lines
    eqn5 = x1*K2(1)- y1 == 0;
    eqn6 = x2*K2(2)- y2 == 0;
    eqn7 = x3*K2(3)- y3 == 0;
    eqn8 = x4*K2(4)- y4 == 0;

    [A,B] = equationsToMatrix([eqn1, eqn2, eqn3,eqn4, eqn5, eqn6,eqn7, eqn8]...
                              , [x1,x2,x3,x4,y1,y2,y3,y4]);
    XY2_solutions = double(linsolve(A,B));
    X2 = XY2_solutions(1:4,:);
    Y2 = XY2_solutions(5:8,:);

    %% Calculate H_vap and H_liq
    V2 = F2*f_01; % Vapour molar flowrate
    L2 = F2-V2; % Liquid molar flowrate
    H2_liq = Cp_liq*X2*T2_guess_c; % Enthalpy of liquid (basis = celcius)
    H2_vap_latent = latent*(Y2);
    H2_vap = Cp_liq*Y2*T2_guess_c +H2_vap_latent ; % Enthalpy of vapour (basis = celcius)
   %H2_vap = Cp_liq*X2*T2_guess_c +H2_vap_latent ; % Enthalpy of vapour (basis = celcius)
    RHS2_enthalpy = V2*H2_vap+L2*H2_liq; % Total enthalpy assuming T_guess 
    LHS2_enthalpy = (Cp_liq*(F2*X_1))*T2_in; % Inlet enthalpy
    per_error = abs((RHS2_enthalpy-LHS2_enthalpy)/RHS2_enthalpy);
    T2_guess_c = (((LHS2_enthalpy/RHS2_enthalpy)*T2_guess_c)-T2_guess_c)*0.1+T2_guess_c;
    
end



fprintf('\n');
fprintf('\n The temperature of drum 2 is %.2f celsius', T2_guess_c);
fprintf('\n X_butane is %.3f', X2(1));
fprintf('\n X_pentane is %.3f', X2(2));
fprintf('\n X_heptane is %.3f', X2(3));
fprintf('\n X_decane is %.3f', X2(4));
fprintf('\n');
% Mass balance check

% Overall mass balance
overall_diff = abs((F2-V2-L2)/F2);
% Component mass balance
component_diff = sum(abs(((F2*z)-(V2*Y2)-(L2*X2))./(F2*z)));
% p*/p=K vs yi/xi=K difference
k_diff = sum(abs((K2-(Y2./X2)')/K2)); 

RR_up = (K2-1).*z'; % Numerator of Rachford-rice
RR_down = 1+(K2-1).*(V2/F2); % Denominator of Rachford-rice
RR = sum(RR_up./RR_down); % Sum of all Rachford-rice terms

% Energy balance
heat_diff = per_error; % Use heat balance error from above directly

% Print balance status
if overall_diff >0.01
    fprintf(2,'\n Overall mass balance does not balance!')
else
    fprintf('\n Overall mass balance OK!')
end

if component_diff >0.01
    fprintf(2,'\n Component mass balance does not balance!')
else
    fprintf('\n Component mass balance OK!')
end

if k_diff >0.01
    fprintf(2,'\n xiKi=yi does not balance!')
else
    fprintf('\n xiKi=yi mass balance OK!')
end

if RR>0.01
    fprintf(2,'\n RR does not balance!')
elseif RR<-0.01
    fprintf(2,'\n RR does not balance!')
else
    fprintf('\n RR balance OK!')
end

if heat_diff >0.01
    fprintf(2,'\n Heat does not balance!')
else
    fprintf('\n Heat mass balance OK!')
end
