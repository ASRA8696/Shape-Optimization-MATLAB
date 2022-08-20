clear
close
clc

L = 1;
Nx = 201;
x = linspace(0,L,Nx);
tfin = 20;

patm = 101325; %[Pa]

pres = 50*patm; %Dimensionless
Tres = 1273; %Dimensionless

pext0 = patm; %Dimensionless
Text0 = 300; %Dimensionless


gamma = 1.4;         % cp/cv in ideal gas for air         
mw_air = 28.97/1000; % molecular weight of air [kg/mol]
R = 8.314;           % Ideal gas constant [J/(mol*K)]

%Derived params
Rair = R/mw_air;     % Mass ideal gas constant for air [J/(kg*K)]
cv = Rair/(gamma-1); % specific heat of air at constant volume [J/(kg*K)]

% Derived reference state variables
eref = cv*Tres;                % specific energy at reference state [J/kg]
rhoref = pres/(Rair*Tres);     % mass density at reference state [kg/m^3]         
aref = sqrt(gamma*Rair*Tres);  % Speed of sound at reference state [m/s]


A1fun = @(x) 1 + 15*(x - 0.5).^2; % 0 < x < 1 [1]
A1 = A1fun(x);

A2fun = @(x) 1 + 15*(x - 0.3).^2; % 0 < x < 1 [1]
A2 = A2fun(x);

Rvals = [1.5, 1.5, 1.4, 1.2, 1, 1, 1.2, 1.5, 1.8, 2.05, 2.2, 2.3, 2.35, 2.4, 2.45, 2.47];
A3 = spline(linspace(0,L,length(Rvals)), Rvals.^2, linspace(0,L,Nx));


%% 
[tvec_realtime1, rho1, p1, T1, v1, Ma1] = ...
    NozzleConservativeFun(tfin, x, A1, pres, pext0, Tres, Text0, ...
                                    gamma, mw_air, R);
                                
[tvec_realtime2, rho2, p2, T2, v2, Ma2] = ...
    NozzleConservativeFun(tfin, x, A2, pres, pext0, Tres, Text0, ...
                                    gamma, mw_air, R);  
                                                             
[tvec_realtime3, rho3, p3, T3, v3, Ma3] = ...
    NozzleConservativeFun(tfin, x, A3, pres, pext0, Tres, Text0, ...
                                    gamma, mw_air, R); 


len1 = length(tvec_realtime1);
m1(1:len1, :) = rho1(1:len1, :).*A1.*v1(1:len1, :);
len2 = length(tvec_realtime2);
m2(1:len2, :) = rho2(1:len2, :).*A2.*v2(1:len2, :);
len3 = length(tvec_realtime3);
m3(1:len3, :) = rho3(1:len3, :).*A3.*v3(1:len3, :);



%% Fixed parameters
g = 9.81; %[m/s^2]
Astar = pi*0.15^2; %[m^2]  30cm diameter
mref = rhoref*Astar*aref;
patm = 101325; %[Pa]


%% Thrust and impulse in function of pres, with pext = 1atm and Tres = 1273 K 
pext = patm;
Tres = 1273; %[K]

preslist = linspace(20*patm, 200*patm, 100);
Thrust1 = nan(size(preslist));
Thrust2 = nan(size(preslist));
Thrust3 = nan(size(preslist));
Isp1 = nan(size(preslist));
Isp2 = nan(size(preslist));
Isp3 = nan(size(preslist));

for i = 1:length(preslist)
    
    pref = preslist(i);
    Tref = Tres;
    
    eref = cv*Tref;                % specific energy at reference state [J/kg]
    rhoref = pref/(Rair*Tref);     % mass density at reference state [kg/m^3]         
    aref = sqrt(gamma*Rair*Tref);  % Speed of sound at reference state [m/s]
    mref = rhoref*Astar*aref;
    
    Thrust1(i) = m1(end, end)*v1(end,end)*mref*aref ...
                + (p1(end, end)*pref - pext)*A1(end)*Astar;

    Thrust2(i) = m2(end, end)*v2(end,end)*mref*aref ...
                + (p2(end, end)*pref - pext)*A2(end)*Astar;

    Thrust3(i) = m3(end, end)*v3(end,end)*mref*aref ...
                + (p3(end, end)*pref - pext)*A3(end)*Astar;


    Isp1(i) = Thrust1(i)/(m1(end,end)*mref*g);
    Isp2(i) = Thrust2(i)/(m2(end,end)*mref*g);
    Isp3(i) = Thrust3(i)/(m3(end,end)*mref*g);
end

figure(1);
subplot(3,1,1)
plot(preslist/patm, Thrust1/1e3, 'r');
hold on
plot(preslist/patm, Thrust2/1e3, 'g');
plot(preslist/patm, Thrust3/1e3, 'b');
hold off
title('Thrust depending on reservoir pressure')
xlabel('Reservoir pressure [atm]');
ylabel('Thrust [kN]');
legend({'Nozzle1','Nozzle2','Nozzle3'});
grid on
figure(2);
subplot(3,1,1)
plot(preslist/patm, Isp1, 'r');
hold on
plot(preslist/patm, Isp2, 'g');
plot(preslist/patm, Isp3, 'b');
hold off
title('Specific impulse depending on reservoir pressure')
xlabel('Reservoir pressure [atm]');
ylabel('Specific impulse [s]');
legend({'Nozzle1','Nozzle2','Nozzle3'});
grid on


%% Thrust and impulse in function of pext, with pres = 100atm and Tres = 1273 K 
pres = 100*patm;
Tres = 1273; %[K]

pextlist = linspace(1e-2*patm, 2*patm, 100);
Thrust1 = nan(size(preslist));
Thrust2 = nan(size(preslist));
Thrust3 = nan(size(preslist));
Isp1 = nan(size(preslist));
Isp2 = nan(size(preslist));
Isp3 = nan(size(preslist));

for i = 1:length(preslist)
    
    pref = pres;
    Tref = Tres;
    pext = pextlist(i);
    
    eref = cv*Tref;                % specific energy at reference state [J/kg]
    rhoref = pref/(Rair*Tref);     % mass density at reference state [kg/m^3]         
    aref = sqrt(gamma*Rair*Tref);  % Speed of sound at reference state [m/s]
    mref = rhoref*Astar*aref;
    
    Thrust1(i) = m1(end, end)*v1(end,end)*mref*aref ...
                + (p1(end, end)*pref - pext)*A1(end)*Astar;

    Thrust2(i) = m2(end, end)*v2(end,end)*mref*aref ...
                + (p2(end, end)*pref - pext)*A2(end)*Astar;

    Thrust3(i) = m3(end, end)*v3(end,end)*mref*aref ...
                + (p3(end, end)*pref - pext)*A3(end)*Astar;


    Isp1(i) = Thrust1(i)/(m1(end,end)*mref*g);
    Isp2(i) = Thrust2(i)/(m2(end,end)*mref*g);
    Isp3(i) = Thrust3(i)/(m3(end,end)*mref*g);
end

figure(1)
subplot(3,1,2)
plot(pextlist/patm, Thrust1/1e3, 'r');
hold on
plot(pextlist/patm, Thrust2/1e3, 'g');
plot(pextlist/patm, Thrust3/1e3, 'b');
hold off
title('Thrust depending on external pressure')
xlabel('External pressure [atm]');
ylabel('Thrust [kN]');
legend({'Nozzle1','Nozzle2','Nozzle3'});
grid on
figure(2)
subplot(3,1,2)
plot(pextlist/patm, Isp1, 'r');
hold on
plot(pextlist/patm, Isp2, 'g');
plot(pextlist/patm, Isp3, 'b');
hold off
title('Specific impulse depending on external pressure')
xlabel('External pressure [atm]');
ylabel('Specific impulse [s]');
legend({'Nozzle1','Nozzle2','Nozzle3'});
grid on




%% Thrust and impulse in function of Tres, with pres = 100atm and pext = 1atm
pres = 100*patm;
pext = 1*patm;

Treslist = linspace(100, 3000, 200);
Thrust1 = nan(size(Treslist));
Thrust2 = nan(size(Treslist));
Thrust3 = nan(size(Treslist));
Isp1 = nan(size(Treslist));
Isp2 = nan(size(Treslist));
Isp3 = nan(size(Treslist));

for i = 1:length(Treslist)
    
    pref = pres;
    Tref = Treslist(i);
    
    eref = cv*Tref;                % specific energy at reference state [J/kg]
    rhoref = pref/(Rair*Tref);     % mass density at reference state [kg/m^3]         
    aref = sqrt(gamma*Rair*Tref);  % Speed of sound at reference state [m/s]
    mref = rhoref*Astar*aref;
    
    Thrust1(i) = m1(end, end)*v1(end,end)*mref*aref ...
                + (p1(end, end)*pref - pext)*A1(end)*Astar;

    Thrust2(i) = m2(end, end)*v2(end,end)*mref*aref ...
                + (p2(end, end)*pref - pext)*A2(end)*Astar;

    Thrust3(i) = m3(end, end)*v3(end,end)*mref*aref ...
                + (p3(end, end)*pref - pext)*A3(end)*Astar;


    Isp1(i) = Thrust1(i)/(m1(end,end)*mref*g);
    Isp2(i) = Thrust2(i)/(m2(end,end)*mref*g);
    Isp3(i) = Thrust3(i)/(m3(end,end)*mref*g);
end

figure(1)
subplot(3,1,3)
plot(Treslist, Thrust1/1e3, 'r');
hold on
plot(Treslist, Thrust2/1e3, 'g');
plot(Treslist, Thrust3/1e3, 'b');
hold off
title('Thrust depending on reservoir temperature')
xlabel('Reservoir temperature [K]');
ylabel('Thrust [kN]');
legend({'Nozzle1','Nozzle2','Nozzle3'});
grid on
figure(2)
subplot(3,1,3)
plot(Treslist, Isp1, 'r');
hold on
plot(Treslist, Isp2, 'g');
plot(Treslist, Isp3, 'b');
hold off
title('Specific impulse depending on reservoir temperature')
xlabel('Reservoir temperature [K]');
ylabel('Specific impulse [s]');
legend({'Nozzle1','Nozzle2','Nozzle3'});
grid on


