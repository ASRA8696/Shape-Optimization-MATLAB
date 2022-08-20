%%%%%%%%%%%%%        SUPERSONIC NOZZLE FLOW            %%%%%%%%%%%%%%



clear
close
clc


L = 1;                  %Nozzle length [m]
Nx = 201;
x = linspace(0,L,Nx);
tfin = 2.5;             % Dimensionless time. 
                        % Represents number of times sound waves 
                        % can travel length L of Nozzle

patm = 101325;         %[Pa]
pres = 100*patm;
pext0 = patm;
Tres = 1273;           %[K]
Text0 = 300;           %[K]

Tref = Tres;
pref = pres;

gamma = 1.4;                   % cp/cv in ideal gas for air         
mw_air = 28.97/1000;           % molecular weight of air [kg/mol]
R = 8.314;                     % Ideal gas constant [J/(mol*K)]

%Derived params
Rair = R/mw_air;               % Mass ideal gas constant for air [J/(kg*K)]
cv = Rair/(gamma-1);           % specific heat of air at constant volume [J/(kg*K)]

% Derived reference state variables
eref = cv*Tref;                % specific energy at reference state [J/kg]
rhoref = pref/(Rair*Tref);     % mass density at reference state [kg/m^3]         
aref = sqrt(gamma*Rair*Tref);  % Speed of sound at reference state [m/s]

%Cross section at throat
Astar = pi*0.15.^2; %30cm diameter

%Dimensionless cross-section profiles
A1fun = @(x) 1 + 15*(x - 0.5).^2; % 0 < x < 1 [1]
A1 = A1fun(x);

A2fun = @(x) 1 + 15*(x - 0.3).^2; % 0 < x < 1 [1]
A2 = A2fun(x);

Rvals = [1.5, 1.5, 1.4, 1.2, 1, 1, 1.2, 1.5, 1.8, 2.05, 2.2, 2.3, 2.35, 2.4, 2.45, 2.47];
A3 = spline(linspace(0,L,length(Rvals)), Rvals.^2, linspace(0,L,Nx));

%% Plot profiles

rad_vec1 = sqrt(A1/pi);
rad_vec1 = rad_vec1./max(rad_vec1); %Normalized radius
lsup1 = +rad_vec1;
linf1 = -rad_vec1;

rad_vec2 = sqrt(A2/pi);
rad_vec2 = rad_vec2./max(rad_vec2); %Normalized radius
lsup2 = +rad_vec2;
linf2 = -rad_vec2;

rad_vec3 = sqrt(A3/pi);
rad_vec3 = rad_vec3./max(rad_vec3); %Normalized radius
lsup3 = +rad_vec3;
linf3 = -rad_vec3;

figure(1);
subplot(3,1,1)
plot(x/L,lsup1, 'r', 'LineWidth', 2);
hold on
plot(x/L,linf1, 'r', 'LineWidth', 2);
ylim([min(linf1), max(lsup1)]);
title('Nozzle 1')
hold off

subplot(3,1,2)
plot(x/L,lsup2, 'g', 'LineWidth', 2);
hold on
plot(x/L,linf2, 'g', 'LineWidth', 2);
ylim([min(linf2), max(lsup2)]);
title('Nozzle 2')
hold off

subplot(3,1,3)
plot(x/L,lsup3, 'b', 'LineWidth', 2);
hold on
plot(x/L,linf3, 'b', 'LineWidth', 2);
ylim([min(linf3), max(lsup3)]);
title('Nozzle 3')
hold off



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

%Re-scale back to dimensional
rho1 = rho1*rhoref;
p1 = p1*pref/patm;  %[atm]
T1 = T1*Tref;
v1 = v1*aref;
m1 = m1*rhoref*Astar*aref;

rho2 = rho2*rhoref;
p2 = p2*pref/patm;  %[atm]
T2 = T2*Tref;
v2 = v2*aref;
m2 = m2*rhoref*Astar*aref;

rho3 = rho3*rhoref;
p3 = p3*pref/patm;  %[atm]
T3 = T3*Tref;
v3 = v3*aref;
m3 = m3*rhoref*Astar*aref;

%% Animations
fig1 = figure(2);
set(fig1, 'Position', [0,0, 1800, 1080]);
for i = 1:length(tvec_realtime1)
    figure(2);
    
    %Profiles
    subplot(5,3,1)
    plot(x/L,lsup1, 'r', 'LineWidth', 2);
    hold on
    plot(x/L,linf1, 'r', 'LineWidth', 2);
    ylim([min(linf1), max(lsup1)]);
    title('Nozzle 1')
    hold off

    subplot(5,3,2)
    plot(x/L,lsup2, 'g', 'LineWidth', 2);
    hold on
    plot(x/L,linf2, 'g', 'LineWidth', 2);
    ylim([min(linf2), max(lsup2)]);
    title('Nozzle 2')
    hold off

    subplot(5,3,3)
    plot(x/L,lsup3, 'b', 'LineWidth', 2);
    hold on
    plot(x/L,linf3, 'b', 'LineWidth', 2);
    ylim([min(linf3), max(lsup3)]);
    title('Nozzle 3')
    hold off

    
    %Nozzle 1
    subplot(5,3,4)
    plot(x, p1(i,:), 'r', 'LineWidth', 2);
    ylim([0, max(max(p1))]);
    title('p1 [atm]');
    grid on  
    
    subplot(5,3,7)
    plot(x, T1(i,:), 'r', 'LineWidth', 2);
    ylim([0, max(max(T1))]);
    title('T1 [K]');
    grid on  
    
    subplot(5,3,10)
    plot(x, Ma1(i,:), 'r', 'LineWidth', 2);
    ylim([0, max(max(Ma1))]);
    title('Ma1');
    grid on  
    
    subplot(5,3,13)
    plot(x, m1(i,:), 'r', 'LineWidth', 2);
    ylim([min(min(m1)), max(max(m1))]);
    title('m1 [kg/s]');
    grid on  
    xlabel('Length [m]');
    
    %Nozzle 2
    subplot(5,3,5)
    plot(x, p2(i,:), 'g', 'LineWidth', 2);
    ylim([0, max(max(p2))]);
    title('p2 [atm]');
    grid on  
    
    subplot(5,3,8)
    plot(x, T2(i,:), 'g', 'LineWidth', 2);
    ylim([0, max(max(T2))]);
    title('T2 [K]');
    grid on  
    
    subplot(5,3,11)
    plot(x, Ma2(i,:), 'g', 'LineWidth', 2);
    ylim([0, max(max(Ma2))]);
    title('Ma2');
    grid on  
    
    subplot(5,3,14)
    plot(x, m2(i,:), 'g', 'LineWidth', 2);
    ylim([min(min(m2)), max(max(m2))]);
    title('m2 [kg/s]');
    grid on  
    xlabel('Length [m]');
    
    %Nozzle 3
    subplot(5,3,6)
    plot(x, p3(i,:), 'b', 'LineWidth', 2);
    ylim([0, max(max(p3))]);
    title('p3 [atm]');
    grid on  
    
    subplot(5,3,9)
    plot(x, T3(i,:), 'b', 'LineWidth', 2);
    ylim([0, max(max(T3))]);
    title('T3 [K]');
    grid on  
    
    subplot(5,3,12)
    plot(x, Ma3(i,:), 'b', 'LineWidth', 2);
    ylim([0, max(max(Ma3))]);
    title('Ma3');
    grid on
    
    subplot(5,3,15)
    plot(x, m3(i,:), 'b', 'LineWidth', 2);
    ylim([min(min(m3)), max(max(m3))]);
    title('m3 [kg/s]');
    grid on
    xlabel('Length [m]');
    
    sgtitle( ['t = ', num2str(tvec_realtime1(i)*1000), ' ms']);
    
    myfram(i) = getframe(fig1);
end
