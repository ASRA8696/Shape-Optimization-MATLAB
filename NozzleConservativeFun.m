%%%%%%%%%%%%%      FULLY EXPANDED FLOW WIHTOUT SHOCK           %%%%%%%%%%%%%%


function [tvec_realtime, rho2, p2, T2, v2, Ma2] = ...
    NozzleConservativeFun(tfin, x, A, pres, pext0, Tres, Text0, gamma, mw_air, R)

%% Parameters

CFLtarget = 0.5; % CFL number to choose time step
dtmax = 1e-2;    % Dimensionless [1]

% Reference state
pref = pres;                   % Reference state pressure [Pa]
Tref = Tres;                   % Reference state temperature [K]

% Derived dimensionless B.C.
Tleft = Tres/Tref;             %Dimensionless temperature
eleft = Tleft;
pleft = pres/pref;             %Dimensionless pressure
rholeft = pleft/Tleft;         %Dimensionless density (from ideal gas law)

% Derived thermo properties
Rair = R/mw_air;               % Mass ideal gas constant for air [J/(kg*K)]
cv = Rair/(gamma-1);           % specific heat of air at constant volume [J/(kg*K)]

% Derived reference state variables
eref = cv*Tref;                % specific energy at reference state [J/kg]
rhoref = pref/(Rair*Tref);     % mass density at reference state [kg/m^3]         
aref = sqrt(gamma*Rair*Tref);  % Speed of sound at reference state [m/s]

% Derived geometrical varaibles 
L = x(end);  %Nozzle length [m]
dx = (x(2)-x(1));  %Use dimensionless dx (with L = 1)
Nx = length(x);


%% Allocate memory
Nt_estim = round(tfin/(dtmax/10));

rho = nan(Nt_estim, Nx);
v = nan(Nt_estim, Nx);
e = nan(Nt_estim, Nx);

%% Initial conditions

%For continuity, use exponential transition between left and right
p0 = ((pres-pext0)*exp(-x/(L/15)) + pext0)/pref;
T0 = ((Tres-Text0)*exp(-x/(L/15)) + Text0)/Tref;
v0 = 1e-3*ones(size(x));
rho0 = p0./T0;
e0 = T0;

rho(1,:) = rho0;
v(1,:) = v0;
e(1,:) = e0;

%% Solve

tval = 0;
n = 1;
rhotemp = rho(1,:);
vtemp = v(1,:);
etemp = e(1,:);

%Variable change for standard conservative form
U1=rhotemp.*A;
U2=rhotemp.*A.*vtemp;
U3=rhotemp.*A.*((etemp/(gamma-1))+(0.5*gamma.*vtemp.*vtemp));

errstat = inf;
dtvec = nan(Nt_estim, 1);
tsimvec = nan(Nt_estim, 1);
tsimvec(1) = 0;
tic;
while tval < tfin
    
     % Copies of solution vectors
     U1_ini=U1;
     U2_ini=U2;
     U3_ini=U3;
    
     % Adaptative time step 
     dt=min(dtmax, min(CFLtarget.*dx./(sqrt(abs(etemp)) + abs(vtemp))));
     dtvec(n) = dt;
     tsimvec(n+1) = tsimvec(n) + dt;
     
     % Defining Flux vectors 
     F1=U2;
     F2=((U2.^2)./(U1))+ ((gamma-1)/gamma).*(U3-(0.5*gamma.*U2.^2)./U1);
     F3=((gamma.*U2.*U3)./U1)-(0.5*gamma*(gamma-1).*U2.^3)./(U1.^2);
     
  
     % Predictor step 
     for i=2:Nx-1
     dF1_dx=(F1(i+1)-F1(i))/dx;
     dF2_dx=(F2(i+1)-F2(i))/dx;
     dF3_dx=(F3(i+1)-F3(i))/dx;
     
     % continuity equation
     dU1_dt_p(i)=-dF1_dx;
     
     % momentum equation
     J2(i) = (1/gamma).*rhotemp(i).*etemp(i).*((A(i+1)-A(i))/dx); %Source term
     dU2_dt_p(i)=-dF2_dx+J2(i);
     
     % energy equation
     dU3_dt_p(i)=-dF3_dx;
     
     %solution update
     U1(i)=U1(i)+dU1_dt_p(i)*dt;
     U2(i)=U2(i)+dU2_dt_p(i)*dt;
     U3(i)=U3(i)+dU3_dt_p(i)*dt;
     end
     
     % Updation of primitive variables and flux vectors
     
     rhotemp = U1./A; 
     vtemp = U2./U1;
     etemp = (gamma-1).*((U3./U1)-((gamma*vtemp.^2)/2));
      
     F1=U2;
     F2=((U2.^2)./(U1))+ ((gamma-1)/gamma).*(U3-(0.5*gamma.*U2.^2)./U1);
     F3=((gamma.*U2.*U3)./U1)-(0.5.*gamma.*(gamma-1).*U2.^3)./(U1.^2);
     
     % Corrector step      
     for j=2:Nx-1
     dF1_dx=(F1(j)-F1(j-1))/dx;
     dF2_dx=(F2(j)-F2(j-1))/dx;
     dF3_dx=(F3(j)-F3(j-1))/dx;
     
     % continuity equation
     dU1_dt_c(j)=-dF1_dx;
     
     % momentum equation
     
     J2(j) = (1/gamma).*rhotemp(j).*etemp(j).*((A(j)-A(j-1))/dx);
     dU2_dt_c(j)=-dF2_dx+J2(j);
     
     % energy equation
     dU3_dt_c(j)=-dF3_dx;
     end
     
     % compute the average time derivative
     dU1_dt_avg=0.5.*(dU1_dt_p+dU1_dt_c);
     dU2_dt_avg=0.5.*(dU2_dt_p+dU2_dt_c);
     dU3_dt_avg=0.5.*(dU3_dt_p+dU3_dt_c);
     
     % Final solution update
     
     for l=2:Nx-1
     U1(l)=U1_ini(l)+dU1_dt_avg(l)*dt;
     U2(l)=U2_ini(l)+dU2_dt_avg(l)*dt;
     U3(l)=U3_ini(l)+dU3_dt_avg(l)*dt;

     
     end

     % Apply boundary conditions 
     % Inlet
     U2(1) = 2*U2(2) - U2(3); %Extrapolation
     % Outlet
     U1(Nx) = 2*U1(Nx-1)-U1(Nx-2);
     U2(Nx) = 2*U2(Nx-1)-U2(Nx-2);
     U3(Nx) = 2*U3(Nx-1)-U3(Nx-2);


     % Updation of flow field variables
     rhotemp = U1./A; 
     vtemp = U2./U1;
     etemp = (gamma-1).*((U3./U1)-((gamma*vtemp.^2)/2));
     
     %Solution for step n+1
     rho(n+1, :) = rhotemp;
     v(n+1, :) = vtemp;
     e(n+1, :) = etemp; 

     errstat = max(abs(v(n+1,:)-v(n,:))) + ...
              max(abs(rho(n+1,:)-rho(n,:))) + ...
                           max(abs(e(n+1,:)-e(n,:)));               

     tval = tval + dt;
     n = n+1;
end
toc;

T = e;
p = rho.*T;
%% Interpolate solution for uniform time-step

Ninterp = 1000;
tvec = linspace(0, tval, Ninterp);
tsimvec = rmmissing(tsimvec);
rho = rmmissing(rho);
p = rmmissing(p);
v = rmmissing(v);
T = rmmissing(T);
Ma = v./sqrt(abs(T));

rho2 = nan(Ninterp,Nx);
p2 = nan(Ninterp,Nx);
v2 = nan(Ninterp,Nx);
T2 = nan(Ninterp,Nx);
Ma2 = nan(Ninterp, Nx);
for i = 1:Nx
    rho2(:,i) = interp1(tsimvec, rho(:,i), tvec);
    p2(:,i) = interp1(tsimvec, p(:,i), tvec);
    v2(:,i) = interp1(tsimvec, v(:,i), tvec);
    T2(:,i) = interp1(tsimvec, T(:,i), tvec);
    Ma2(:,i) = interp1(tsimvec, Ma(:,i), tvec);
end

% Get real time values
tvec_realtime = tvec*L/aref;

end