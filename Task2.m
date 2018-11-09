%% Task 2 Problem Statement
% Compute power delivered by the PV cell. Compute PV cell conversion
% efficiency. Finally, determine power output. Use the curve fit for I0_IV
% as found in Task 1 to accomplish this. 

%% I0_IV Curve (Task 1)
% Data Input & Linear Curve Fit
data = [28.85, 0.379; 30.492, 0.3737; 43.279, 0.3526; 50.49, 0.3316;...
    60.65, 0.311];
x = data(:,1)+273.15; %convert oC to K
y = data(:,2);

X = [ones(length(x),1), x];

b = X\y;

T_pv = linspace(20+273.15,100+273.15,81);

% Plot Fit
Voc = @(t) b(2).*t + b(1);

%For Kelvin
kb = 1.3807*10^-23; %[J/K]
qe = 1.602*10^-19; %[C]
I0_IV = @(t) (exp(Voc(t).*qe./(kb.*t))-1).^-1;

%% Initializing Constants
pressure = 163; %[kPa]
rc = 15; %concentration ratio, unitless
I_d = 1080; %[W/m^2]
P = rc*I_d; %[W/m^2]
L_x = 0.1; %[m]
L_y = 0.1; %[m]
A = L_x*L_y;

Vg = 1.1; %[V]

T_source=6000; %Temperature of the sun (source)

%% Eq. 14.16: Solving for phi
phi = @(T_source) (pi^4*kb*T_source./(2.404*15*P)).^-1;
% Accepts constant or array of source temperature.
% Returns constant phi value or array of phi vals.

%% Eq. 14.25: Solving for phi_g
sigma = @(x) x.^2./(exp(x)-1);

X = qe*Vg./(kb.*T_source);
integ = integral(sigma, X, Inf);
phi_i = phi(T_source);
phi_g = phi_i*0.416*integ;

IV = phi_g.*qe*A;

%% P = IV 
RL = 0.0070; %[Ohm]
fg = 2.65*10^14; %265 [THz]
RS = 0;
TK1 = 20+273; %K
TK2 = 80+273; %K

syms IL
eqn1 = kb*TK1/qe*log(I0_IV(TK1).^-1-I0_IV(TK1).^-1*IL/IV+1) == IL*RL;
IL_op1 = double(solve(eqn1,IL));

eqn2 = kb*TK2/qe*log(I0_IV(TK2).^-1-I0_IV(TK2).^-1*IL/IV+1) == IL*RL;
IL_op2 = double(solve(eqn2,IL));

VL = @(IL, t) kb*t/qe*log(I0_IV(t).^-1-I0_IV(t).^-1.*IL./IV+1);
VL_op1 = VL(IL_op1,TK1);
VL_op2 = VL(IL_op2,TK2);

%% Plot Resistor and PV - Operating Point
IL_arr = linspace(0,IV, ceil(IV)*10);
plot(IL_arr, IL_arr*RL); %Resistor
hold on;
plot(IL_arr, VL(IL_arr,TK1)); %At operating temp. 1
hold on;
plot(IL_arr, VL(IL_arr,TK2)); %At operating temp. 2
xlabel("I_L [A]");
ylabel("V_L [V]");
%title("VL over IL Performance");
legend("Resistor","PV Cell @ 20 ^oC", "PV Cell @ 80 ^oC","Location","northwest");

%% Power and Efficiency Results

P1 = IL_op1*VL_op1;
eff1 = P1/(P*A);
P2 = IL_op2*VL_op2;
eff2 = P2/(P*A);
fprintf("At 20 oC...\n" ...
    + "    Power = " + P1 + " W \n" ...
    + "    Operating point: IL = " + IL_op1 + " A , VL = " + VL_op1 + " V \n"...
    + "    Efficiency = " + eff1 + "\n\n");
fprintf("At 80 oC...\n" ...
    + "    Power = " + P2 + " W \n" ...
    + "    Operating point: IL = " + IL_op2 + " A , VL = " + VL_op2 + " V \n"...
    + "    Efficiency = " + eff2 + "\n");