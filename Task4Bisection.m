%% Task 4 Problem Statement
% Task 2 finds power and performance of PV cell
% Task 3 finds power and params for TE module
% Combine programs from Task 2 and 3 in order to determine operating point
% and performance for HYBRID system (PV + TE). 
tic;
tol = 1e-5;
%% Initialize Constants

T_pv = [170, 0, 230];

ID = 1080; %[W/m^2]
a = 0.0017; %[V/K], from Task 3
lA = 0.032; %[W/(cm K)], from Task 3
lB = 0.021; %[W/(cm K)], from Task 3
pA = 0.0020; %[ohm cm], from Task 3
pB = 0.0030; %[ohm cm], from Task 3
A_cont = 300; %[W/K]
Rs = 0; %[ohm]
Vg = 1.1; %[V]
rc = 15; %unitless
Ly = 0.1; %[m] from 10 cm
Lz = 0.1; %[m] from 10 cm
A_pv = Ly*Lz; %[m^2]
RL_pv = 0.0070; %[ohm]
RL_te = 0.10; %[ohm]

T_source=6000; %Temperature of the sun (source)
P = rc*ID; %[W/m^2]

Tsat = 90.2; %[K]
n = 12; %unitless
Uev = 25; %[W/(m^2*K)]
tw = 0.004; %[m], 4 mm
k_wEff = 6.5; %[W/(mK)]

%% I0_IV Curve (Task 1 for PV Cell Task 2)
% Data Input & Linear Curve Fit
data = [28.85, 0.379; 30.492, 0.3737; 43.279, 0.3526; 50.49, 0.3316;...
    60.65, 0.311];
x = data(:,1)+273.15; %convert oC to K
y = data(:,2);
X = [ones(length(x),1), x];
b = X\y;
% Plot Fit
Voc = @(t) b(2).*t + b(1);
%For Kelvin
kb = 1.3807*10^-23; %[J/K]
qe = 1.602*10^-19; %[C]
I0_IV = @(t) (exp(Voc(t).*qe./(kb.*t))-1).^-1;

%% Step 2: Task 2 (PV Cell Analysis)
phi = @(T_source) (pi^4*kb*T_source./(2.404*15*P)).^-1;
sigma = @(x) x.^2./(exp(x)-1);
X = qe*Vg./(kb.*T_source);
integ = integral(sigma, X, Inf);
phi_i = phi(T_source);
phi_g = phi_i*0.416*integ;
IV = phi_g.*qe*A_pv;
i = 1;

errorEQ = [1000 1000 1000];

while abs(errorEQ(2))>tol || (errorEQ(1)*errorEQ(3)<0)
    midT_pv = (T_pv(1)+T_pv(3))/2;
    T_pv(2) = midT_pv;
    
    syms IL
    eqn1 = kb.*T_pv(1)./qe.*log(I0_IV(T_pv(1)).^-1-I0_IV(T_pv(1)).^-1.*IL./IV+1) == IL.*RL_pv;
    eqn2 = kb.*T_pv(2)./qe.*log(I0_IV(T_pv(2)).^-1-I0_IV(T_pv(2)).^-1.*IL./IV+1) == IL.*RL_pv;
    eqn3 = kb.*T_pv(3)./qe.*log(I0_IV(T_pv(3)).^-1-I0_IV(T_pv(3)).^-1.*IL./IV+1) == IL.*RL_pv;
    IL_op(1) = double(solve(eqn1,IL));
    IL_op(2) = double(solve(eqn2,IL));
    IL_op(3) = double(solve(eqn3,IL));

    VL = @(IL, t) kb*t/qe.*log(I0_IV(t).^-1-I0_IV(t).^-1.*IL./IV+1);
    VL_op = VL(IL_op,T_pv);

    PL = IL_op.*VL_op; %DELIVERABLE: POWER OUTPUT
    eff_pv = PL/(P*A_pv);

    Q_pv = P*A_pv-PL; %DELIVERABLE: WASTE HEAT TRANSFER RATE
    Qplot(i) = Q_pv(2);
%% Step 3: Task 3 (TE Module Analysis)

    TH = T_pv(2) - Q_pv/A_cont;
    RA_min = (sqrt(lA*pA)+sqrt(lB*pB))^2;
    Z = a^2/RA_min;

    error = [1000 1000 1000]; %initialize error
    TC = [90 0 120]; %LB Midpoint UB

    while abs(error(2))>tol || (error(1)*error(2)<0)
        midTC = (TC(1)+TC(3))/2;
        TC(2) = midTC;

        m_maxEff = sqrt(1+0.5.*(TH+TC).*Z);
        R = RL_te./(n.*m_maxEff);
        A = RA_min./R;

        % (ii) Computing QC_boil
        QC_boil = (k_wEff/tw + Uev)*(TC-Tsat)*A_pv; %[W]

        % (iii) Compute efficiency using Eq. 11
        eff = (TH-TC)./TH.*((1+m_maxEff).^2./(m_maxEff.*Z.*TH) + 1 + 1./(2*m_maxEff).*(1+TC./TH)).^-1;

        % (iv)  Computer power output Wdot using Eq. 10
        R_batt = n*R;
        Wdot = (n^2*a^2.*(TH-TC).^2.*RL_te)./(RL_te+R_batt).^2; %[W]

        % (v) Compute heat input and rejection
        QH = Wdot./eff; %[W]
        QC_TE = QH - Wdot; %[W]

        % (vi) Computer error - run if statement to find optimal TC and QH for
        %       minimal error
        error = QC_TE-QC_boil;
        if error(1)*error(2)<0
            TC(3) = TC(2);
        elseif error(1)*error(2)>0 && error(1)>0
            TC(1) = TC(2);
        elseif error(1)*error(2)>0 && error(1)<0
            TC(3) = TC(2);
        end
    end

    error_result = error(2);
    TC_result = TC(2);
    QH_result = QH(2);
    QC_result = QC_boil(2);
    A_batt = n*A(2);
    Wdot_result = Wdot(2);
    eff_result = eff(2);
    QHplot(i) = QH_result;
%% Step 4: Compute Error EQ
% Power, Operating Point, PV cell efficiency
    
    errorEQ = Q_pv-QH;%_result;
    if errorEQ(1)*errorEQ(2)<0
        T_pv(3) = T_pv(2);
    elseif errorEQ(1)*errorEQ(2)>0 && errorEQ(1)>0
        T_pv(1) = T_pv(2);
    elseif errorEQ(1)*errorEQ(2)>0 && errorEQ(1)<0
        T_pv(3) = T_pv(2);
    end
    
    PL_result = PL(2);
    Q_pv_result = Q_pv(2);
    eff_pv_result = eff_pv(2);
    errorEQ_result = errorEQ(2);
    IL_pv_result = IL_op(2);
    VL_pv_result = VL_op(2);
    T_pv_final = T_pv(2);
    
    T_pvplot(i) = T_pv_final;
    i= i +1;
end

syms IL_te
heatRateEqn = QH_result == n*a*IL_te*TH(2) + A_batt*(TH(2)-TC_result)-0.5*IL_te^2*R_batt(2);

IL_result = double(solve(heatRateEqn, IL_te)); 
VL_result = n*a*(TH(2)-TC_result)-R_batt(2)*IL_result; 
VL_true = VL_result(1); %[V]
IL_true = IL_result(1); %[A]

%% Plots
%Error plot 
plot(T_pvplot, Qplot)
hold on;
plot(T_pvplot, QHplot);
legend("Q_p_v","Q_H_,_t_e",'Location','southeast');
xlabel("T_p_v [K]");
ylabel("Heat Transfer [W]");

%% Results
P_tot = PL_result + Wdot_result;
eff_tot = P_tot/(P*A_pv);

fprintf("PV Cell Results: \n" ...
    + "    Operating Point: VL = " + VL_pv_result + " V, IL = " + IL_pv_result  + " A \n" ...
    + "    Operating Point Temperature: T_pv = " + T_pv(2) + " K \n" ...
    + "    Power = " + PL_result + " W \n" ...
    + "    Efficiency = " + eff_pv_result + "\n \n");

fprintf("Thermoelectric Module Results: \n" ...
    + "    Operating Point: VL = " + VL_op(2) + " V, IL = " + IL_op(2)  + " A \n" ...
    + "    Power Output: Wdot = " + Wdot_result + " W \n" ...
    + "    Efficiency = " + eff_result + "\n \n");

fprintf("Combined power output of hybrid device: \n" ...
    + "    Total Power Output = " + (P_tot) + " W \n" ...
    + "    Efficiency = " + eff_tot + "\n \n");

toc; %Find program runtime