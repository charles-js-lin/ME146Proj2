%% Task 5 Problem Statement
% a) Iterate over a range of ID from [720, 1080]
% b) Iterate over a range of rc from [10, 18]

tic;
tol = 1e-5;
%% Initialize Constants
results = [];
ID_range = linspace(720, 1080, 10);
rc_range = linspace(10, 18, 9);
% Iterate over 10 points for 5(a), 9 points for 5(b)
for h=1:10
T_pv = [170, 0, 240];
% 5(a)
ID = ID_range(h);%[W/m^2]
rc = 15; %unitless

% 5(b)
% ID = 1080;
% rc = rc_range(h);

a = 0.0017; %[V/K], from Task 3
lA = 0.032; %[W/(cm K)], from Task 3
lB = 0.021; %[W/(cm K)], from Task 3
pA = 0.0020; %[ohm cm], from Task 3
pB = 0.0030; %[ohm cm], from Task 3
A_cont = 300; %[W/K]
Rs = 0; %[ohm]
Vg = 1.1; %[V]
Ly = 0.1; %[m] from 10 cm
Lz = 0.1; %[m] from 10 cm
A_pv = Ly*Lz; %[m^2]
%Tw_in = 0.004; %[m] from 4 mm, Task 3
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
        else
            TC(1) = TC(2);
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
    if i>=90
        break
    end
end

syms IL_te
heatRateEqn = QH_result == n*a*IL_te*TH(2) + A_batt*(TH(2)-TC_result)-0.5*IL_te^2*R_batt(2);

IL_result = double(solve(heatRateEqn, IL_te)); 
VL_result = n*a*(TH(2)-TC_result)-R_batt(2)*IL_result; 
VL_true = VL_result(1); %[V]
IL_true = IL_result(1); %[A]

%% Results
P_tot = PL_result + Wdot_result;
eff_tot = P_tot/(P*A_pv);

% 5(a)
col = [ID; eff_pv_result; eff_result; eff_tot; P_tot];
disp(ID);

% 5(b)
% col = [rc; eff_pv_result; eff_result; eff_tot; P_tot];
% disp(rc);

results = [results col];

end
%% Plotting 5(a): Variation in ID
results_ID = [720	760	800	840	880	920	960	1000	1040	1080;...
0.120243207	0.126923381	0.133603536	0.140283589	0.146963128	0.153640291	0.160307132	0.166931865	0.173396299	0.179304313;...
0.142863119	0.147008975	0.150958348	0.154724322	0.15831879	0.161752576	0.165036272	0.168182356	0.171205564	0.174183114;...
0.245927994	0.255273475	0.264393325	0.273302617	0.282014897	0.290541145	0.298886905	0.307039218	0.314904727	0.322255628;...
26.5602234	29.10117612	31.72719905	34.43612979	37.22596642	40.09467807	43.03971426	46.05588263	49.12513744	52.20541181];

subplot(3,1,1);
plot(results(1,:),results(2,:));
ylabel("PV Eff.");
title("Component & Total Efficiency over ID (720 to 1080)");
subplot(3,1,2);
plot(results(1,:),results(3,:));
ylabel("TE Eff.");
subplot(3,1,3);
plot(results(1,:),results(4,:));
xlabel("ID [W/m^2]");
ylabel("Total Eff.");

figure;
plot(results(1,:),results(5,:));
xlabel("ID [W/m^2]");
ylabel("Total Power [W]");
title("Total Power over ID (720 to 1080)");

%% Plotting 5(b): Variation in rc
% 
% results_rc = [10	11	12	13	14	15	16	17	18;...
% 0.120243207	0.132267509	0.144291435	0.156309166	0.16824432	0.179304313	0.184836665	0.178849182	0.165745873;...
% 0.142863119	0.150183531	0.156900869	0.16308345	0.168796666	0.174183114	0.179856849	0.186624746	0.193995284;...
% 0.245927994	0.262586628	0.278552841	0.293901167	0.308641912	0.322255628	0.33144936	0.332096236	0.327574222;...
% 26.5602234	31.19529139	36.10044822	41.26372385	46.66665708	52.20541181	57.27444937	60.97286891	63.68042874];
% 
% 
% subplot(3,1,1);
% plot(results(1,:),results(2,:));
% ylabel("PV Eff.");
% title("Component & Total Efficiency over rc (10 to 18)");
% subplot(3,1,2);
% plot(results(1,:),results(3,:));
% ylabel("TE Eff.");
% subplot(3,1,3);
% plot(results(1,:),results(4,:));
% xlabel("rc");
% ylabel("Total Eff.");
% 
% figure;
% plot(results(1,:),results(5,:));
% xlabel("rc");
% ylabel("Total Power [W]");
% title("Total Power over rc (10 to 18)");

toc; %Find program runtime