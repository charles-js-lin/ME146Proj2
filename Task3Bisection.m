%% Task 3 Problem Statement
% Assume: n = 12, R_L,te = 0.10 ohms, TH = 50 oC. All necessary constants
%           given.
% Find: TC, power output of the module Wdot, QHdot, QCdot, current to the
%           load IL, voltage across the load VL
tic
%% Assumed Constants and Initial Calcs
tol = 1e-5;
a = 0.0017; %[V/K], mean Seebeck coefficient
pA = 0.0020; %[ohm*cm]
pB = 0.0030; %[ohm*cm]
lA = 0.032; %[W/(cm*K)]
lB = 0.021; %[W/(cm*K)]
A_pv = 0.01; %[m^2]
k_wEff = 6.5; %[W/(mK)]
Tsat = 95.0; %[K]
Uev = 25; %[W/(m^2*K)]
tw = 0.004; %[m], 4 mm

RA_min = (sqrt(lA*pA)+sqrt(lB*pB))^2;
Z = a^2/RA_min;

n = 12; %12 thermoelectric pairs
R_lte = 0.10; %[ohms]
TH = 50 + 273.15; %[K]

%% TC and QHdot Computation
error = [1000 1000 1000]; %initialize error
TC = [90 0 120]; %LB Midpoint UB

while abs(error(2))>tol || (error(1)*error(3)<0)
    midpoint = (TC(1)+TC(3))/2;
    TC(2) = midpoint;
    
    m_maxEff = sqrt(1+0.5.*(TH+TC).*Z);
    R = R_lte./(n.*m_maxEff);
    A = RA_min./R;

    % (ii) Computing QC_boil
    QC_boil = (k_wEff/tw + Uev)*(TC-Tsat)*A_pv; %[W]

    % (iii) Compute efficiency using Eq. 11
    eff = (TH-TC)./TH.*((1+m_maxEff).^2./(m_maxEff.*Z*TH) + 1 + 1./(2*m_maxEff).*(1+TC./TH)).^-1;

    % (iv)  Computer power output Wdot using Eq. 10
    R_batt = n*R;
    Wdot = (n^2*a^2.*(TH-TC).^2.*R_lte)./(R_lte+R_batt).^2; %[W]

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

%% IL and VL Computation
% Use Eqs. (8) and (9)

syms IL
heatRateEqn = QH_result == n*a*IL*TH + A_batt*(TH-TC_result)-0.5*IL^2*R_batt(2);
%We get 2 roots for IL because heatRateEqn is 2nd order eqn
IL_result = double(solve(heatRateEqn, IL)); 

%Solve for VL using both IL roots
%Select true VL and IL based on most reasonable VL_result
VL_result = n*a*(TH-TC_result)-R_batt(2)*IL_result; 
VL_true = VL_result(1); %[V]
IL_true = IL_result(1); %[A]

%% Results
fprintf("TC = " + TC_result + " K \n" ...
    + "Power output of module Wdot = " + Wdot_result + " W \n" ...
    + "Heat rate QH = " + QH_result + " W \n" ...
    + "Contact surface heat rate QC = " + QC_result + " W \n" ...
    + "Current to the load IL = " + IL_true + " A \n" ...
    + "Voltage across the load VL = " + VL_true + " V \n \n" ...
    + "Error in QC = " + error_result + "\n");
  toc