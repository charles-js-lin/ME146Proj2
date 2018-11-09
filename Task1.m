%% Data Input & Linear Curve Fit
data = [28.85, 0.379; 30.492, 0.3737; 43.279, 0.3526; 50.49, 0.3316;...
    60.65, 0.311];
x = data(:,1)+273.15; %convert oC to K
y = data(:,2);

X = [ones(length(x),1), x];

b = X\y;

T_pv = linspace(20+273.15,100+273.15,81);

%% Plot Fit
Voc = @(t) b(2).*t + b(1);
plot(T_pv, Voc(T_pv)); 

hold on;
plot(x,y,'o');
xlabel("Temperature [K]");
ylabel("V_o_c [V]");
title("Linear Fit");

%% Plot I0_IV
kb = 1.3807*10^-23; %[J/K]
qe = 1.602*10^-19; %[C]

%For Kelvin
I0_IV = @(t) (exp(Voc(t).*qe./(kb.*t))-1).^-1;
figure;
plot(T_pv, I0_IV(T_pv));
xlabel("Temperature [K]");
ylabel("I_0/I_V Ratio");
title("I0\_IV Ratio over Temp.");