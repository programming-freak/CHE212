clc; clear;

%% Loading the Data
filename = 'v = 1.5.txt';
fid = fopen(filename, 'r');
data = textscan(fid, '%s %s', 'Delimiter', '\t');
fclose(fid);
time_str = data{1};
T_str = data{2};
T = str2double(erase(T_str, ','));
valid = ~isnan(T);
time_str = time_str(valid);
T = T(valid);
time_dt = datetime(time_str, 'InputFormat', 'HH:mm:ss.SSS');
t = seconds(time_dt - time_dt(1));
t = double(t);

% TEMPERATURE (K)

T_exp = T + 273;
T0 = T_exp(1);
Tinf = 32.6;
Tinf_K = Tinf + 273.15;
theta = (T_exp - Tinf_K) ./ (T0 - Tinf_K);

% MATERIAL + GEOMETRY

rho = 2230;
cp = 830;
V = 17.4e-6;
m = rho * V;

% Geometry
d1 = 0.017; L1 = 0.013;
d2 = 0.022; L2 = 0.025;
d3 = 0.0113;
r3 = d3/2;

L3 = 0.122 - r3;

A = pi*d1*L1 + pi*d2*L2 + pi*d3*L3 + 2*pi*r3^2;

% LINEAR FIT 

% Valid theta
valid_theta = theta > 0;

t_lin = t(valid_theta);
theta_lin = theta(valid_theta);

ln_theta = log(theta_lin);

%  LINEAR FIT 

valid_theta = theta > 0;

t_lin = t(valid_theta);
theta_lin = theta(valid_theta);

% Remove duplicates
[theta_unique, idx] = unique(theta_lin, 'stable');

t_unique = t_lin(idx);
ln_theta_unique = log(theta_unique);

% REMOVE LAST 3 POINTS
n_remove = 2;

t_fit = t_unique(2:end-n_remove);
ln_fit = ln_theta_unique(2:end-n_remove);

% Linear regression
p = polyfit(t_fit, ln_fit, 1);

tau = -1/p(1);
h_eff = (m * cp) / (A * tau);

fprintf('\n= LINEARIZED MODEL =\n');
fprintf('tau = %.3f s\n', tau);
fprintf('h_eff = %.3f W/m^2-K\n', h_eff);


%  NONLINEAR ODE FIT 

sigma = 5.67e-8;
eps = 0.75;

mask = theta > 0.05;  

t_fit_ode = t(mask);
T_fit_data = T_exp(mask);

w = 1 + 3*(t_fit_ode / max(t_fit_ode)).^2;

% Model
model = @(h, tspan) solve_ode(h, tspan, T0, Tinf_K, m, cp, A, sigma, eps);

% Objective with weights
objective = @(h) sum( w .* (model(h, t_fit_ode) - T_fit_data).^2 );

h0 = 10;

h_conv = fminsearch(objective, h0);

fprintf('\n= ODE MODEL =\n');
fprintf('h_conv = %.3f W/m^2-K\n', h_conv);% Model handle with fitted Tinf
model = @(tspan) solve_ode(h_conv, tspan, T0, Tinf_K, m, cp, A, sigma, eps);

% PLOT 1: ln(theta) vs time 

figure;

plot(t_unique, ln_theta_unique, 'bo', 'DisplayName', 'Unique Data'); hold on;

% Only plot fit region
plot(t_fit, polyval(p, t_fit), 'r-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');

xlabel('Time (s)');
ylabel('ln(\theta)');
legend;
grid on;
title('ln(\theta) vs Time');

% PLOT 2: T vs time 

T_fit = model(t);

figure;

plot(t, T_exp, 'bo', 'DisplayName', 'Data'); hold on;
plot(t, T_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'ODE Fit');

xlabel('Time (s)');
ylabel('Temperature (K)');
legend;
grid on;
title('Temperature vs Time');

% THEORETICAL h (CB cylinders + Ranz-Marshall tip)

% AIR PROPERTIES
rho_air = 1.225;        % kg/m^3
mu_air  = 1.85e-5;     % Pa.s
k_air   = 0.026;       % W/m-K
Pr      = 0.71;

V = 1.5;   

% GEOMETRY 
% Cylinders
d1 = 0.017;  L1 = 0.013;
d2 = 0.022;  L2 = 0.025;
d3 = 0.0113; 
r3 = d3/2;
L3 = 0.122 - r3;   % cylinder part excluding hemisphere

% AREAS
A1 = pi*d1*L1;
A2 = pi*d2*L2;
A3 = pi*d3*L3;
A_tip = 2*pi*r3^2;

% 1. CYLINDERS → Churchill-Bernstein

% Reynolds numbers
Re1 = rho_air * V * d1 / mu_air;
Re2 = rho_air * V * d2 / mu_air;
Re3 = rho_air * V * d3 / mu_air;

% Churchill-Bernstein function
CB = @(Re) 0.3 + ...
    (0.62*sqrt(Re)*Pr^(1/3)) / (1 + (0.4/Pr)^(2/3))^(1/4) * ...
    (1 + (Re/282000)^(5/8))^(4/5);

Nu1 = CB(Re1);
Nu2 = CB(Re2);
Nu3 = CB(Re3);

h1 = Nu1 * k_air / d1;
h2 = Nu2 * k_air / d2;
h3 = Nu3 * k_air / d3;

% HEMISPHERE → Ranz-Marshall (sphere)

Re_tip = rho_air * V * d3 / mu_air;

Nu_tip = 2 + 0.6 * sqrt(Re_tip) * Pr^(1/3);

h_tip = Nu_tip * k_air / d3;

% AREA-WEIGHTED AVERAGE

A_total = A1 + A2 + A3 + A_tip;

h_avg = (h1*A1 + h2*A2 + h3*A3 + h_tip*A_tip) / A_total;

% OUTPUT

fprintf('\n===== THEORETICAL h (MULTI-SECTION) =====\n');
fprintf('h1 = %.2f W/m^2-K\n', h1);
fprintf('h2 = %.2f W/m^2-K\n', h2);
fprintf('h3 = %.2f W/m^2-K\n', h3);
fprintf('h_tip = %.2f W/m^2-K\n', h_tip);

fprintf('\n>>> AREA-WEIGHTED h_avg = %.2f W/m^2-K\n', h_avg);

% FUNCTION

function T_model = solve_ode(h, tspan, T0, Tinf_K, m, cp, A, sigma, eps)

    ode = @(t, T) ( ...
        h*A*(Tinf_K - T) .* (abs(T - Tinf_K)/50).^0.25 + ...
        sigma*eps*A*(Tinf_K^4 - T.^4) ...
        ) / (m*cp);

    t_range = [tspan(1), tspan(end)];

    [t_sol, T_sol] = ode45(ode, t_range, T0);

    T_model = interp1(t_sol, T_sol, tspan, 'linear');

end