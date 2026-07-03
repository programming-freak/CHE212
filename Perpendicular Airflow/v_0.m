clc; clear;

%% Loading the Data
filename = 'v = 0.txt';
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

%  TEMPERATURE (K)

T_exp = T + 273;
T0 = T_exp(1);
Tinf = 32.5;
Tinf_K = Tinf + 273.15;
theta = (T_exp - Tinf_K) ./ (T0 - Tinf_K);

%  MATERIAL + GEOMETRY

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

%  LINEAR FIT 

% Valid theta
valid_theta = theta > 0;

t_lin = t(valid_theta);
theta_lin = theta(valid_theta);

ln_theta = log(theta_lin);
[theta_unique, idx] = unique(theta_lin, 'stable');

t_unique = t_lin(idx);
ln_theta_unique = log(theta_unique);
mask = t_unique < 400;  

t_fit = t_unique(mask);
ln_fit = ln_theta_unique(mask);

% Linear regression
p = polyfit(t_fit, ln_fit, 1);

slope = p(1);
tau = -1/slope;

% h effective
h_eff = (m * cp) / (A * tau);

fprintf('\n= LINEARIZED MODEL =\n');
fprintf('tau = %.3f s\n', tau);
fprintf('h_eff = %.3f W/m^2-K\n', h_eff);

%  NONLINEAR ODE FIT 

sigma = 5.67e-8;
eps = 0.75;

% THETA-BASED CUTOFF
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

%%  PLOT 2: T vs time 

T_fit = model(t);

figure;

plot(t, T_exp, 'bo', 'DisplayName', 'Data'); hold on;
plot(t, T_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'ODE Fit');

xlabel('Time (s)');
ylabel('Temperature (K)');
legend;
grid on;
title('Temperature vs Time');

%% =======================
% NATURAL CONVECTION (v = 0)
%% =======================

% Air properties (film temperature ~320 K approx)
rho = 1.225;
mu = 1.9e-5;
k = 0.028;
Pr = 0.71;

nu = mu / rho;

% Temperatures (use average from your data)
T_s = mean(T_exp(1:20));     % initial hot temp
T_inf = 32.5 + 273.15;

T_film = (T_s + T_inf)/2;

% Thermal expansion coefficient
beta = 1 / T_film;

% Geometry
D = A/((L1 + L2 + L3 + r3)*pi) ;  % effective diameter

% Gravity
g = 9.81;

% Grashof
Gr = g * beta * (T_s - T_inf) * D^3 / nu^2;

% Rayleigh
Ra = Gr * Pr;

% Churchill-Chu
Nu = (0.60 + (0.387 * Ra^(1/6)) / ...
    (1 + (0.559/Pr)^(9/16))^(8/27))^2;

% h
h_nat = Nu * k / D;

fprintf('\n NATURAL CONVECTION \n');
fprintf('Ra = %.3e\n', Ra);
fprintf('Nu = %.3f\n', Nu);
fprintf('h_nat = %.3f W/m^2-K\n', h_nat);

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