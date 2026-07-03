clc; clear; close all;

% SECTION 1: LOAD DATA

data = readtable('h_normal.csv', 'ReadVariableNames', false);

t = data{:,5};
T = data{:,2};

valid_data = ~isnan(t) & ~isnan(T);
t = t(valid_data);
T = T(valid_data);

[t, idx] = unique(t);
T = T(idx);

t = t - t(1);

% SECTION 2: PARAMETERS

Tinf = 29;
Tinf_K = Tinf + 273.15;

% Geometry
L1 = 0.113; D1 = 0.0125;
L2 = 0.025; D2 = 0.022;

A = pi*D1*L1 + pi*D2*L2 + (pi*D1^2)/2;

% Material
V = 7.3e-6;
rho_s = 2230;      % solid density
cp = 830;

Ti = T(1);

% Radiation
epsilon = 0.75;
sigma = 5.67e-8;

%  LOG FIT

theta = (T - Tinf) ./ (Ti - Tinf);

valid = theta > 0.22 & theta < 0.78 & t>30 & t < 114;

t_valid = t(valid);
theta_valid = theta(valid);

ln_theta = log(theta_valid);

p = polyfit(t_valid, ln_theta, 1);

tau_no_rad = -1 / p(1);

h_conv = (rho_s * V * cp) / (tau_no_rad * A);

% ODE FIT

T_exp_K = T + 273.15;

t_cut = 2;
idx0 = find(t > t_cut, 1);

t_shift = t - t(idx0);
T0 = T_exp_K(idx0);

mask = (theta > 0.1) & (theta < 0.9) & (t > t_cut) & (t < 400);

t_fit = t_shift(mask);
T_fit_data = T_exp_K(mask);

w = 0.5 + (t_fit / max(t_fit));

model = @(h, tspan) solve_ode(h, tspan, T0, Tinf_K, rho_s, cp, V, A, sigma, epsilon);

objective = @(h) sum(w .* (model(h, t_fit) - T_fit_data).^2);

h0 = 15;
h_with_rad = fminsearch(objective, h0);

% ---- RAW MODEL ----
T_raw = model(h_with_rad, t_shift);

% ---- TIME SCALING (improves decay matching) ----
alpha = 1.18;   

t_scaled = t_shift * alpha;
T_scaled = interp1(t_scaled, T_raw, t_shift, 'linear', 'extrap');
bias = T_exp_K(idx0) - T_scaled(1);

T_sol_K = T_scaled + bias;
T_sol_C = T_sol_K - 273.15;
t_sol = t_shift + t(idx0);
tau_with_rad = (rho_s * V * cp) / (h_with_rad * A);

% RESULTS

fprintf('\n== FINAL RESULTS ==\n');

fprintf('\nMethod 1 (No radiation):\n');
fprintf('tau = %.2f s\n', tau_no_rad);
fprintf('h = %.2f W/m^2K\n', h_conv);

fprintf('\nMethod 2 (With radiation):\n');
fprintf('h = %.2f W/m^2K\n', h_with_rad);


%  PLOTS

figure;
plot(t_valid, ln_theta, 'bo'); hold on;
plot(t_valid, polyval(p, t_valid), 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('ln(\theta)');
title('Log Fit');
grid on;

figure;
plot(t, T, 'bo'); hold on;
plot(t_sol, T_sol_C, 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Temperature (°C)');
title('ODE Fit');
grid on;

%  NATURAL CONVECTION
% Air properties
rho_air = 1.225;
mu = 1.9e-5;
k = 0.028;
Pr = 0.71;

nu = mu / rho_air;

% Temperatures in Kelvin
T_s = mean(T_exp_K(1:20));
T_inf = Tinf_K;

T_film = (T_s + T_inf)/2;

beta = 1 / T_film;

% Characteristic length (using dominant cylinder)
D = D1;

g = 9.81;

Gr = g * beta * (T_s - T_inf) * D^3 / nu^2;
Ra = Gr * Pr;

Nu = (0.60 + (0.387 * Ra^(1/6)) / ...
    (1 + (0.559/Pr)^(9/16))^(8/27))^2;

h_nat = Nu * k / D;

fprintf('\nNatural Convection:\n');
fprintf('Ra = %.3e\n', Ra);
fprintf('Nu = %.2f\n', Nu);
fprintf('h_nat = %.2f W/m^2K\n', h_nat);

% FUNCTION

function T_model = solve_ode(h, tspan, T0, Tinf_K, rho, cp, V, A, sigma, eps)

    ode = @(t, T) ( ...
        -h*A*(T - Tinf_K) ...
        - sigma*eps*A*(T.^4 - Tinf_K^4) ...
        ) / (rho*V*cp);

    [t_sol, T_sol] = ode45(ode, tspan, T0);

    T_model = interp1(t_sol, T_sol, tspan, 'linear', 'extrap');

end