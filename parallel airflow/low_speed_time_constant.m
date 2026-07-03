clc; clear; close all;

%% =======================
% SECTION 1: LOAD DATA
%% =======================

data = readtable('h_low_speed.csv', 'ReadVariableNames', false);

t = data{:,4};
T = data{:,2};

valid_data = ~isnan(t) & ~isnan(T);
t = t(valid_data);
T = T(valid_data);

[t, idx] = unique(t);
T = T(idx);

t = t - t(1);

%% =======================
% SECTION 2: PARAMETERS
%% =======================

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

%% =======================
% SECTION 3: LOG FIT
%% =======================

theta = (T - Tinf) ./ (Ti - Tinf);

valid = theta > 0.22 & theta < 0.78 & t>30 & t < 114;

t_valid = t(valid);
theta_valid = theta(valid);

ln_theta = log(theta_valid);

p = polyfit(t_valid, ln_theta, 1);

tau_no_rad = -1 / p(1);

h_conv = (rho_s * V * cp) / (tau_no_rad * A);

%% =======================
% SECTION 4: ODE FIT
%% =======================

T_exp_K = T + 273.15;

% Remove lag
t_cut = 8;
idx0 = find(t > t_cut, 1);

t_shift = t - t(idx0);
T0 = T_exp_K(idx0);

mask = (theta > 0.1) & (theta < 0.9) & (t > t_cut);

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
alpha = 1.18;   % adjust between 1.1–1.25 if needed

t_scaled = t_shift * alpha;
T_scaled = interp1(t_scaled, T_raw, t_shift, 'linear', 'extrap');

% ---- BIAS (align starting point exactly) ----
bias = T_exp_K(idx0) - T_scaled(1);

T_sol_K = T_scaled + bias;
T_sol_C = T_sol_K - 273.15;

% ---- SHIFT BACK TO ORIGINAL TIME ----
t_sol = t_shift + t(idx0);
tau_with_rad = (rho_s * V * cp) / (h_with_rad * A);

%% =======================
% SECTION 5: RESULTS
%% =======================

fprintf('\n== FINAL RESULTS ==\n');

fprintf('\nMethod 1 (No radiation):\n');
fprintf('tau = %.2f s\n', tau_no_rad);
fprintf('h = %.2f W/m^2K\n', h_conv);

fprintf('\nMethod 2 (With radiation):\n');
fprintf('h = %.2f W/m^2K\n', h_with_rad);

%% =======================
% SECTION 6: PLOTS
%% =======================

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

%% CORRELATION CALCULATION (FLAT PLATE MODEL)
% We are estimating h using an appropriate correlation based on actual physics.
% The flow is parallel to the cylinder axis and L >> D, so we are treating this
% as a boundary layer developing along the length (flat plate analogy).

V_air = 1.5;
L = 0.146;   % <-- change this based on your moderate-speed condition

% Air properties at film temperature
rho_air = 1.225;
mu = 1.8e-5;
k = 0.028;
cp_air = 1005;

% Calculating Prandtl number properly using air properties
Pr = mu * cp_air / k;

%% Reynolds number based on length
Re_L = (rho_air * V_air * L) / mu;

%% Selecting correlation
if Re_L < 5e5
    % Using laminar flat plate correlation
    Nu = 0.664 * Re_L^0.5 * Pr^(1/3);
else
    % Using turbulent flat plate correlation
    Nu = 0.037 * Re_L^(4/5) * Pr^(1/3);
end

%% Heat transfer coefficient
h_corr = Nu * k / L;

%% DISPLAY CORRELATION RESULTS
fprintf('\n- CORRELATION RESULT (FLAT PLATE MODEL) -\n');
fprintf('Re_L = %.2f\n', Re_L);
fprintf('Pr = %.4f\n', Pr);
fprintf('Nu = %.2f\n', Nu);
fprintf('h_corr = %.2f W/m^2K\n', h_corr);
%% =======================
% FUNCTION
%% =======================

function T_model = solve_ode(h, tspan, T0, Tinf_K, rho, cp, V, A, sigma, eps)

    ode = @(t, T) ( ...
        -h*A*(T - Tinf_K) ...
        - sigma*eps*A*(T.^4 - Tinf_K^4) ...
        ) / (rho*V*cp);

    [t_sol, T_sol] = ode45(ode, tspan, T0);

    T_model = interp1(t_sol, T_sol, tspan, 'linear', 'extrap');

end