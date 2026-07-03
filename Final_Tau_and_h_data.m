clc; clear; close all;

% DATA

tau_perpendicular = [217.136 120.951 118.233 113.349 110.214 95.345 90.458 69.776];
v_perpendicular = linspace(0, 3.5, 8);
h_eff_perpendicular = [21.963 39.429 40.233 42.073 43.270 50.018 52.720 68.347];
h_conv_perpendicular = [17.775 37.227 37.337 39.140 39.985 46.888 48.607 64.969];
h_corr_perpendicular = [9.564 22.62 29.16 35.091 40.621 45.54 50.03 54.20];

tau_parallel = [97.34 77.74 63.56 54.07];
v_parallel = [0 1.45 2.5 3.55];
h_eff_parallel = [21.65 27.11 33.16 38.98];
h_conv_parallel = [13.83 21.03 28.08 35.87];
h_corr_parallel = [9.77 13.44 17.35 20.53];

% 1. h vs Velocity (MAIN PLOT)

figure;
plot(v_perpendicular, h_conv_perpendicular, 'ro-', 'LineWidth', 2, 'MarkerSize', 6); hold on;
plot(v_parallel, h_conv_parallel, 'bo-', 'LineWidth', 2, 'MarkerSize', 6);

xlabel('Velocity (m/s)');
ylabel('h (W/m^2K)');
title('Heat Transfer Coefficient vs Velocity');
legend('Perpendicular Flow', 'Parallel Flow', 'Location', 'northwest');
grid on;

% 2. Experimental vs Correlation
figure;

subplot(1,2,1)
plot(v_perpendicular, h_conv_perpendicular, 'ro-', 'LineWidth', 2); hold on;
plot(v_perpendicular, h_corr_perpendicular, 'bo--', 'LineWidth', 2);
xlabel('Velocity (m/s)');
ylabel('h (W/m^2K)');
title('Perpendicular Flow');
legend('Experimental (ODE)', 'Correlation');
grid on;

subplot(1,2,2)
plot(v_parallel, h_conv_parallel, 'ro-', 'LineWidth', 2); hold on;
plot(v_parallel, h_corr_parallel, 'bo--', 'LineWidth', 2);
xlabel('Velocity (m/s)');
ylabel('h (W/m^2K)');
title('Parallel Flow');
legend('Experimental (ODE)', 'Correlation');
grid on;

% 3. Time Constant vs Velocity
figure;
plot(v_perpendicular, tau_perpendicular, 'ro-', 'LineWidth', 2); hold on;
plot(v_parallel, tau_parallel, 'bo-', 'LineWidth', 2);

xlabel('Velocity (m/s)');
ylabel('\tau (s)');
title('Time Constant vs Velocity');
legend('Perpendicular Flow', 'Parallel Flow');
grid on;

% 4. h_eff vs h_conv
figure;

subplot(1,2,1)
plot(v_perpendicular, h_eff_perpendicular, 'go-', 'LineWidth', 2); hold on;
plot(v_perpendicular, h_conv_perpendicular, 'ro-', 'LineWidth', 2);
xlabel('Velocity (m/s)');
ylabel('h (W/m^2K)');
title('Perpendicular Flow');
legend('h_{eff}', 'h_{conv}');
grid on;

subplot(1,2,2)
plot(v_parallel, h_eff_parallel, 'go-', 'LineWidth', 2); hold on;
plot(v_parallel, h_conv_parallel, 'ro-', 'LineWidth', 2);
xlabel('Velocity (m/s)');
ylabel('h (W/m^2K)');
title('Parallel Flow');
legend('h_{eff}', 'h_{conv}');
grid on;

% 5. Ratio Plot (VERY POWERFUL)

% Interpolate parallel data to match perpendicular velocities
h_parallel_interp = interp1(v_parallel, h_conv_parallel, v_perpendicular, 'linear', 'extrap');

ratio = h_conv_perpendicular ./ h_parallel_interp;

figure;
plot(v_perpendicular, ratio, 'ko-', 'LineWidth', 2);

xlabel('Velocity (m/s)');
ylabel('h_{perp} / h_{parallel}');
title('Effect of Flow Orientation');
grid on;
