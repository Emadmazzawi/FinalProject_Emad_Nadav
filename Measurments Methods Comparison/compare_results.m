% =========================================================================
% SCRIPT: Compare Multiple Measurement Days / Conditions (Bandgap Only)
% Compares Bandgap (Tauc) against Jain (1992) Theory down to -60 C
% Includes Hardcoded PL measurements (Converted from nm to eV).
% Comments: English Only
% =========================================================================
clear all; close all; clc;

%% --- 1. DEFINE FILES TO COMPARE ---
data_files = {'results_vacuum1.mat', 'results_nitro.mat', 'non_results.mat', 'results_nitro_vacuum.mat'};
run_labels = {'Vacuum 1', 'Nitrogen', 'Cooling Only', 'Vacuum + Nitro'};

% Colors exactly matching the noise analysis: Green, Red, Blue, Purple
colors  = {[0 0.7 0], [1 0 0], [0 0.4 1], [0.6 0 0.8]}; 
markers = {'o', '^', 'v', 'd'};
num_runs = length(data_files);

%% --- 2. INITIALIZE COMPARISON FIGURE ---
h_fig_eg = figure('Color', 'w', 'Name', 'Bandgap Comparison', 'units', 'normalized', 'outerposition', [0.2 0.1 0.6 0.7]);
hold on; grid on; box on;
set(gca, 'FontSize', 16); % Increased axis font for poster
xlabel('Temperature (°C)', 'FontSize', 20, 'FontWeight', 'bold'); 
ylabel('Band Gap Energy (eV)', 'FontSize', 20, 'FontWeight', 'bold'); 
title('GaAs Bandgap: Theory vs. Ellipsometry & PL (Up to -60°C)', 'FontSize', 24, 'FontWeight', 'bold');

%% --- 3. LOAD ELLIPSOMETRY DATA AND PLOT ---
loaded_n_doping = 3.8e18; % Default fallback value

for i = 1:num_runs
    if isfile(data_files{i})
        data = load(data_files{i});
        T = data.temps_C;
        Eg = data.extracted_Eg;
        if isfield(data, 'n_doping'), loaded_n_doping = data.n_doping; end
        
        valid_eg = ~isnan(Eg) & (T >= -60) & (T <= 25);
        
        set(0, 'CurrentFigure', h_fig_eg);
        if any(valid_eg)
            plot(T(valid_eg), Eg(valid_eg), ['-' markers{i}], 'Color', colors{i}, 'MarkerFaceColor', colors{i}, ...
                'MarkerSize', 10, 'DisplayName', run_labels{i}, 'LineWidth', 2);
        end
    else
        warning('File not found: %s', data_files{i});
    end
end

%% --- 4. HARDCODED PL MEASUREMENTS ---
% Actual experimental PL peak values from lab measurements
T_pl  = [20, 10, 0, -10, -20, -30, -40, -50, -60, -70, -80, -90, -100];
meas_Pk = [851.32, 849.05, 846.94, 844.86, 842.93, 840.99, 839.20, 837.55, 835.92, 834.4, 832.8, 831.1, 829.4];

% Convert PL peak wavelength (nm) to Energy (eV) using E = 1240 / lambda
Eg_pl = 1240 ./ meas_Pk;

% Ensure they are column vectors
T_pl = T_pl(:);
Eg_pl = Eg_pl(:);

% Filter to match the plot range (down to -60C)
valid_pl = ~isnan(Eg_pl) & (T_pl >= -60) & (T_pl <= 25);

set(0, 'CurrentFigure', h_fig_eg);
if any(valid_pl)
    % Plot PL with a distinct Pentagram (Star) marker and Orange color
    plot(T_pl(valid_pl), Eg_pl(valid_pl), '-p', 'Color', [1 0.5 0], 'MarkerFaceColor', [1 0.5 0], ...
        'MarkerSize', 14, 'DisplayName', 'PL Lab Measurements', 'LineWidth', 2.5);
end

%% --- 5. PLOT JAIN (1992) THEORY ---
set(0, 'CurrentFigure', h_fig_eg);
T_range = linspace(-65, 25, 200) + 273.15;

% Varshni Parameters for Intrinsic GaAs
Eg_0 = 1.519; alpha_v = 5.405e-4; beta_v = 204;
Eg_intrinsic = Eg_0 - (alpha_v .* T_range.^2) ./ (T_range + beta_v);

% Jain Corrections (Burstein-Moss & BGN)
m_e = 0.067; m_hh = 0.51; 
m_vc_star = (m_e * m_hh) / (m_e + m_hh) * 9.11e-31; 
Delta_BM = (( (1.054e-34)^2 / (2 * m_vc_star) ) * (3*pi^2 * loaded_n_doping*1e6)^(2/3)) / 1.602e-19; 
N_norm = loaded_n_doping / 1e18;
Delta_BGN = (62 * (N_norm)^(1/3) + 7.4 * (N_norm)^(1/4)) / 1000; 
Eg_theory = Eg_intrinsic + Delta_BM - Delta_BGN;

% Plot theory lines behind the scatter data
plot(T_range-273.15, Eg_intrinsic, '--k', 'LineWidth', 2, 'DisplayName', 'Varshni (Intrinsic)');
plot(T_range-273.15, Eg_theory, '-k', 'LineWidth', 3, 'DisplayName', sprintf('Jain 1992 Theory (N_D=%.1e)', loaded_n_doping));

% Text box for theory parameters
str_info = {
    sprintf('\\Delta E_{BM} = +%.3f eV', Delta_BM),
    sprintf('\\Delta E_{BGN} = -%.3f eV', Delta_BGN)
};
text(-60, min(Eg_intrinsic), str_info, 'FontSize', 16, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'VerticalAlignment', 'bottom');

%% --- 6. FINALIZE FIGURE ---
set(0, 'CurrentFigure', h_fig_eg);
xlim([-65 25]); 

% Dynamically adjust Y-limits based on the theory and data range
min_y = min([Eg_intrinsic, min(Eg_pl(valid_pl))]) - 0.02;
max_y = max([Eg_theory, max(Eg_pl(valid_pl))]) + 0.05;
if isempty(max_y), max_y = 1.5; end % Fallback
ylim([min_y max_y]);

% Ensure the legend is properly placed and large enough for a poster
h_leg = legend('show');
set(h_leg, 'Location', 'northwest', 'FontSize', 16);