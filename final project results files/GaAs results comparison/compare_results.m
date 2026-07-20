% =========================================================================
% SCRIPT: Compare Multiple Measurement Conditions (Bandgap Analysis)
% PURPOSE: Compares Bandgap (Eg) extracted via Tauc plot against:
%          1. Independent Photoluminescence (PL) lab measurements.
%          2. Jain (1992) Theory for heavily doped n-GaAs.
% NOTE: Analysis is strictly limited to T >= -60°C due to defect absorption.
% =========================================================================
clear all; close all; clc;

%% === STEP 1: DEFINE DATA SOURCES & VISUAL STYLES ===
% We define the files saved from the Master Script runs.
data_files = {'results_vacuum1.mat', 'results_nitro.mat', 'non_results.mat', 'results_nitro_vacuum.mat'};

% Presentation labels corresponding to the files
run_labels = {'Vacuum Only', 'Nitrogen Flow', 'Cooling Only', 'Vacuum + Nitro'};

% Colors are strictly matched to the Noise Analysis script for visual consistency:
% Green (Vacuum), Red (Nitrogen), Blue (Cooling), Purple (Combined)
colors  = {[0 0.7 0], [1 0 0], [0 0.4 1], [0.6 0 0.8]}; 
markers = {'o', '^', 'v', 'd'};
num_runs = length(data_files);

%% === STEP 2: INITIALIZE COMPARISON FIGURE ===
% Create a large, high-resolution figure optimized for poster/presentation viewing.
h_fig_eg = figure('Color', 'w', 'Name', 'Figure 7: Bandgap Comparison', 'units', 'normalized', 'outerposition', [0.2 0.1 0.6 0.7]);
hold on; grid on; box on;
set(gca, 'FontSize', 16); 

xlabel('Temperature (°C)', 'FontSize', 20, 'FontWeight', 'bold'); 
ylabel('Band Gap Energy (eV)', 'FontSize', 20, 'FontWeight', 'bold'); 
title('GaAs Bandgap: Theory vs. Ellipsometry & PL (Up to -60°C)', 'FontSize', 24, 'FontWeight', 'bold');

%% === STEP 3: LOAD ELLIPSOMETRY DATA AND PLOT ===
% LOOP PURPOSE: 
% Iterates through the saved `.mat` files. It extracts the temperature (T) 
% and Bandgap (Eg) arrays. It filters out temperatures below -60°C because 
% sub-bandgap trap/defect absorption distorts the Tauc plot linearity.

loaded_n_doping = 3.8e18; % Default fallback doping value [cm^-3]

for i = 1:num_runs
    if isfile(data_files{i})
        data = load(data_files{i});
        T    = data.temps_C;
        Eg   = data.extracted_Eg;
        
        % Update doping concentration if it was saved in the .mat file
        if isfield(data, 'n_doping')
            loaded_n_doping = data.n_doping; 
        end
        
        % Logical mask: Only keep valid numbers within the reliable physical range
        valid_eg = ~isnan(Eg) & (T >= -60) & (T <= 25);
        
        % Plot the valid data points
        set(0, 'CurrentFigure', h_fig_eg);
        if any(valid_eg)
            plot(T(valid_eg), Eg(valid_eg), ['-' markers{i}], 'Color', colors{i}, ...
                'MarkerFaceColor', colors{i}, 'MarkerSize', 10, ...
                'DisplayName', run_labels{i}, 'LineWidth', 2);
        end
    else
        warning('File not found in current directory: %s', data_files{i});
    end
end

%% === STEP 4: INTEGRATE INDEPENDENT PL MEASUREMENTS ===
% Photoluminescence (PL) serves as an independent physical validation.
% We input the hardcoded peak emission wavelengths obtained from the lab.

T_pl    = [20, 10, 0, -10, -20, -30, -40, -50, -60, -70, -80, -90, -100];
meas_Pk = [851.32, 849.05, 846.94, 844.86, 842.93, 840.99, 839.20, 837.55, 835.92, 834.4, 832.8, 831.1, 829.4];

% Convert PL peak wavelength (nm) to Energy (eV) using the Planck-Einstein relation:
% E = hc / lambda. In practical units: E(eV) = 1240 / lambda(nm).
Eg_pl = 1240 ./ meas_Pk;

% Ensure they are column vectors for consistent plotting
T_pl  = T_pl(:);
Eg_pl = Eg_pl(:);

% Filter to match the Ellipsometry plot range (down to -60°C) for a fair comparison
valid_pl = ~isnan(Eg_pl) & (T_pl >= -60) & (T_pl <= 25);

set(0, 'CurrentFigure', h_fig_eg);
if any(valid_pl)
    % Plot PL with a distinct Pentagram (Star) marker and Orange color
    plot(T_pl(valid_pl), Eg_pl(valid_pl), '-p', 'Color', [1 0.5 0], ...
        'MarkerFaceColor', [1 0.5 0], 'MarkerSize', 14, ...
        'DisplayName', 'PL Lab Measurements', 'LineWidth', 2.5);
end

%% === STEP 5: COMPUTE & PLOT JAIN (1992) THEORY ===
% We calculate the theoretical expected bandgap for heavily doped n-GaAs.

set(0, 'CurrentFigure', h_fig_eg);
T_range = linspace(-65, 25, 200) + 273.15;

% 1. Varshni Model for Intrinsic GaAs baseline
Eg_0 = 1.519; alpha_v = 5.405e-4; beta_v = 204;
Eg_intrinsic = Eg_0 - (alpha_v .* T_range.^2) ./ (T_range + beta_v);

% 2. Jain Model Quantum Corrections (Burstein-Moss & Bandgap Narrowing)
m_e = 0.067; m_hh = 0.51; 
m_vc_star = (m_e * m_hh) / (m_e + m_hh) * 9.11e-31; 

Delta_BM  = (( (1.054e-34)^2 / (2 * m_vc_star) ) * (3*pi^2 * loaded_n_doping*1e6)^(2/3)) / 1.602e-19; 
N_norm    = loaded_n_doping / 1e18;
Delta_BGN = (62 * (N_norm)^(1/3) + 7.4 * (N_norm)^(1/4)) / 1000; 

% Net Theoretical Bandgap
Eg_theory = Eg_intrinsic + Delta_BM - Delta_BGN;

% Plot theory lines behind the scatter data for clear visibility
plot(T_range-273.15, Eg_intrinsic, '--k', 'LineWidth', 2, 'DisplayName', 'Varshni (Intrinsic)');
plot(T_range-273.15, Eg_theory, '-k', 'LineWidth', 3, 'DisplayName', sprintf('Jain 1992 Theory (N_D=%.1e)', loaded_n_doping));

% Text box explicitly declaring the calculated quantum shifts
str_info = {
    sprintf('\\Delta E_{BM} = +%.3f eV', Delta_BM),
    sprintf('\\Delta E_{BGN} = -%.3f eV', Delta_BGN)
};
text(-60, min(Eg_intrinsic), str_info, 'FontSize', 16, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'VerticalAlignment', 'bottom');

%% === STEP 6: FINALIZE FIGURE LIMITS & LEGEND ===
set(0, 'CurrentFigure', h_fig_eg);
xlim([-65 25]); 

% Dynamically adjust Y-limits based on the theory and data range to prevent cutoff
min_y = min([Eg_intrinsic, min(Eg_pl(valid_pl))]) - 0.02;
max_y = max([Eg_theory, max(Eg_pl(valid_pl))]) + 0.05;

if isempty(max_y)
    max_y = 1.5; % Fallback limit if data arrays are empty
end 
ylim([min_y max_y]);

% Ensure the legend is properly placed and large enough for a poster
h_leg = legend('show');
set(h_leg, 'Location', 'northwest', 'FontSize', 14);

fprintf('\nComparison graph generated successfully.\n');