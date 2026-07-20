% =========================================================================
% MASTER SCRIPT: GaAs Optical Characterization (Nitrogen-Vacuum)
% Range: 20°C to -110°C  
% MODELS INCLUDED: Varshni & Jain (1992) with Burstein-Moss & BGN
% =========================================================================
clear all; close all; clc;

%% === STEP 1: USER-DEFINED VARIABLES & SETUP ===
% NOTE: Ensure all raw data Excel files (.xlsx) are located in the same 
% directory as this MATLAB script before running.

show_individual_tauc = false;  % Toggle (true/false): Pop up separate Tauc plots
d_thick              = 350e-4; % Sample thickness [cm]
n_doping             = 3.8e18; % Doping concentration [cm^-3]

% Temperature range: 20°C down to -90°C in 5-degree increments
% WARNING: Energy gap extraction at T < -60°C may be distorted due to 
% dominant absorption from sub-bandgap defects (traps).
temps_C = 20:-5:-90; 
temps_K = temps_C + 273.15; 
N       = length(temps_C);

% Labels and colors for plotting
legend_labels = arrayfun(@(t) sprintf('%d °C', t), temps_C, 'UniformOutput', false);
cmap          = flipud(jet(N)); 
extracted_Eg  = nan(1, N); % Array to store calculated Bandgap energies

%% === STEP 2: PHYSICAL CONSTANTS ===
% Fundamental constants that do not change between experimental runs.

h_const = 6.626e-34; % Planck's constant [J*s]
c_const = 3e8;       % Speed of light [m/s]
kb      = 8.6173e-5; % Boltzmann constant [eV/K]

%% === STEP 3: FILE DISCOVERY ===
% LOOP PURPOSE: 
% This loop iterates through the defined temperature array. For each temperature, 
% it searches the current directory for matching Excel files based on a specific 
% naming convention (e.g., '_20c' for positive, '_m20c' for negative). 
% It sorts them into optical data or ellipsometry data.

opt_files   = cell(1, N);
ellip_files = cell(1, N);

for i = 1:N
    T = temps_C(i);
    
    % Format string matching
    if T >= 0
        t_str = sprintf('_%dc', T); 
    else
        t_str = sprintf('_m%dc', abs(T)); 
    end
    
    found = dir(['*' t_str '*.xlsx']);
    for k = 1:length(found)
        fname = found(k).name;
        if startsWith(fname, '~$'), continue; end % Skip open temp files
        
        if contains(fname, '_data')
            ellip_files{i} = fname;
        else
            opt_files{i} = fname; 
        end
    end
end

%% === STEP 4: INITIALIZE MAIN FIGURES ===
h_fig_ellip     = figure('Color', 'w', 'Name', 'Figure 1: Psi & Delta', 'units', 'normalized', 'outerposition', [0 0 1 1]);
h_fig_opt       = figure('Color', 'w', 'Name', 'Figure 2: Absorption', 'units', 'normalized', 'outerposition', [0 0 1 1]);
h_fig_nk        = figure('Color', 'w', 'Name', 'Figure 3: n and k Index', 'units', 'normalized', 'outerposition', [0 0 1 1]);
h_fig_tauc      = figure('Color', 'w', 'Name', 'Figure 4: Combined Tauc Plots', 'units', 'normalized', 'outerposition', [0 0 1 1]);
h_fig_alpha_log = figure('Color', 'w', 'Name', 'Figure 5: Alpha (Log) Poster', 'units', 'normalized', 'outerposition', [0 0 1 1]);

% Tracking arrays to ensure only successfully processed data appears in legends
ellip_plotted = false(1, N);
opt_plotted   = false(1, N);
tauc_plotted  = false(1, N);
alpha_plotted = false(1, N); 

%% === STEP 5: MAIN PROCESSING & CALCULATION ===
% LOOP PURPOSE:
% This is the core analysis loop. It extracts the raw data using helper 
% functions, plots the optical properties, and performs the Tauc method 
% linear fitting to extract the experimental bandgap (Eg).

for i = 1:N
    has_opt   = ~isempty(opt_files{i});
    has_ellip = ~isempty(ellip_files{i});
    
    % --- 5A. Process Ellipsometry Data (Psi & Delta) ---
    if has_ellip
        [L_ellip, psi, delta] = extract_ellip(ellip_files{i});
        ellip_plotted(i) = true;
        
        set(0, 'CurrentFigure', h_fig_ellip);
        subplot(2,1,1); hold on; plot(L_ellip, psi, 'LineWidth', 1.5, 'Color', cmap(i,:));
        subplot(2,1,2); hold on; plot(L_ellip, delta, 'LineWidth', 1.5, 'Color', cmap(i,:));
    end
    
    % --- 5B. Process Optical Properties & Tauc Analysis ---
    if has_opt
        [L_opt, n_idx, k_ext, alpha, E_ev] = extract_opt(opt_files{i});
        opt_plotted(i)   = true;
        alpha_plotted(i) = true; 
        
        % Plot Absorption
        set(0, 'CurrentFigure', h_fig_opt);
        subplot(1,2,1); hold on; plot(L_opt, alpha, 'LineWidth', 1.5, 'Color', cmap(i,:));
        subplot(1,2,2); hold on; plot(E_ev, max(alpha, 1), 'LineWidth', 1.5, 'Color', cmap(i,:)); 
        
        % Plot Refractive Index (n) & Extinction Coefficient (k)
        set(0, 'CurrentFigure', h_fig_nk);
        subplot(2,1,1); hold on; plot(L_opt, n_idx, 'LineWidth', 1.5, 'Color', cmap(i,:));
        subplot(2,1,2); hold on; plot(L_opt, k_ext, 'LineWidth', 1.5, 'Color', cmap(i,:));
        
        % Plot Alpha Log for Poster Display
        set(0, 'CurrentFigure', h_fig_alpha_log); hold on;
        plot(E_ev, max(alpha, 1), 'LineWidth', 2.5, 'Color', cmap(i,:)); 
        
        % --- Tauc Analysis Calculations ---
        y_tauc = (alpha .* E_ev).^2;
        
        % Calculate derivative to find the steepest slope (absorption edge).
        % We use 'smoothdata' with a 'gaussian' window to filter out high-frequency 
        % numerical noise before derivation. The gaussian shape is ideal because it 
        % smooths the curve smoothly without shifting the position of the peak.
        deriv  = diff(smoothdata(y_tauc, 'gaussian', 5)) ./ diff(E_ev);
        E_mid  = (E_ev(1:end-1) + E_ev(2:end)) / 2;
        
        % Restrict fit to the physical GaAs range
        search_range = (E_mid > 1.35 & E_mid < 1.70); 
        
        if any(search_range)
            tauc_plotted(i) = true;
            [~, max_idx_rel] = max(deriv(search_range));
            search_indices = find(search_range);
            center_idx = search_indices(max_idx_rel);
            
            % 17-point fitting window around the steepest slope
            fit_win = (center_idx - 8) : (center_idx + 8);
            fit_win = fit_win(fit_win > 0 & fit_win <= length(E_ev));
            
            % Linear fit (y = mx + b)
            p = polyfit(E_ev(fit_win), y_tauc(fit_win), 1);
            Eg = -p(2) / p(1); % Find x-intercept
            extracted_Eg(i) = Eg;
            
            % Extrapolation line for visualization
            x_ext = linspace(Eg, E_ev(center_idx) + 0.05, 50);
            y_ext = polyval(p, x_ext);
            
            % Plot Combined Tauc Data
            set(0, 'CurrentFigure', h_fig_tauc); hold on;
            plot(E_ev, y_tauc, 'LineWidth', 2, 'Color', cmap(i,:));
            plot(x_ext, y_ext, '--', 'Color', cmap(i,:), 'LineWidth', 1.2, 'HandleVisibility', 'off');
            plot(Eg, 0, 'x', 'Color', cmap(i,:), 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
            
            % Individual Tauc Plots (if enabled)
            if show_individual_tauc
                fig_ind = figure('Color', 'w', 'Name', sprintf('Tauc Plot - %s', legend_labels{i}));
                hold on;
                plot(E_ev, y_tauc, 'b', 'LineWidth', 2);
                plot(x_ext, y_ext, 'r--', 'LineWidth', 2);
                plot(Eg, 0, 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
                grid on; box on; xlim([1.35 1.70]);
                xlabel('Energy (eV)'); ylabel('(\alpha E)^2');
                title(sprintf('Tauc Plot: %s | Extracted E_g = %.4f eV', legend_labels{i}, Eg));
            end
        end
    end
end

%% === STEP 6: FINALIZE & FORMAT MAIN FIGURES ===
% Figure 1: Psi & Delta
set(0, 'CurrentFigure', h_fig_ellip);
subplot(2,1,1); grid on; box on; xlabel('Wavelength (nm)'); ylabel('\Psi (deg)'); title('Psi (\Psi)'); xlim([250 1000]);
subplot(2,1,2); grid on; box on; xlabel('Wavelength (nm)'); ylabel('\Delta (deg)'); title('Delta (\Delta)'); xlim([250 1000]);
h_leg1 = legend(legend_labels(ellip_plotted), 'NumColumns', 1);
set(h_leg1, 'Position', [0.90 0.3 0.08 0.4], 'Units', 'normalized');

% Figure 2: Absorption
set(0, 'CurrentFigure', h_fig_opt);
subplot(1,2,1); grid on; box on; xlabel('Wavelength (nm)'); ylabel('\alpha (cm^{-1})'); title('Absorption vs Wavelength'); xlim([400 1000]);
subplot(1,2,2); grid on; box on; yscale("log"); xlabel('Energy (eV)'); ylabel('\alpha (cm^{-1})'); title('Absorption (Log)'); xlim([1.3 2.5]);
h_leg2 = legend(legend_labels(opt_plotted), 'NumColumns', 1);
set(h_leg2, 'Position', [0.90 0.3 0.08 0.4], 'Units', 'normalized');

% Figure 3: n and k Indices
set(0, 'CurrentFigure', h_fig_nk);
subplot(2,1,1); grid on; box on; xlabel('Wavelength (nm)'); ylabel('n'); title('n Index'); xlim([400 1000]);
subplot(2,1,2); grid on; box on; xlabel('Wavelength (nm)'); ylabel('k'); title('k Index'); xlim([400 1000]);
h_leg_nk = legend(legend_labels(opt_plotted), 'NumColumns', 1);
set(h_leg_nk, 'Position', [0.90 0.3 0.08 0.4], 'Units', 'normalized');

% Figure 4: Combined Tauc Plots
set(0, 'CurrentFigure', h_fig_tauc);
grid on; box on; xlabel('Energy (eV)'); ylabel('(\alpha E)^2'); title('Combined Tauc Plots'); xlim([1.35 1.70]); ylim([0 1.2e9]);
legend(legend_labels(tauc_plotted), 'Location', 'eastoutside');

% Figure 5: Alpha Log (Poster View)
set(0, 'CurrentFigure', h_fig_alpha_log);
grid on; box on; yscale("log"); 
xlabel('Energy (eV)', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('\alpha (cm^{-1})', 'FontSize', 24, 'FontWeight', 'bold');
title('Absorption Coefficient (Log Scale) for n-GaAs: 20°C to -100°C', 'FontSize', 28, 'FontWeight', 'bold');
xlim([1.3 1.8]); 
h_leg_alpha = legend(legend_labels(alpha_plotted), 'NumColumns', 1, 'FontSize', 12);
set(h_leg_alpha, 'Location', 'eastoutside'); set(gca, 'FontSize', 16);

%% === STEP 7: BANDGAP THEORY COMPARISON (VARSHNI & JAIN 1992) ===
h_fig_eg_model = figure('Color', 'w', 'Name', 'Figure 6: Bandgap Comparison', 'units', 'normalized', 'outerposition', [0 0 1 1]);
hold on; 
T_range = linspace(25, -115, 200) + 273.15;

% --- Varshni Model for Intrinsic GaAs ---
% Formula: Eg(T) = Eg(0) - (alpha * T^2) / (T + beta)
% Where:
%   Eg(0) = Band gap at 0 Kelvin
%   alpha = Empirical constant A (temperature coefficient) specific to the semiconductor
%   beta  = Empirical constant B (closely related to the Debye temperature)
%
% Constants specific to GaAs:
Eg_0    = 1.519;    % Bandgap at 0K [eV]
alpha_v = 5.405e-4; % Constant A [eV/K]
beta_v  = 204;      % Constant B [K]

Eg_intrinsic = Eg_0 - (alpha_v .* T_range.^2) ./ (T_range + beta_v);

% --- Jain (1992) Model for Heavily Doped n-GaAs ---
% This model calculates the net effective bandgap shift in heavily doped 
% semiconductors. It accounts for two competing quantum effects:
% 1. Burstein-Moss Shift (Delta_BM): Widens the bandgap as free electrons 
%    fill the bottom of the conduction band.
% 2. Bandgap Narrowing (Delta_BGN): Shrinks the bandgap due to many-body 
%    interactions between charge carriers and dopant ions.

m_e = 0.067; m_hh = 0.51; 
m_vc_star = (m_e * m_hh) / (m_e + m_hh) * 9.11e-31; 

% Physical Shifts Calculations
Delta_BM  = (( (1.054e-34)^2 / (2 * m_vc_star) ) * (3*pi^2 * n_doping*1e6)^(2/3)) / 1.602e-19; 
N_norm    = n_doping / 1e18;
Delta_BGN = (62 * (N_norm)^(1/3) + 7.4 * (N_norm)^(1/4)) / 1000; 

% Net Bandgap Model Calculation
Eg_theory = Eg_intrinsic + Delta_BM - Delta_BGN;

% Plotting Theoretical vs Experimental
plot(T_range-273.15, Eg_intrinsic, '--k', 'LineWidth', 1.5, 'DisplayName', 'Varshni (Intrinsic)');
plot(T_range-273.15, Eg_theory, '-k', 'LineWidth', 2, 'DisplayName', sprintf('Jain 1992 (Doped, N_D=%.1e)', n_doping));
plot(temps_C(~isnan(extracted_Eg)), extracted_Eg(~isnan(extracted_Eg)), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Experiment (Tauc)');

grid on; box on; 
xlabel('Temperature (°C)'); ylabel('Band Gap Energy (eV)'); title('GaAs Bandgap Model vs. Experiment');
legend('Location', 'northeast', 'FontSize', 12); 
xlim([-110 25]); ylim([min(Eg_intrinsic)-0.02 max(Eg_theory)+0.05]);

% Theoretical Annotation Box
str_info = {
    sprintf('Jain 1992 Model Details:'),
    sprintf('N_D = %.1e cm^{-3}', n_doping),
    sprintf('\\Delta E_{BM} (Burstein-Moss) = +%.3f eV', Delta_BM),
    sprintf('\\Delta E_{BGN} (BGN) = -%.3f eV', Delta_BGN),
    sprintf('Total Shift = %+.3f eV', Delta_BM - Delta_BGN)
};
text(-105, min(Eg_intrinsic), str_info, 'FontSize', 11, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'VerticalAlignment', 'bottom');

%% === STEP 8: CONSOLE OUTPUT & DATA EXPORT ===
fprintf('\n===================================================\n        FINAL BAND GAP (Eg) RESULTS (TAUC)         \n===================================================\n');
fprintf('%-25s | %-15s\n', 'Temperature Condition', 'Band Gap (eV)');
fprintf('---------------------------------------------------\n');
for i = 1:N
    if ~isnan(extracted_Eg(i))
        fprintf('%-25s | %.4f eV\n', legend_labels{i}, extracted_Eg(i));
    else
        fprintf('%-25s | %s\n', legend_labels{i}, 'Fit Failed'); 
    end
end
fprintf('===================================================\n');

save_filename = 'results_nitro_vacuum.mat'; 
save(save_filename, 'temps_C', 'extracted_Eg', 'n_doping');
fprintf('Results successfully saved to %s\n', save_filename);

%% === STEP 9: AUTO-SAVE ALL FIGURES ===
% LOOP PURPOSE: 
% Iterates through all currently generated MATLAB figure windows. It fetches 
% the figure title, sanitizes the name, and saves them neatly into a local directory.

fprintf('\nSaving all figures... (This may take a minute, please wait)\n');
save_folder = 'Saved_Figures';
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

all_figs = findobj('Type', 'figure');
for k = 1:length(all_figs)
    current_fig = all_figs(k);
    
    fig_title = current_fig.Name;
    if isempty(fig_title)
        fig_title = sprintf('Figure_%d', current_fig.Number); 
    end
    
    % Ensure filename is safe for Windows/Mac
    safe_filename = regexprep(fig_title, '[\\/:*?"<>|]', '_'); 
    safe_filename = strrep(safe_filename, ' ', '_');
    
    % Print progress to the console so we know it's not frozen
    fprintf('  Saving: %s...\n', safe_filename);
    
    % 1. Save as MATLAB .fig file (for future editing)
    fig_path = fullfile(save_folder, [safe_filename '.fig']);
    savefig(current_fig, fig_path);
    
    % 2. Save as high-res PNG (Great for PowerPoint presentations!)
    png_path = fullfile(save_folder, [safe_filename '.png']);
    exportgraphics(current_fig, png_path, 'Resolution', 300);
end
fprintf('Done! %d figures saved to folder "%s".\n\n', length(all_figs), save_folder);

%% =========================================================================
%                  HELPER FUNCTIONS (DATA EXTRACTION)
% =========================================================================

function [L, n, k, a, E] = extract_opt(filename)
    % Extracts optical variables (Wavelength, n, k, Alpha, Energy) 
    % from a standard J.A. Woollam formatted Excel file.
    data = readtable(filename, 'VariableNamingRule', 'preserve'); 
    raw = data{:, :};
    
    if iscell(raw)
        proc = nan(size(raw));
        for r = 1:size(raw,1)
            for c = 1:size(raw,2)
                if isnumeric(raw{r,c})
                    proc(r,c) = raw{r,c};
                else
                    proc(r,c) = str2double(string(raw{r,c}));
                end
            end
        end
        raw = proc;
    end
    
    L = raw(:, 1); 
    n = raw(:, 2); 
    k = raw(:, 3); 
    a = raw(:, 4); 
    E = raw(:, 5);
end

function [L, psi, delta] = extract_ellip(filename)
    % Extracts ellipsometry variables (Wavelength, Psi, Delta) 
    % from an Excel file, converting cell formatting if needed.
    data = readtable(filename, 'VariableNamingRule', 'preserve'); 
    raw = data{:, :};
    
    if iscell(raw)
        proc = nan(size(raw));
        for r = 1:size(raw,1)
            for c = 1:size(raw,2)
                if isnumeric(raw{r,c})
                    proc(r,c) = raw{r,c};
                else
                    proc(r,c) = str2double(string(raw{r,c}));
                end
            end
        end
        raw = proc;
    end
    
    L     = raw(:, 1); 
    psi   = raw(:, 3); 
    delta = raw(:, 4); 
end