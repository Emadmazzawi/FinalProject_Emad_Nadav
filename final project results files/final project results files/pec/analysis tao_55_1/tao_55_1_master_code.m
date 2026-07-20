% =========================================================================
% MASTER SCRIPT: PEC Optical Characterization (Ta3N5)
% Range: 20°C to -110°C | No Theoretical Bandgap Models Included
% Extracted Bandgap (~2.0 - 2.1 eV) vs. Temperature (Indirect Bandgap)
% =========================================================================
clear all; close all; clc;

%% === STEP 1: USER-DEFINED VARIABLES & SETUP ===
% Toggle to true to pop up a separate fitting window for each temperature.
show_individual_tauc = false; 

% Temperature range: 20°C down to -110°C in 10-degree steps.
temps_C = 20:-10:-110; 
N       = length(temps_C);

legend_labels = arrayfun(@(t) sprintf('%d C', t), temps_C, 'UniformOutput', false);
cmap          = flipud(jet(N)); 
extracted_Eg  = nan(1, N); % Array to store calculated Bandgap energies

%% === STEP 2: PHYSICAL CONSTANTS ===
h_const = 6.626e-34; % Planck's constant [J*s]
c_const = 3e8;       % Speed of light [m/s]
kb      = 8.6173e-5; % Boltzmann constant [eV/K]

%% === STEP 3: FILE DISCOVERY ===
% LOOP PURPOSE: 
% Automatically scan the current directory for WVASE Excel files matching 
% the requested temperatures and sort them into Optical and Ellipsometry files.

opt_files   = cell(1, N);
ellip_files = cell(1, N);

for i = 1:N
    T = temps_C(i);
    
    if T >= 0
        t_str = sprintf('_%dc', T); 
    else
        t_str = sprintf('_m%dc', abs(T)); 
    end
    
    found = dir(['*' t_str '*.xlsx']);
    for k = 1:length(found)
        fname = found(k).name;
        if startsWith(fname, '~$'), continue; end % Skip temp files
        
        if contains(fname, '_data')
            ellip_files{i} = fname;
        else
            opt_files{i} = fname; 
        end
    end
end

%% === STEP 4: INITIALIZE MAIN FIGURES ===
% Pre-allocate the figure windows for the dashboard display.
h_fig_ellip     = figure('Color', 'w', 'Name', 'Figure 1: Psi & Delta', 'units', 'normalized', 'outerposition', [0 0 1 1]);
h_fig_opt       = figure('Color', 'w', 'Name', 'Figure 2: Optical Constants', 'units', 'normalized', 'outerposition', [0 0 1 1]);
h_fig_tauc      = figure('Color', 'w', 'Name', 'Figure 3: Tauc Plots', 'units', 'normalized', 'outerposition', [0 0 1 1]);
h_fig_alpha_log = figure('Color', 'w', 'Name', 'Figure 4: Alpha (Log) Poster', 'units', 'normalized', 'outerposition', [0 0 1 1]);

% Tracking arrays for dynamic legends
ellip_plotted = false(1, N);
opt_plotted   = false(1, N);
tauc_plotted  = false(1, N);
alpha_plotted = false(1, N);

%% === STEP 5: MAIN PROCESSING & TAUC ANALYSIS ===
% LOOP PURPOSE:
% Extracts raw data, plots the optical spectra, and performs a Tauc Plot 
% analysis specifically tailored for INDIRECT bandgap materials.

for i = 1:N
    has_opt   = ~isempty(opt_files{i});
    has_ellip = ~isempty(ellip_files{i});
    
    % --- 5A. Process Ellipsometry Data ---
    if has_ellip
        [L_ellip, psi, delta] = extract_ellip(ellip_files{i});
        ellip_plotted(i) = true;
        
        set(0, 'CurrentFigure', h_fig_ellip);
        subplot(1,2,1); hold on; plot(L_ellip, psi, 'LineWidth', 1.5, 'Color', cmap(i,:));
        subplot(1,2,2); hold on; plot(L_ellip, delta, 'LineWidth', 1.5, 'Color', cmap(i,:));
    end
    
    % --- 5B. Process Optical Properties ---
    if has_opt
        [L_opt, n_idx, k_ext, alpha, E_ev] = extract_opt(opt_files{i});
        opt_plotted(i)   = true;
        alpha_plotted(i) = true; 
        
        set(0, 'CurrentFigure', h_fig_opt);
        subplot(2,2,1); hold on; plot(L_opt, alpha, 'LineWidth', 1.5, 'Color', cmap(i,:));
        subplot(2,2,2); hold on; plot(E_ev, max(alpha, 1), 'LineWidth', 1.5, 'Color', cmap(i,:)); 
        subplot(2,2,3); hold on; plot(L_opt, n_idx, 'LineWidth', 1.5, 'Color', cmap(i,:));
        subplot(2,2,4); hold on; plot(L_opt, k_ext, 'LineWidth', 1.5, 'Color', cmap(i,:));
        
        % Alpha (Log) vs. Energy (for poster presentation)
        set(0, 'CurrentFigure', h_fig_alpha_log); hold on;
        plot(E_ev, max(alpha, 1), 'LineWidth', 2.5, 'Color', cmap(i,:)); 
        
        % --- 5C. TAUC ANALYSIS (INDIRECT BANDGAP) ---
        % Ta3N5 is an indirect bandgap material. Therefore, the Tauc relation 
        % requires taking the square root: (alpha * E)^0.5.
        y_tauc = sqrt(alpha .* E_ev);
        
        % Derivative to find the linear absorption edge (using gaussian smoothing)
        deriv = diff(smoothdata(y_tauc, 'gaussian', 5)) ./ diff(E_ev);
        E_mid = (E_ev(1:end-1) + E_ev(2:end)) / 2;
        
        % Search window explicitly set around the expected ~2.1 eV bandgap for Ta3N5
        search_range = (E_mid > 1.8 & E_mid < 2.5); 
        
        if any(search_range)
            tauc_plotted(i) = true;
            [~, max_idx_rel] = max(deriv(search_range));
            search_indices = find(search_range);
            center_idx = search_indices(max_idx_rel);
            
            % 17-point window for linear fitting (y = mx + b)
            fit_win = (center_idx - 8) : (center_idx + 8);
            fit_win = fit_win(fit_win > 0 & fit_win <= length(E_ev));
            
            p  = polyfit(E_ev(fit_win), y_tauc(fit_win), 1);
            Eg = -p(2) / p(1); % Extract x-intercept
            extracted_Eg(i) = Eg;
            
            % Extrapolation line for plotting
            x_ext = linspace(Eg, E_ev(center_idx) + 0.1, 50);
            y_ext = polyval(p, x_ext);
            
            % Add to Combined Tauc Figure
            set(0, 'CurrentFigure', h_fig_tauc); hold on;
            plot(E_ev, y_tauc, 'LineWidth', 2, 'Color', cmap(i,:));
            plot(x_ext, y_ext, '--', 'Color', cmap(i,:), 'LineWidth', 1.2, 'HandleVisibility', 'off');
            plot(Eg, 0, 'x', 'Color', cmap(i,:), 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
            
            % Individual Tauc Plots
            if show_individual_tauc
                fig_ind = figure('Color', 'w', 'Name', sprintf('Tauc Plot - %s', legend_labels{i}));
                hold on;
                plot(E_ev, y_tauc, 'b', 'LineWidth', 2);
                plot(x_ext, y_ext, 'r--', 'LineWidth', 2);
                plot(Eg, 0, 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
                grid on; box on; xlim([1.8 2.5]);
                xlabel('Energy (eV)'); ylabel('(\alpha E)^{1/2}');
                title(sprintf('Tauc Plot (Indirect): %s | Extracted E_g = %.4f eV', legend_labels{i}, Eg));
            end
        end
    end
end

%% === STEP 6: FINALIZE MAIN FIGURES FORMATTING ===
% Apply final cosmetic settings, labels, and limits to all figures.

% Figure 1
set(0, 'CurrentFigure', h_fig_ellip);
subplot(1,2,1); grid on; box on; xlabel('Wavelength (nm)'); ylabel('\Psi (deg)'); title('Psi (\Psi)'); xlim([250 1000]);
subplot(1,2,2); grid on; box on; xlabel('Wavelength (nm)'); ylabel('\Delta (deg)'); title('Delta (\Delta)'); xlim([250 1000]);
if any(ellip_plotted)
    h_leg1 = legend(legend_labels(ellip_plotted), 'NumColumns', 1);
    set(h_leg1, 'Position', [0.90 0.3 0.08 0.4], 'Units', 'normalized');
end

% Figure 2
set(0, 'CurrentFigure', h_fig_opt);
subplot(2,2,1); grid on; box on; xlabel('Wavelength (nm)'); ylabel('\alpha (cm^{-1})'); title('Absorption vs Wavelength'); xlim([400 1000]);
subplot(2,2,2); grid on; box on; yscale("log"); xlabel('Energy (eV)'); ylabel('\alpha (cm^{-1})'); title('Absorption (Log)'); xlim([1.5 3.0]);
subplot(2,2,3); grid on; box on; xlabel('Wavelength (nm)'); ylabel('n'); title('n Index'); xlim([400 1000]);
subplot(2,2,4); grid on; box on; xlabel('Wavelength (nm)'); ylabel('k'); title('k Index'); xlim([400 1000]);
if any(opt_plotted)
    h_leg2 = legend(legend_labels(opt_plotted), 'NumColumns', 1);
    set(h_leg2, 'Position', [0.90 0.3 0.08 0.4], 'Units', 'normalized');
end

% Figure 3
set(0, 'CurrentFigure', h_fig_tauc);
grid on; box on; xlabel('Energy (eV)'); ylabel('(\alpha E)^{1/2}'); title('Combined Tauc Plots (Indirect)'); 
xlim([1.8 2.5]); 
if any(tauc_plotted)
    legend(legend_labels(tauc_plotted), 'Location', 'eastoutside');
end

% Figure 4 (Alpha Log Poster)
set(0, 'CurrentFigure', h_fig_alpha_log);
grid on; box on; yscale("log"); 
xlabel('Energy (eV)', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('\alpha (cm^{-1})', 'FontSize', 24, 'FontWeight', 'bold');
title('Absorption Coefficient (Log Scale): 20°C to -110°C', 'FontSize', 28, 'FontWeight', 'bold');
xlim([1.7 2.8]); 
if any(alpha_plotted)
    h_leg_alpha = legend(legend_labels(alpha_plotted), 'NumColumns', 1, 'FontSize', 12);
    set(h_leg_alpha, 'Location', 'eastoutside');
end
set(gca, 'FontSize', 16);

%% === STEP 7: EXPERIMENTAL BANDGAP VS. TEMPERATURE (ISOLATED) ===
% Plots the pure experimental trend since no theoretical model is evaluated here.

h_fig_eg_trend = figure('Color', 'w', 'Name', 'Figure 5: Extracted Bandgap Trend', 'units', 'normalized', 'outerposition', [0.2 0.2 0.6 0.6]);
hold on; 

valid_idx = ~isnan(extracted_Eg);
valid_temps = temps_C(valid_idx);
valid_Eg = extracted_Eg(valid_idx);

if any(valid_idx)
    plot(valid_temps, valid_Eg, '-ob', 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Experimental E_g (Tauc)');
    
    % Annotate exact values near the points
    for i = 1:length(valid_temps)
        text(valid_temps(i), valid_Eg(i) + 0.002, sprintf('%.3f', valid_Eg(i)), ...
            'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end

grid on; box on; 
xlabel('Temperature (°C)', 'FontSize', 16, 'FontWeight', 'bold'); 
ylabel('Band Gap Energy (eV)', 'FontSize', 16, 'FontWeight', 'bold'); 
title('Extracted PEC Bandgap vs. Temperature', 'FontSize', 20, 'FontWeight', 'bold');
xlim([-120 30]); set(gca, 'FontSize', 14);
legend('Location', 'best', 'FontSize', 14);

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

save_filename = 'results_PEC_bandgap.mat'; 
save(save_filename, 'temps_C', 'extracted_Eg');
fprintf('Results successfully saved to %s\n', save_filename);

%% === STEP 9: AUTO-SAVE ALL FIGURES ===
% Automatically saves all unclosed figures as both editable .fig and high-res .png
fprintf('\nSaving all figures... (This may take a minute, please wait)\n');
save_folder = 'Saved_Figures_PEC';
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

% Use findall to robustly fetch all figure handles
all_figs = findall(0, 'Type', 'figure');
for k = 1:length(all_figs)
    current_fig = all_figs(k);
    
    % Ensure figure wasn't manually closed during loop execution
    if ~isvalid(current_fig)
        continue;
    end
    
    fig_title = current_fig.Name;
    if isempty(fig_title)
        fig_title = sprintf('Figure_%d', current_fig.Number); 
    end
    
    safe_filename = regexprep(fig_title, '[\\/:*?"<>|]', '_'); 
    safe_filename = strrep(safe_filename, ' ', '_');
    
    fprintf('  Saving: %s...\n', safe_filename);
    
    % Save Native MATLAB Figure
    fig_path = fullfile(save_folder, [safe_filename '.fig']);
    savefig(current_fig, fig_path);
    
    % Save High-Res PNG for Presentation
    if isvalid(current_fig)
        png_path = fullfile(save_folder, [safe_filename '.png']);
        exportgraphics(current_fig, png_path, 'Resolution', 300);
    end
end
fprintf('Done! Figures saved to folder "%s".\n\n', save_folder);

%% =========================================================================
%                  HELPER FUNCTIONS (DATA EXTRACTION)
% =========================================================================

function [L, n, k, a, E] = extract_opt(filename)
    data = readtable(filename, 'VariableNamingRule', 'preserve'); raw = data{:, :};
    if iscell(raw)
        proc = nan(size(raw));
        for r=1:size(raw,1)
            for c=1:size(raw,2)
                if isnumeric(raw{r,c})
                    proc(r,c)=raw{r,c};
                else
                    proc(r,c)=str2double(string(raw{r,c}));
                end
            end
        end
        raw = proc;
    end
    L = raw(:, 1); n = raw(:, 2); k = raw(:, 3); a = raw(:, 4); E = raw(:, 5);
end

function [L, psi, delta] = extract_ellip(filename)
    data = readtable(filename, 'VariableNamingRule', 'preserve'); raw = data{:, :};
    if iscell(raw)
        proc = nan(size(raw));
        for r=1:size(raw,1)
            for c=1:size(raw,2)
                if isnumeric(raw{r,c})
                    proc(r,c)=raw{r,c};
                else
                    proc(r,c)=str2double(string(raw{r,c}));
                end
            end
        end
        raw = proc;
    end
    L = raw(:, 1); psi = raw(:, 3); delta = raw(:, 4); 
end