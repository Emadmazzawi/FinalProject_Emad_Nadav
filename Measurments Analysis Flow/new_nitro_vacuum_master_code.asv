% =========================================================================
% MASTER SCRIPT: GaAs Optical Characterization (Nitrogen-Vacuum)
% Range: 20°C to -110°C  
% MODELS INCLUDED: Varshni & Jain (1992) with Burstein-Moss & BGN
% =========================================================================

clear all; close all; clc;

%% --- USER CONTROLS ---
% Set to true to pop up separate windows for each temperature's Tauc plot
show_individual_tauc = false; 

%% --- GLOBAL PARAMETERS ---
h_const = 6.626e-34;    
c_const = 3e8;          
kb = 8.6173e-5;         
d_thick = 350e-4;       
n_doping = 3.8e18;      % Doping concentration [cm^-3]

% Set temperature range (20C down to -110C in 5-degree steps)
temps_C = 20:-5:-90; 
temps_K = temps_C + 273.15; 
N = length(temps_C);

legend_labels = arrayfun(@(t) sprintf('%d °C', t), temps_C, 'UniformOutput', false);
cmap = flipud(jet(N)); 
extracted_Eg = nan(1, N);

%% --- 0. SMART FILE DISCOVERY ---
opt_files = cell(1, N);
ellip_files = cell(1, N);

for i = 1:N
    T = temps_C(i);
    % Format string matching: positive temps end in 'c', negative in 'm[val]c'
    if T >= 0
        t_str = sprintf('_%dc', T); 
    else
        t_str = sprintf('_m%dc', abs(T)); 
    end
    
    found = dir(['*' t_str '*.xlsx']);
    for k = 1:length(found)
        fname = found(k).name;
        if startsWith(fname, '~$'), continue; end % Skip open Excel temp files
        if contains(fname, '_data')
            ellip_files{i} = fname;
        else
            opt_files{i} = fname; 
        end
    end
end

%% --- 1. INITIALIZE MAIN FIGURES ---
h_fig_ellip     = figure('Color', 'w', 'Name', 'Figure 1: Psi & Delta', 'units', 'normalized', 'outerposition', [0 0 1 1]);
h_fig_opt       = figure('Color', 'w', 'Name', 'Figure 2: Absorption', 'units', 'normalized', 'outerposition', [0 0 1 1]);
h_fig_nk        = figure('Color', 'w', 'Name', 'Figure 3: n and k Index', 'units', 'normalized', 'outerposition', [0 0 1 1]);
h_fig_tauc      = figure('Color', 'w', 'Name', 'Figure 4: Combined Tauc Plots', 'units', 'normalized', 'outerposition', [0 0 1 1]);
h_fig_alpha_log = figure('Color', 'w', 'Name', 'Figure 5: Alpha (Log) Poster', 'units', 'normalized', 'outerposition', [0 0 1 1]);

% Tracking arrays to only label successful plots in legends
ellip_plotted = false(1, N);
opt_plotted   = false(1, N);
tauc_plotted  = false(1, N);
alpha_plotted = false(1, N); 

%% --- 2. MAIN PROCESSING LOOP ---
for i = 1:N
    has_opt = ~isempty(opt_files{i});
    has_ellip = ~isempty(ellip_files{i});
    
    % Process Ellipsometry Data (Psi & Delta)
    if has_ellip
        [L_ellip, psi, delta] = extract_ellip(ellip_files{i});
        ellip_plotted(i) = true;
        
        set(0, 'CurrentFigure', h_fig_ellip);
        subplot(2,1,1); hold on; plot(L_ellip, psi, 'LineWidth', 1.5, 'Color', cmap(i,:));
        subplot(2,1,2); hold on; plot(L_ellip, delta, 'LineWidth', 1.5, 'Color', cmap(i,:));
    end
    
    % Process Optical Properties (n, k, alpha, Tauc)
    if has_opt
        [L_opt, n_idx, k_ext, alpha, E_ev] = extract_opt(opt_files{i});
        opt_plotted(i) = true;
        alpha_plotted(i) = true; 
        
        % Plot Absorption vs Wavelength & Energy
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
        
        % --- Tauc Analysis for Direct Bandgap ---
        y_tauc = (alpha .* E_ev).^2;
        % Compute derivative to find the linear absorption edge
        deriv = diff(smoothdata(y_tauc, 'gaussian', 5)) ./ diff(E_ev);
        E_mid = (E_ev(1:end-1) + E_ev(2:end)) / 2;
        
        % Restrict search window to appropriate physical range for GaAs
        search_range = (E_mid > 1.35 & E_mid < 1.70); 
        
        if any(search_range)
            tauc_plotted(i) = true;
            [~, max_idx_rel] = max(deriv(search_range));
            search_indices = find(search_range);
            center_idx = search_indices(max_idx_rel);
            
            % Define fitting window around the steepest point
            fit_win = (center_idx - 8) : (center_idx + 8);
            fit_win = fit_win(fit_win > 0 & fit_win <= length(E_ev));
            
            % Linear fit (y = mx + b)
            p = polyfit(E_ev(fit_win), y_tauc(fit_win), 1);
            Eg = -p(2) / p(1); % Intercept on x-axis
            extracted_Eg(i) = Eg;
            
            % Generate extrapolation line for plotting
            x_ext = linspace(Eg, E_ev(center_idx) + 0.05, 50);
            y_ext = polyval(p, x_ext);
            
            % Add to Combined Tauc Figure
            set(0, 'CurrentFigure', h_fig_tauc); hold on;
            plot(E_ev, y_tauc, 'LineWidth', 2, 'Color', cmap(i,:));
            plot(x_ext, y_ext, '--', 'Color', cmap(i,:), 'LineWidth', 1.2, 'HandleVisibility', 'off');
            plot(Eg, 0, 'x', 'Color', cmap(i,:), 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
            
            % Pop up individual Tauc plots if requested by user
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

%% --- 3. FINALIZE MAIN FIGURES ---
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
grid on; box on;
yscale("log"); 
xlabel('Energy (eV)', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('\alpha (cm^{-1})', 'FontSize', 24, 'FontWeight', 'bold');
title('Absorption Coefficient (Log Scale) for n-GaAs: 20°C to -100°C', 'FontSize', 28, 'FontWeight', 'bold');
xlim([1.3 1.8]); 
h_leg_alpha = legend(legend_labels(alpha_plotted), 'NumColumns', 1, 'FontSize', 12);
set(h_leg_alpha, 'Location', 'eastoutside');
set(gca, 'FontSize', 16);

%% --- 4. BANDGAP THEORY COMPARISON (VARSHNI & JAIN 1992) ---
h_fig_eg_model = figure('Color', 'w', 'Name', 'Figure 6: Bandgap Comparison', 'units', 'normalized', 'outerposition', [0 0 1 1]);
hold on; 

T_range = linspace(25, -115, 200) + 273.15;

% Varshni Parameters for Intrinsic GaAs
Eg_0 = 1.519; alpha_v = 5.405e-4; beta_v = 204;
Eg_intrinsic = Eg_0 - (alpha_v .* T_range.^2) ./ (T_range + beta_v);

% Jain (1992) Model for Heavily Doped n-GaAs
m_e = 0.067; m_hh = 0.51; 
m_vc_star = (m_e * m_hh) / (m_e + m_hh) * 9.11e-31; 

% Burstein-Moss Shift Calculation (Delta E_BM)
Delta_BM = (( (1.054e-34)^2 / (2 * m_vc_star) ) * (3*pi^2 * n_doping*1e6)^(2/3)) / 1.602e-19; 

% Bandgap Narrowing Calculation (Delta E_BGN)
N_norm = n_doping / 1e18;
Delta_BGN = (62 * (N_norm)^(1/3) + 7.4 * (N_norm)^(1/4)) / 1000; 

% Net Bandgap Model Calculation
Eg_theory = Eg_intrinsic + Delta_BM - Delta_BGN;

plot(T_range-273.15, Eg_intrinsic, '--k', 'LineWidth', 1.5, 'DisplayName', 'Varshni (Intrinsic)');
plot(T_range-273.15, Eg_theory, '-k', 'LineWidth', 2, 'DisplayName', sprintf('Jain 1992 (Doped, N_D=%.1e)', n_doping));
plot(temps_C(~isnan(extracted_Eg)), extracted_Eg(~isnan(extracted_Eg)), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Experiment (Tauc)');

grid on; box on; 
xlabel('Temperature (°C)'); ylabel('Band Gap Energy (eV)'); title('GaAs Bandgap Model vs. Experiment');
legend('Location', 'northeast', 'FontSize', 12); 
xlim([-110 25]); ylim([min(Eg_intrinsic)-0.02 max(Eg_theory)+0.05]);

% Add Theoretical Annotation Box
str_info = {
    sprintf('Jain 1992 Model Details:'),
    sprintf('N_D = %.1e cm^{-3}', n_doping),
    sprintf('\\Delta E_{BM} (Burstein-Moss) = +%.3f eV', Delta_BM),
    sprintf('\\Delta E_{BGN} (BGN) = -%.3f eV', Delta_BGN),
    sprintf('Total Shift = %+.3f eV', Delta_BM - Delta_BGN)
};
text(-105, min(Eg_intrinsic), str_info, 'FontSize', 11, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'VerticalAlignment', 'bottom');

%% --- 5. PRINT RESULTS TABLE ---
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

%% --- 6. EXPORT RESULTS FOR COMPARISON ---
save_filename = 'results_nitro_vacuum.mat'; 
save(save_filename, 'temps_C', 'extracted_Eg', 'n_doping');
fprintf('Results successfully saved to %s\n', save_filename);

%% --- 7. AUTO-SAVE ALL FIGURES ---
fprintf('\nSaving all figures... ');

% Create a dedicated folder for saving the figures
save_folder = 'Saved_Figures';
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

% Find all open figure objects
all_figs = findobj('Type', 'figure');
for k = 1:length(all_figs)
    current_fig = all_figs(k);
    
    % Retrieve the figure name to use as a filename
    fig_title = current_fig.Name;
    if isempty(fig_title)
        fig_title = sprintf('Figure_%d', current_fig.Number); 
    end
    
    % Sanitize the filename by removing invalid characters
    safe_filename = regexprep(fig_title, '[\\/:*?"<>|]', '_'); 
    safe_filename = strrep(safe_filename, ' ', '_');
    
    % Save as a native MATLAB figure file (.fig)
    fig_path = fullfile(save_folder, [safe_filename '.fig']);
    savefig(current_fig, fig_path);
end
fprintf('Done! %d figures saved to folder "%s".\n', length(all_figs), save_folder);

%% =========================================================================
%                  INTERNAL DATA EXTRACTION FUNCTIONS
% =========================================================================
function [L, n, k, a, E] = extract_opt(filename)
    % Function to read Excel files with 'Wavelength', 'n', 'k', 'Alpha', 'Energy'
    data = readtable(filename, 'VariableNamingRule', 'preserve'); 
    raw = data{:, :};
    
    % Convert cell array to double array if mixed types are present
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
    
    % Assign column variables (Standard J.A. Woollam .xlsx format)
    L = raw(:, 1); 
    n = raw(:, 2); 
    k = raw(:, 3); 
    a = raw(:, 4); 
    E = raw(:, 5);
end

function [L, psi, delta] = extract_ellip(filename)
    % Function to read Excel files with 'Wavelength', 'AoI', 'Psi', 'Delta'
    data = readtable(filename, 'VariableNamingRule', 'preserve'); 
    raw = data{:, :};
    
    % Convert cell array to double array if mixed types are present
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
    
    % Assign column variables
    L = raw(:, 1); 
    psi = raw(:, 3); 
    delta = raw(:, 4); 
end