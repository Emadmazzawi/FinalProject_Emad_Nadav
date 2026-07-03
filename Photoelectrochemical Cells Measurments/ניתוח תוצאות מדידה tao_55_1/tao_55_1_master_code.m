% =========================================================================
% MASTER SCRIPT: PEC Optical Characterization
% Range: 20 C to -110 C | English Comments | Table Output
% Extracted Bandgap (~2.1 eV) vs. Temperature (No Theoretical Models)
% =========================================================================
clear all; close all; clc;

%% --- USER CONTROLS ---
show_individual_tauc = true; % Set to true to pop up separate windows for each temp

%% --- GLOBAL PARAMETERS ---
h_const = 6.626e-34;    
c_const = 3e8;          
kb = 8.6173e-5;         

% SET TEMPERATURE RANGE (20 C to -110 C in 10-degree steps)
temps_C = 20:-10:-110; 
N = length(temps_C);
legend_labels = arrayfun(@(t) sprintf('%d C', t), temps_C, 'UniformOutput', false);
cmap = flipud(jet(N)); 

extracted_Eg = nan(1, N);

%% --- 0. SMART FILE DISCOVERY ---
opt_files = cell(1, N);
ellip_files = cell(1, N);

for i = 1:N
    T = temps_C(i);
    if T >= 0, t_str = sprintf('_%dc', T); else, t_str = sprintf('_m%dc', abs(T)); end
    
    found = dir(['*' t_str '*.xlsx']);
    for k = 1:length(found)
        fname = found(k).name;
        if startsWith(fname, '~$'), continue; end 
        if contains(fname, '_data'), ellip_files{i} = fname;
        else, opt_files{i} = fname; end
    end
end

%% --- 1. INITIALIZE MAIN FIGURES ---
h_fig_ellip     = figure('Color', 'w', 'Name', 'Psi & Delta', 'units', 'normalized', 'outerposition', [0 0 1 1]);
h_fig_opt       = figure('Color', 'w', 'Name', 'Optical Constants', 'units', 'normalized', 'outerposition', [0 0 1 1]);
h_fig_tauc      = figure('Color', 'w', 'Name', 'Tauc Plots', 'units', 'normalized', 'outerposition', [0 0 1 1]);
h_fig_alpha_log = figure('Color', 'w', 'Name', 'Alpha (Log) Poster', 'units', 'normalized', 'outerposition', [0 0 1 1]);

ellip_plotted = false(1, N);
opt_plotted   = false(1, N);
tauc_plotted  = false(1, N);
alpha_plotted = false(1, N);

%% --- 2. MAIN PROCESSING LOOP ---
for i = 1:N
    has_opt = ~isempty(opt_files{i});
    has_ellip = ~isempty(ellip_files{i});
    
    if has_ellip
        [L_ellip, psi, delta] = extract_ellip(ellip_files{i});
        ellip_plotted(i) = true;
        
        set(0, 'CurrentFigure', h_fig_ellip);
        subplot(1,2,1); hold on; plot(L_ellip, psi, 'LineWidth', 1.5, 'Color', cmap(i,:));
        subplot(1,2,2); hold on; plot(L_ellip, delta, 'LineWidth', 1.5, 'Color', cmap(i,:));
    end
    
    if has_opt
        [L_opt, n_idx, k_ext, alpha, E_ev] = extract_opt(opt_files{i});
        opt_plotted(i) = true;
        alpha_plotted(i) = true; 
        
        set(0, 'CurrentFigure', h_fig_opt);
        subplot(2,2,1); hold on; plot(L_opt, alpha, 'LineWidth', 1.5, 'Color', cmap(i,:));
        subplot(2,2,2); hold on; plot(E_ev, max(alpha, 1), 'LineWidth', 1.5, 'Color', cmap(i,:)); 
        subplot(2,2,3); hold on; plot(L_opt, n_idx, 'LineWidth', 1.5, 'Color', cmap(i,:));
        subplot(2,2,4); hold on; plot(L_opt, k_ext, 'LineWidth', 1.5, 'Color', cmap(i,:));
        
        % Alpha (Log) vs. Energy (for poster)
        set(0, 'CurrentFigure', h_fig_alpha_log); hold on;
        plot(E_ev, max(alpha, 1), 'LineWidth', 2.5, 'Color', cmap(i,:)); 
        
        % --- TAUC ANALYSIS ---
        % Assuming direct bandgap material for PEC: (alpha * E)^2
        y_tauc = (alpha .* E_ev).^2;
        deriv = diff(smoothdata(y_tauc, 'gaussian', 5)) ./ diff(E_ev);
        E_mid = (E_ev(1:end-1) + E_ev(2:end)) / 2;
        
        % TAUC SEARCH LIMIT: Updated for ~2.1 eV expected bandgap
        search_range = (E_mid > 1.8 & E_mid < 2.5); 
        
        if any(search_range)
            tauc_plotted(i) = true;
            [~, max_idx_rel] = max(deriv(search_range));
            search_indices = find(search_range);
            center_idx = search_indices(max_idx_rel);
            fit_win = (center_idx - 8) : (center_idx + 8);
            fit_win = fit_win(fit_win > 0 & fit_win <= length(E_ev));
            
            p = polyfit(E_ev(fit_win), y_tauc(fit_win), 1);
            Eg = -p(2) / p(1);
            extracted_Eg(i) = Eg;
            
            x_ext = linspace(Eg, E_ev(center_idx) + 0.1, 50);
            y_ext = polyval(p, x_ext);
            
            % Plot to Combined Tauc Figure
            set(0, 'CurrentFigure', h_fig_tauc); hold on;
            plot(E_ev, y_tauc, 'LineWidth', 2, 'Color', cmap(i,:));
            plot(x_ext, y_ext, '--', 'Color', cmap(i,:), 'LineWidth', 1.2, 'HandleVisibility', 'off');
            plot(Eg, 0, 'x', 'Color', cmap(i,:), 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
            
            if show_individual_tauc
                fig_ind = figure('Color', 'w', 'Name', sprintf('Tauc Plot - %s', legend_labels{i}));
                hold on;
                plot(E_ev, y_tauc, 'b', 'LineWidth', 2);
                plot(x_ext, y_ext, 'r--', 'LineWidth', 2);
                plot(Eg, 0, 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
                grid on; box on; xlim([1.8 2.5]);
                xlabel('Energy (eV)'); ylabel('(\alpha E)^2');
                title(sprintf('Tauc Plot: %s | Extracted E_g = %.4f eV', legend_labels{i}, Eg));
            end
        end
    end
end

%% --- 3. FINALIZE MAIN FIGURES ---
set(0, 'CurrentFigure', h_fig_ellip);
subplot(1,2,1); grid on; box on; xlabel('Wavelength (nm)'); ylabel('\Psi (deg)'); title('Psi (\Psi)'); xlim([250 1000]);
subplot(1,2,2); grid on; box on; xlabel('Wavelength (nm)'); ylabel('\Delta (deg)'); title('Delta (\Delta)'); xlim([250 1000]);
if any(ellip_plotted)
    h_leg1 = legend(legend_labels(ellip_plotted), 'NumColumns', 1);
    set(h_leg1, 'Position', [0.90 0.3 0.08 0.4], 'Units', 'normalized');
end

set(0, 'CurrentFigure', h_fig_opt);
subplot(2,2,1); grid on; box on; xlabel('Wavelength (nm)'); ylabel('\alpha (cm^{-1})'); title('Absorption vs Wavelength'); xlim([400 1000]);
subplot(2,2,2); grid on; box on; yscale("log"); xlabel('Energy (eV)'); ylabel('\alpha (cm^{-1})'); title('Absorption (Log)'); xlim([1.5 3.0]);
subplot(2,2,3); grid on; box on; xlabel('Wavelength (nm)'); ylabel('n'); title('n Index'); xlim([400 1000]);
subplot(2,2,4); grid on; box on; xlabel('Wavelength (nm)'); ylabel('k'); title('k Index'); xlim([400 1000]);
if any(opt_plotted)
    h_leg2 = legend(legend_labels(opt_plotted), 'NumColumns', 1);
    set(h_leg2, 'Position', [0.90 0.3 0.08 0.4], 'Units', 'normalized');
end

set(0, 'CurrentFigure', h_fig_tauc);
grid on; box on; xlabel('Energy (eV)'); ylabel('(\alpha E)^2'); title('Combined Tauc Plots'); 
xlim([1.8 2.5]); % Updated limits for PEC bandgap
if any(tauc_plotted)
    legend(legend_labels(tauc_plotted), 'Location', 'eastoutside');
end

%% --- 3.1 FINALIZE ALPHA LOG POSTER FIGURE ---
set(0, 'CurrentFigure', h_fig_alpha_log);
grid on; box on;
yscale("log"); 
xlabel('Energy (eV)', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('\alpha (cm^{-1})', 'FontSize', 24, 'FontWeight', 'bold');
title('Absorption Coefficient (Log Scale): 20 C to -110 C', 'FontSize', 28, 'FontWeight', 'bold');
xlim([1.7 2.8]); % Adjusted to center the absorption edge
if any(alpha_plotted)
    h_leg_alpha = legend(legend_labels(alpha_plotted), 'NumColumns', 1, 'FontSize', 12);
    set(h_leg_alpha, 'Location', 'eastoutside');
end
set(gca, 'FontSize', 16);

%% --- 4. EXPERIMENTAL BANDGAP VS. TEMPERATURE ---
h_fig_eg_trend = figure('Color', 'w', 'Name', 'Extracted Bandgap Trend', 'units', 'normalized', 'outerposition', [0.2 0.2 0.6 0.6]);
hold on; 

valid_idx = ~isnan(extracted_Eg);
valid_temps = temps_C(valid_idx);
valid_Eg = extracted_Eg(valid_idx);

if any(valid_idx)
    plot(valid_temps, valid_Eg, '-ob', 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Experimental E_g (Tauc)');
    
    % Annotate points
    for i = 1:length(valid_temps)
        text(valid_temps(i), valid_Eg(i) + 0.002, sprintf('%.3f', valid_Eg(i)), ...
            'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end

grid on; box on; 
xlabel('Temperature (C)', 'FontSize', 16, 'FontWeight', 'bold'); 
ylabel('Band Gap Energy (eV)', 'FontSize', 16, 'FontWeight', 'bold'); 
title('Extracted Bandgap vs. Temperature', 'FontSize', 20, 'FontWeight', 'bold');
xlim([-120 30]);
set(gca, 'FontSize', 14);
legend('Location', 'best', 'FontSize', 14);

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

%% --- 6. EXPORT RESULTS ---
save_filename = 'results_PEC_bandgap.mat'; 
save(save_filename, 'temps_C', 'extracted_Eg');
fprintf('Results successfully saved to %s\n', save_filename);

%% --- 7. AUTO-SAVE ALL FIGURES ---
fprintf('\nSaving all figures... ');
save_folder = 'Saved_Figures_PEC';
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
    
    safe_filename = regexprep(fig_title, '[\\/:*?"<>|]', '_'); 
    safe_filename = strrep(safe_filename, ' ', '_');
    
    fig_path = fullfile(save_folder, [safe_filename '.fig']);
    savefig(current_fig, fig_path);
end
fprintf('Done! %d figures saved to folder "%s".\n', length(all_figs), save_folder);

%% =========================================================================
%                  INTERNAL DATA EXTRACTION FUNCTIONS
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