% =========================================================================
% MASTER SCRIPT: GaAs Optical Characterization (Unified Web Export V3)
% =========================================================================
clear all; close all; clc;
%% --- GLOBAL PARAMETERS ---
d_thick = 350e-4;       
n_doping = 3.8e18;      
temps_C = 20:-5:-110; 
temps_K = temps_C + 273.15; 
N = length(temps_C);
legend_labels = arrayfun(@(t) sprintf('%d °C', t), temps_C, 'UniformOutput', false);
cmap = flipud(jet(N)); 
extracted_Eg = nan(1, N);

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

save_folder = 'Web_Export';
if ~exist(save_folder, 'dir'), mkdir(save_folder); end
sample_name = 'gaas_n36'; 

%% --- EXTRACT 20C BASELINE FOR LIVE OVERLAY ---
idx_20 = find(temps_C == 20);
has_20_ellip = false; has_20_opt = false;
if ~isempty(idx_20)
    if ~isempty(ellip_files{idx_20})
        [L_ellip_20, psi_20, delta_20] = extract_ellip(ellip_files{idx_20});
        has_20_ellip = true;
    end
    if ~isempty(opt_files{idx_20})
        [L_opt_20, n_idx_20, k_ext_20, alpha_20, E_ev_20] = extract_opt(opt_files{idx_20});
        has_20_opt = true;
    end
end

%% --- INITIALIZE GLOBAL FIGURES ---
h_fig_tauc = figure('Visible', 'off'); hold on;
h_fig_eg_model = figure('Visible', 'off'); hold on;
h_fig_global_alpha = figure('Visible', 'off'); hold on;

opt_plotted = false(1, N);
ellip_plotted = false(1, N);
tauc_plotted = false(1, N);

%% --- MAIN LOOP ---
for i = 1:N
    has_opt = ~isempty(opt_files{i});
    has_ellip = ~isempty(ellip_files{i});
    
    if temps_C(i) >= 0, t_str = sprintf('%dC', temps_C(i)); else, t_str = sprintf('m%dC', abs(temps_C(i))); end
    prefix = fullfile(save_folder, sprintf('%s_%s', sample_name, t_str));
    
    if has_ellip
        [L_ellip, psi, delta] = extract_ellip(ellip_files{i});
        ellip_plotted(i) = true;
        
        fig_temp_ellip = figure('Visible', 'off');
        subplot(2,1,1); hold on; grid on; box on;
        if has_20_ellip, plot(L_ellip_20, psi_20, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2); end
        plot(L_ellip, psi, 'b', 'LineWidth', 2); title('Psi (\Psi) Spectrum'); xlim([300 1000]);
        
        subplot(2,1,2); hold on; grid on; box on;
        if has_20_ellip, plot(L_ellip_20, delta_20, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2); end
        plot(L_ellip, delta, 'r', 'LineWidth', 2); title('Delta (\Delta) Spectrum'); xlim([300 1000]);
        saveas(fig_temp_ellip, [prefix '_psi.png']);
        close(fig_temp_ellip);
    end
    
    if has_opt
        [L_opt, n_idx, k_ext, alpha, E_ev] = extract_opt(opt_files{i});
        opt_plotted(i) = true;
        
        % Global Alpha collection
        set(0, 'CurrentFigure', h_fig_global_alpha);
        plot(L_opt, alpha, 'LineWidth', 1.5, 'Color', cmap(i,:));
        
        % Individual n,k export
        fig_temp_nk = figure('Visible', 'off');
        subplot(2,1,1); hold on; grid on; box on;
        if has_20_opt, plot(L_opt_20, n_idx_20, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2); end
        plot(L_opt, n_idx, 'g', 'LineWidth', 2); title('Refractive Index (n)');
        
        subplot(2,1,2); hold on; grid on; box on;
        if has_20_opt, plot(L_opt_20, k_ext_20, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2); end
        plot(L_opt, k_ext, 'm', 'LineWidth', 2); title('Extinction Coefficient (k)');
        saveas(fig_temp_nk, [prefix '_nk.png']);
        close(fig_temp_nk);
        
        % Individual Alpha Export with Live Overlay
        fig_temp_alpha = figure('Visible', 'off'); hold on; grid on; box on;
        if has_20_opt, plot(L_opt_20, alpha_20, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5); end
        plot(L_opt, alpha, 'Color', [0 0.6 0.3], 'LineWidth', 2);
        xlabel('Wavelength (nm)'); ylabel('\alpha (cm^{-1})'); title(sprintf('Absorption Coefficient (\\alpha) | %d°C', temps_C(i)));
        xlim([400 1000]); saveas(fig_temp_alpha, [prefix '_alpha.png']);
        close(fig_temp_alpha);
        
        % Tauc calculations
        y_tauc = (alpha .* E_ev).^2;
        deriv = diff(smoothdata(y_tauc, 'gaussian', 5)) ./ diff(E_ev);
        E_mid = (E_ev(1:end-1) + E_ev(2:end)) / 2;
        search_range = (E_mid > 1.35 & E_mid < 1.70); 
        
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
            
            x_ext = linspace(Eg, E_ev(center_idx) + 0.05, 50);
            y_ext = polyval(p, x_ext);
            
            set(0, 'CurrentFigure', h_fig_tauc);
            plot(E_ev, y_tauc, 'LineWidth', 2, 'Color', cmap(i,:));
            
            fig_temp_tauc = figure('Visible', 'off'); hold on; grid on; box on;
            plot(E_ev, y_tauc, 'k', 'LineWidth', 2);
            plot(x_ext, y_ext, 'r--', 'LineWidth', 2);
            plot(Eg, 0, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
            title(sprintf('Tauc Plot | Eg = %.3f eV', Eg)); xlim([1.35 1.70]);
            saveas(fig_temp_tauc, [prefix '_tauc.png']);
            close(fig_temp_tauc);
        end
    end
end

%% --- FINALIZE & SAVE GLOBAL FIGURES ---
set(0, 'CurrentFigure', h_fig_tauc);
grid on; box on; xlabel('Energy (eV)'); ylabel('(\alpha E)^2'); title('Combined Tauc Plots (20°C to -110°C)'); xlim([1.35 1.70]);
saveas(h_fig_tauc, fullfile(save_folder, 'global_tauc.png'));

set(0, 'CurrentFigure', h_fig_global_alpha);
grid on; box on; xlabel('Wavelength (nm)'); ylabel('\alpha (cm^{-1})'); title('Global Absorption Coefficient (\alpha) vs. Wavelength'); xlim([400 1000]);
saveas(h_fig_global_alpha, fullfile(save_folder, 'global_alpha.png'));

set(0, 'CurrentFigure', h_fig_eg_model);
T_range = linspace(25, -115, 200) + 273.15;
Eg_0 = 1.519; alpha_v = 5.405e-4; beta_v = 204;
Eg_intrinsic = Eg_0 - (alpha_v .* T_range.^2) ./ (T_range + beta_v);
m_e = 0.067; m_hh = 0.51; m_vc_star = (m_e * m_hh) / (m_e + m_hh) * 9.11e-31; 
Delta_BM = (( (1.054e-34)^2 / (2 * m_vc_star) ) * (3*pi^2 * n_doping*1e6)^(2/3)) / 1.602e-19; 
N_norm = n_doping / 1e18; Delta_BGN = (62 * (N_norm)^(1/3) + 7.4 * (N_norm)^(1/4)) / 1000; 
Eg_theory = Eg_intrinsic + Delta_BM - Delta_BGN;
plot(T_range-273.15, Eg_intrinsic, '--k', 'LineWidth', 1.5, 'DisplayName', 'Varshni (Intrinsic)');
plot(T_range-273.15, Eg_theory, '-k', 'LineWidth', 2, 'DisplayName', 'Jain 1992 (Doped)');
valid_mask = ~isnan(extracted_Eg) & (temps_C >= -70);
plot(temps_C(valid_mask), extracted_Eg(valid_mask), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Experiment (Tauc)');
grid on; box on; xlabel('Temperature (°C)'); ylabel('Band Gap Energy (eV)'); title('GaAs Bandgap Model vs. Experiment');
legend('Location', 'northeast'); xlim([-80 25]);
saveas(h_fig_eg_model, fullfile(save_folder, 'global_bandgap_model.png'));

%% --- UNIFIED NOISE ANALYSIS (2 CLEAR PLOTS) ---
filename_noise = 'ntype36_vacuum_nitro_raw_data.xlsx';
wl_noise_min = 250; wl_noise_max = 1000; 
noise_metric_delta = nan(1, N); noise_metric_psi = nan(1, N);
h_fig_delta = figure('Visible', 'off', 'Position', [100 100 900 400]);
ax_d1 = subplot(1,2,1); hold on; grid on; box on; xlabel('Wavelength (nm)'); ylabel('\Delta (deg)'); title('Delta Spectra Overlay');
ax_d2 = subplot(1,2,2); hold on; grid on; box on; xlabel('Temperature (°C)'); ylabel('RMS Noise'); title('Delta RMS Noise Trend');
h_fig_psi = figure('Visible', 'off', 'Position', [100 100 900 400]);
ax_p1 = subplot(1,2,1); hold on; grid on; box on; xlabel('Wavelength (nm)'); ylabel('\Psi (deg)'); title('Psi Spectra Overlay');
ax_p2 = subplot(1,2,2); hold on; grid on; box on; xlabel('Temperature (°C)'); ylabel('RMS Noise'); title('Psi RMS Noise Trend');

try 
    rawDataCell = readcell(filename_noise); [numRows, numCols] = size(rawDataCell);
    for i = 1:N
        col_start = (i - 1) * 5 + 1;
        if col_start + 4 <= numCols
            wl_col = rawDataCell(:, col_start + 1); psi_col = rawDataCell(:, col_start + 3); delta_col = rawDataCell(:, col_start + 4);
            wl_raw = []; psi_raw = []; delta_raw = [];
            for k = 1:numRows
                v_wl = wl_col{k}; v_p = psi_col{k}; v_d = delta_col{k};
                if isnumeric(v_wl) && isnumeric(v_p) && isnumeric(v_d) && ~isnan(v_wl) && ~isnan(v_p) && ~isnan(v_d)
                    wl_raw(end+1,1) = v_wl; psi_raw(end+1,1) = v_p; delta_raw(end+1,1) = v_d;
                end
            end
            idx = (wl_raw >= wl_noise_min) & (wl_raw <= wl_noise_max);
            if sum(idx) > 10
                w = wl_raw(idx); p = psi_raw(idx); d = delta_raw(idx);
                t_d = smoothdata(d, 'sgolay', 45); r_d = d - t_d; noise_metric_delta(i) = rms(r_d);
                t_p = smoothdata(p, 'sgolay', 45); r_p = p - t_p; noise_metric_psi(i) = rms(r_p);
                plot(ax_d1, w, d, 'Color', cmap(i,:)); plot(ax_p1, w, p, 'Color', cmap(i,:));
            end
        end
    end
    bar(ax_d2, temps_C, noise_metric_delta, 'FaceColor', 'flat'); if ~isempty(ax_d2.Children), ax_d2.Children(1).CData = flipud(cmap); end
    plot(ax_d2, temps_C, noise_metric_delta, '-ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); xlim(ax_d2, [-120 30]);
    bar(ax_p2, temps_C, noise_metric_psi, 'FaceColor', 'flat'); if ~isempty(ax_p2.Children), ax_p2.Children(1).CData = flipud(cmap); end
    plot(ax_p2, temps_C, noise_metric_psi, '-ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); xlim(ax_p2, [-120 30]);
    saveas(h_fig_delta, fullfile(save_folder, 'global_noise_delta.png')); saveas(h_fig_psi, fullfile(save_folder, 'global_noise_psi.png'));
catch
    warning('Noise file unreadable. Skipping noise plots.');
end
close all;
function [L, n, k, a, E] = extract_opt(filename)
    data = readtable(filename, 'VariableNamingRule', 'preserve'); raw = data{:, :};
    if iscell(raw), proc = nan(size(raw)); for r=1:size(raw,1), for c=1:size(raw,2), if isnumeric(raw{r,c}), proc(r,c)=raw{r,c}; else, proc(r,c)=str2double(string(raw{r,c})); end; end; end; raw = proc; end
    L = raw(:, 1); n = raw(:, 2); k = raw(:, 3); a = raw(:, 4); E = raw(:, 5);
end
function [L, psi, delta] = extract_ellip(filename)
    data = readtable(filename, 'VariableNamingRule', 'preserve'); raw = data{:, :};
    if iscell(raw), proc = nan(size(raw)); for r=1:size(raw,1), for c=1:size(raw,2), if isnumeric(raw{r,c}), proc(r,c)=raw{r,c}; else, proc(r,c)=str2double(string(raw{r,c})); end; end; end; raw = proc; end
    L = raw(:, 1); psi = raw(:, 3); delta = raw(:, 4); 
end