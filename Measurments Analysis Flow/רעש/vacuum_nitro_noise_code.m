% =========================================================================
% UNIFIED NOISE ANALYSIS: vacuum and nitrogen (20C to -110C)
% Evaluates RMS noise in Ellipsometry data (Psi and Delta) across temps.
% Includes Phase Unwrapping to prevent artificial noise spikes.
% =========================================================================

clear all; close all; clc;

%% --- 1. SETTINGS ---
filename = 'ntype36_vacuum_nitro_raw_data.xlsx';

% Create the temperature array: 20, 15, 10, ..., -110 (27 steps)
temps_C = 20:-5:-110; 
N = length(temps_C);

% Analysis wavelength range (Visible to NIR)
wl_noise_min = 200;  
wl_noise_max = 1000; 

% Metrics storage arrays
noise_metric_delta = nan(1, N);
noise_metric_psi   = nan(1, N);

% Smooth colormap transition from Red (Hot) to Blue (Cold)
% Index 1 (20C) will be Red, Index 27 (-110C) will be Blue
cmap = flipud(jet(N));

%% --- 2. SETUP FIGURES ---
% Figure 1: Global Summary
h_fig_main = figure('Color', 'w', 'Name', 'Figure 1: Global Summary', 'units', 'normalized', 'outerposition', [0.05 0.5 0.4 0.45]);

ax1 = subplot(2, 2, 1); hold on; grid on; box on;
xlabel('Wavelength (nm)'); ylabel('Delta Unwrapped (deg)'); title('Delta Spectra Overlay');

ax2 = subplot(2, 2, 2); hold on; grid on; box on;
xlabel('Temperature (C)'); ylabel('RMS Noise'); title('Delta Noise Trend');

ax3 = subplot(2, 2, 3); hold on; grid on; box on;
xlabel('Wavelength (nm)'); ylabel('Psi (deg)'); title('Psi Spectra Overlay');

ax4 = subplot(2, 2, 4); hold on; grid on; box on;
xlabel('Temperature (C)'); ylabel('RMS Noise'); title('Psi Noise Trend');

% Figure 2: Delta Residual Matrix (5x6 grid for 27 plots)
h_fig_delta_grid = figure('Color', 'w', 'Name', 'Figure 2: Delta Noise Matrix', 'units', 'normalized', 'outerposition', [0.5 0.05 0.45 0.9]);

% Figure 3: Psi Residual Matrix (5x6 grid for 27 plots)
h_fig_psi_grid = figure('Color', 'w', 'Name', 'Figure 3: Psi Noise Matrix', 'units', 'normalized', 'outerposition', [0.05 0.05 0.45 0.4]);

%% --- 3. DATA EXTRACTION & MATRIX PLOTTING ---
disp(['Analyzing ', num2str(N), ' measurement blocks...']);

try
    rawDataCell = readcell(filename);
catch
    error('Could not read %s. Please ensure the file is closed in Excel.', filename);
end

[numRows, numCols] = size(rawDataCell);

for i = 1:N
    % Column math for 5-column blocks (1, 6, 11, 16...)
    col_start = (i - 1) * 5 + 1;
    
    if col_start + 4 <= numCols
        wl_col    = rawDataCell(:, col_start + 1);
        psi_col   = rawDataCell(:, col_start + 3);
        delta_col = rawDataCell(:, col_start + 4);
        
        wl_raw = []; psi_raw = []; delta_raw = [];
        
        for k = 1:numRows
            v_wl = wl_col{k}; 
            v_p  = psi_col{k}; 
            v_d  = delta_col{k};
            
            if isnumeric(v_wl) && isnumeric(v_p) && isnumeric(v_d) && ...
               ~isnan(v_wl) && ~isnan(v_p) && ~isnan(v_d)
                wl_raw(end+1,1)    = v_wl;
                psi_raw(end+1,1)   = v_p;
                delta_raw(end+1,1) = v_d;
            end
        end
        
        % Filter by wavelength range
        idx = (wl_raw >= wl_noise_min) & (wl_raw <= wl_noise_max);
        
        if sum(idx) > 10
            w = wl_raw(idx); 
            p = psi_raw(idx); 
            d = delta_raw(idx);
            
            % --- FIX: UNWRAP DELTA TO PREVENT PHASE JUMP ARTIFACTS ---
            d_unwrapped = unwrap(d * pi / 180) * 180 / pi;
            
            % --- MATH: ISOLATE RESIDUALS (JITTER) ---
            t_d = smoothdata(d_unwrapped, 'sgolay', 45);
            r_d = d_unwrapped - t_d;
            noise_metric_delta(i) = rms(r_d);
            
            t_p = smoothdata(p, 'sgolay', 45);
            r_p = p - t_p;
            noise_metric_psi(i) = rms(r_p);
            
            % --- PLOT ON MAIN SUMMARY (Fig 1) ---
            set(0, 'CurrentFigure', h_fig_main);
            plot(ax1, w, d_unwrapped, 'Color', cmap(i,:)); 
            plot(ax3, w, p, 'Color', cmap(i,:));
            
            % --- PLOT ON DELTA MATRIX (Fig 2) ---
            set(0, 'CurrentFigure', h_fig_delta_grid);
            subplot(5, 6, i); hold on; grid on; box on;
            plot(w, r_d, 'Color', cmap(i,:), 'LineWidth', 0.5);
            yline(0, 'k-', 'LineWidth', 1);
            title(sprintf('%d C (RMS: %.3f)', temps_C(i), noise_metric_delta(i)), 'FontSize', 8);
            if i > 24, xlabel('nm'); end % Add X-labels only to the bottommost plots
            xlim([wl_noise_min, wl_noise_max]);
            ylim([-7 7]); 
            
            % --- PLOT ON PSI MATRIX (Fig 3) ---
            set(0, 'CurrentFigure', h_fig_psi_grid);
            subplot(5, 6, i); hold on; grid on; box on;
            plot(w, r_p, 'Color', cmap(i,:), 'LineWidth', 0.5);
            yline(0, 'k-', 'LineWidth', 1);
            title(sprintf('%d C (RMS: %.3f)', temps_C(i), noise_metric_psi(i)), 'FontSize', 8);
            if i > 24, xlabel('nm'); end
            xlim([wl_noise_min, wl_noise_max]);
            ylim([-3 3]);
            
            fprintf('Processed %d C (Delta RMS: %.4f)\n', temps_C(i), noise_metric_delta(i));
        end
    end
end

%% --- 4. FINALIZE SUMMARY BARS (Fig 1) ---
set(0, 'CurrentFigure', h_fig_main);

% Delta Bar
bar(ax2, temps_C, noise_metric_delta, 'FaceColor', 'flat');
% Flip colormap again just for the bars so -110 gets Blue and 20 gets Red
if ~isempty(ax2.Children), ax2.Children(1).CData = flipud(cmap); end
plot(ax2, temps_C, noise_metric_delta, '-ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');

% Psi Bar
bar(ax4, temps_C, noise_metric_psi, 'FaceColor', 'flat');
if ~isempty(ax4.Children), ax4.Children(1).CData = flipud(cmap); end
plot(ax4, temps_C, noise_metric_psi, '-ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');

disp(['Analysis complete. ', num2str(N), ' measurement blocks mapped successfully.']);