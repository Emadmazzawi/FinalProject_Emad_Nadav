%% --- ADVANCED NOISE ANALYSIS: NORMALIZED RMS vs. TEMP (POSTER - 4 METHODS) ---
clear all; close all; clc;

% --- 1. SETTINGS: FILES AND LABELS ---
% Removed ONLY Vacuum 2. Kept the other 4.
filenames = {
    'ntype36_vacuum1_raw_data.xlsx', ...
    'ntype36_nitrogin_only_raw_data.xlsx', ...
    'ntype36_cooling_only_raw_data.xlsx', ...
    'ntype36_vacuum_nitro_raw_data.xlsx'
};
labels = {'Vacuum 1', 'Nitrogen', 'Cooling Only', 'Vacuum + Nitro'}; 
N = length(filenames);

% --- WAVELENGTH RANGE ---
wl_min = 200; 
wl_max = 1000;

% Define the Temperature Range
target_temps = -110:5:20; 
num_T = length(target_temps);

% Initialize RMS Matrices
rms_delta_all    = nan(N, num_T); 
rms_psi_all      = nan(N, num_T); 
rms_combined_all = nan(N, num_T); 

% Colors and Markers (Green, Red, Blue, Purple)
colors = [0 0.7 0; 1 0 0; 0 0.4 1; 0.6 0 0.8];
markers = {'o', '^', 'v', 'd'};

% --- 2. DATA PROCESSING ---
for i = 1:N
    current_file = filenames{i};
    fprintf('\n----------------------------------------\n');
    fprintf('Processing file: %s\n', current_file);
    
    if ~isfile(current_file)
        warning('File not found: %s. Skipping...', current_file);
        continue;
    end
    
    rawDataCell = readcell(current_file);
    [numRows, numCols] = size(rawDataCell);
    
    num_blocks = floor(numCols / 5);
    
    for b = 0 : num_blocks - 1
        col_temp  = b * 5 + 1; 
        col_wl    = b * 5 + 2;
        col_psi   = b * 5 + 4;
        col_delta = b * 5 + 5;
        
        % --- ROBUST TEMPERATURE EXTRACTION ---
        T_vals = [];
        for k = 1:min(500, numRows)
            val = rawDataCell{k, col_temp};
            if isnumeric(val) && ~isnan(val)
                T_vals(end+1) = val;
            elseif ischar(val) || isstring(val)
                str_val = lower(string(val)); 
                str_val = regexprep(str_val, 'm\s*(\d+\.?\d*)', '-$1');
                str_val = replace(str_val, char(8211), '-');
                str_val = replace(str_val, char(8722), '-');
                nums = regexp(str_val, '-?\d+\.?\d*', 'match');
                if ~isempty(nums)
                    parsed = str2double(nums{end});
                    if ~isnan(parsed)
                        T_vals(end+1) = parsed;
                    end
                end
            end
        end
        
        T_vals = T_vals(T_vals >= -150 & T_vals <= 100);
        if isempty(T_vals), continue; end
        
        block_temp = mode(T_vals);
        [min_diff, target_idx] = min(abs(target_temps - block_temp));
        
        if min_diff <= 2.5 
            W_raw = []; P_raw = []; D_raw = [];
            for k = 1:numRows
                v_w = rawDataCell{k, col_wl};
                v_p = rawDataCell{k, col_psi};
                v_d = rawDataCell{k, col_delta};
                
                if isnumeric(v_w) && isnumeric(v_p) && isnumeric(v_d) ...
                   && ~isnan(v_w) && ~isnan(v_p) && ~isnan(v_d)
                    W_raw(end+1,1) = v_w;
                    P_raw(end+1,1) = v_p;
                    D_raw(end+1,1) = v_d;
                end
            end
            
            valid = (W_raw >= wl_min) & (W_raw <= wl_max);
            w = W_raw(valid); p = P_raw(valid); d = D_raw(valid);
            
            if length(w) > 10
                % --- NOISE MATH (RMS) ---
                trend_d = smoothdata(d, 'sgolay', 45);
                rms_delta = rms(d - trend_d);
                rms_delta_all(i, target_idx) = rms_delta;
                
                trend_p = smoothdata(p, 'sgolay', 45);
                rms_psi = rms(p - trend_p);
                rms_psi_all(i, target_idx) = rms_psi;
                
                % --- COMBINED NOISE SCORE ---
                rms_combined_all(i, target_idx) = (rms_delta + rms_psi) / 2;
            end
        end
    end
end

% --- 2.5 BASELINE NORMALIZATION (RELATIVE TO 20C) ---
fprintf('\nApplying Baseline Normalization (subtracting 20C noise)...\n');
idx_20C = find(target_temps == 20);

rms_delta_norm    = nan(N, num_T);
rms_psi_norm      = nan(N, num_T);
rms_combined_norm = nan(N, num_T);

for i = 1:N
    base_d = rms_delta_all(i, idx_20C);
    base_p = rms_psi_all(i, idx_20C);
    base_c = rms_combined_all(i, idx_20C);
    
    % Subtract the baseline from all measurements in this run
    rms_delta_norm(i, :)    = rms_delta_all(i, :) - base_d;
    rms_psi_norm(i, :)      = rms_psi_all(i, :) - base_p;
    rms_combined_norm(i, :) = rms_combined_all(i, :) - base_c;
end

% --- 3. PLOT NORMALIZED NOISE VS. TEMPERATURE ---
h_fig_temp_delta = figure('Color', 'w', 'Name', 'Normalized Delta Noise', 'units', 'normalized', 'outerposition', [0.05 0.2 0.3 0.6]);
hold on; grid on; box on; set(gca, 'FontSize', 16); 

h_fig_temp_psi = figure('Color', 'w', 'Name', 'Normalized Psi Noise', 'units', 'normalized', 'outerposition', [0.35 0.2 0.3 0.6]);
hold on; grid on; box on; set(gca, 'FontSize', 16); 

h_fig_temp_comb = figure('Color', 'w', 'Name', 'Normalized Combined Noise', 'units', 'normalized', 'outerposition', [0.65 0.2 0.3 0.6]);
hold on; grid on; box on; set(gca, 'FontSize', 16); 

for i = 1:N
    % Plot Normalized Delta
    set(0, 'CurrentFigure', h_fig_temp_delta);
    plot(target_temps, rms_delta_norm(i, :), ['-' markers{i}], 'Color', colors(i,:), ...
        'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', colors(i,:), 'DisplayName', labels{i});
    
    % Plot Normalized Psi
    set(0, 'CurrentFigure', h_fig_temp_psi);
    plot(target_temps, rms_psi_norm(i, :), ['-' markers{i}], 'Color', colors(i,:), ...
        'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', colors(i,:), 'DisplayName', labels{i});
        
    % Plot Normalized Combined
    set(0, 'CurrentFigure', h_fig_temp_comb);
    plot(target_temps, rms_combined_norm(i, :), ['-' markers{i}], 'Color', colors(i,:), ...
        'LineWidth', 4, 'MarkerSize', 11, 'MarkerFaceColor', colors(i,:), 'DisplayName', labels{i});
end

% Format Delta Figure 
set(0, 'CurrentFigure', h_fig_temp_delta);
xlabel('Temperature (°C)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Added RMS Error (rel. to 20°C)', 'FontSize', 20, 'FontWeight', 'bold'); 
title(sprintf('NORMALIZED Delta (\\Delta) RMS\n(%d-%d nm)', wl_min, wl_max), 'FontSize', 24, 'FontWeight', 'bold');
xlim([-115 25]); 
yline(0, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off'); 
legend('Location', 'northwest', 'FontSize', 18); 

% Format Psi Figure 
set(0, 'CurrentFigure', h_fig_temp_psi);
xlabel('Temperature (°C)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Added RMS Error (rel. to 20°C)', 'FontSize', 20, 'FontWeight', 'bold'); 
title(sprintf('NORMALIZED Psi (\\Psi) RMS\n(%d-%d nm)', wl_min, wl_max), 'FontSize', 24, 'FontWeight', 'bold');
xlim([-115 25]); 
yline(0, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');
legend('Location', 'northwest', 'FontSize', 18); 

% Format Combined Figure 
set(0, 'CurrentFigure', h_fig_temp_comb);
xlabel('Temperature (°C)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Added Average RMS (rel. to 20°C)', 'FontSize', 20, 'FontWeight', 'bold'); 
title(sprintf('NORMALIZED COMBINED NOISE\n(%d-%d nm)', wl_min, wl_max), 'FontSize', 24, 'FontWeight', 'bold');
xlim([-115 25]); 
yline(0, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');
legend('Location', 'northwest', 'FontSize', 18); 

fprintf('\nAnalysis complete. 4 NORMALIZED POSTER charts generated (Vacuum 2 removed).\n');

% =========================================================================
% NEW: 4-QUADRANT COMPARISON FIGURE (Delta & Psi per condition)
% =========================================================================
fprintf('\nGenerating 4-Quadrant Comparison Figure...\n');
h_fig_quad = figure('Color', 'w', 'Name', 'Delta & Psi Noise Comparison', 'units', 'normalized', 'outerposition', [0.1 0.1 0.8 0.8]);

for i = 1:N
    subplot(2, 2, i);
    hold on; grid on; box on; set(gca, 'FontSize', 12);
    
    % --- Left Y-Axis: Delta Noise ---
    yyaxis left
    plot(target_temps, rms_delta_all(i, :), '-o', 'Color', [0 0.4470 0.7410], 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', [0 0.4470 0.7410]);
    ylabel('\Delta RMS Noise (deg)', 'FontWeight', 'bold');
    
    % --- Right Y-Axis: Psi Noise ---
    yyaxis right
    plot(target_temps, rms_psi_all(i, :), '-^', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', [0.8500 0.3250 0.0980]);
    ylabel('\Psi RMS Noise (deg)', 'FontWeight', 'bold');
    
    % Format Subplot
    title(labels{i}, 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('Temperature (°C)', 'FontWeight', 'bold');
    xlim([-115 25]);
    
    % Add custom legend for clarity
    legend({'\Delta Noise', '\Psi Noise'}, 'Location', 'northwest', 'FontSize', 10);
end

% Add a main title for the entire figure
sgtitle(sprintf('Noise Trends Analysis per Measurement Condition\nWavelength Range: %d-%d nm', wl_min, wl_max), 'FontSize', 20, 'FontWeight', 'bold');