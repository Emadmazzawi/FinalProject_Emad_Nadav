% =========================================================================
% SCRIPT: ADVANCED NOISE ANALYSIS (NORMALIZED RMS vs. TEMP)
% PURPOSE: Evaluates the environmental control methods by calculating the 
%          added RMS scattering noise as temperatures drop to -110°C.
%          Reproduces the Combined Noise analysis graph from the project.
% =========================================================================
clear all; close all; clc;

%% === STEP 1: SETTINGS, FILES, AND LABELS ===
% Define the 4 measurement setups we are comparing to justify the final design.
filenames = {
    'ntype36_vacuum1_raw_data.xlsx', ...
    'ntype36_nitrogin_only_raw_data.xlsx', ...
    'ntype36_cooling_only_raw_data.xlsx', ...
    'ntype36_vacuum_nitro_raw_data.xlsx'
};

% Presentation labels corresponding to the files
labels = {'Vacuum Only', 'Nitrogen Flow', 'Cooling Only', 'Vacuum + Nitro'}; 
N = length(filenames);

% --- WAVELENGTH RANGE ---
wl_min = 200; 
wl_max = 1000;

% Define the Temperature Range: 20°C down to -110°C in 5°C steps
target_temps = -110:5:20; 
num_T = length(target_temps);

% --- PRE-ALLOCATION ---
% Initialize matrices with 'nan' to hold the RMS noise for each condition and temp.
rms_delta_all    = nan(N, num_T); 
rms_psi_all      = nan(N, num_T); 
rms_combined_all = nan(N, num_T); 

% Colors (Green, Red, Blue, Purple) and Markers explicitly matched to presentation
colors = [0 0.7 0; 1 0 0; 0 0.4 1; 0.6 0 0.8];
markers = {'o', '^', 'v', 'd'};

%% === STEP 2: DATA PROCESSING AND NOISE EXTRACTION ===
% LOOP PURPOSE: 
% Reads raw data, robustly extracts the temperature from WVASE headers, applies 
% phase unwrap, and computes the Savitzky-Golay RMS noise for each block.

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
    
    % WVASE exports data in 5-column blocks per temperature measurement
    num_blocks = floor(numCols / 5);
    
    for b = 0 : num_blocks - 1
        col_temp  = b * 5 + 1; 
        col_wl    = b * 5 + 2;
        col_psi   = b * 5 + 4;
        col_delta = b * 5 + 5;
        
        % --- 2A. ROBUST TEMPERATURE EXTRACTION ---
        % WVASE Excel exports often have messy headers (e.g., text, weird hyphens, 
        % or string characters instead of clean minus signs). This regex block acts 
        % as a parser to find and convert the true numeric temperature of the block.
        T_vals = [];
        for k = 1:min(500, numRows)
            val = rawDataCell{k, col_temp};
            if isnumeric(val) && ~isnan(val)
                T_vals(end+1) = val;
            elseif ischar(val) || isstring(val)
                str_val = lower(string(val)); 
                str_val = regexprep(str_val, 'm\s*(\d+\.?\d*)', '-$1'); % Fix 'm100' to '-100'
                str_val = replace(str_val, char(8211), '-');            % Fix en-dash
                str_val = replace(str_val, char(8722), '-');            % Fix minus sign
                nums = regexp(str_val, '-?\d+\.?\d*', 'match');
                if ~isempty(nums)
                    parsed = str2double(nums{end});
                    if ~isnan(parsed)
                        T_vals(end+1) = parsed;
                    end
                end
            end
        end
        
        % Filter unrealistic numbers and find the most common value (mode) as the block's Temp
        T_vals = T_vals(T_vals >= -150 & T_vals <= 100);
        if isempty(T_vals), continue; end
        
        block_temp = mode(T_vals);
        
        % Map the extracted temperature to our predefined target array (within 2.5°C margin)
        [min_diff, target_idx] = min(abs(target_temps - block_temp));
        
        if min_diff <= 2.5 
            % --- 2B. DATA CLEANING ---
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
            
            % Slice to desired wavelength range
            valid = (W_raw >= wl_min) & (W_raw <= wl_max);
            w = W_raw(valid); p = P_raw(valid); d = D_raw(valid);
            
            if length(w) > 10
                % --- 2C. NOISE MATH (RMS) WITH SAVITZKY-GOLAY ---
                % Prevent artificial noise spikes from phase boundary wrapping
                d_unwrapped = unwrap(d * pi / 180) * 180 / pi;
                
                % Delta Noise: Subtract 45-point 2nd-order polynomial fit
                trend_d = smoothdata(d_unwrapped, 'sgolay', 45);
                rms_delta = rms(d_unwrapped - trend_d);
                rms_delta_all(i, target_idx) = rms_delta;
                
                % Psi Noise: Subtract 45-point 2nd-order polynomial fit
                trend_p = smoothdata(p, 'sgolay', 45);
                rms_psi = rms(p - trend_p);
                rms_psi_all(i, target_idx) = rms_psi;
                
                % Combined Noise Score (Average of Psi and Delta errors)
                rms_combined_all(i, target_idx) = (rms_delta + rms_psi) / 2;
            end
        end
    end
end

%% === STEP 3: BASELINE NORMALIZATION (RELATIVE TO 20°C) ===
% WHY NORMALIZE? Every machine has inherent electronic/optical noise even at 
% room temperature. By subtracting the noise measured at 20°C from all other 
% temperatures, we isolate ONLY the *added* physical scattering noise caused 
% by condensation or ice formation at low temperatures.

fprintf('\nApplying Baseline Normalization (subtracting 20°C noise)...\n');
idx_20C = find(target_temps == 20);

rms_delta_norm    = nan(N, num_T);
rms_psi_norm      = nan(N, num_T);
rms_combined_norm = nan(N, num_T);

for i = 1:N
    % Fetch the baseline noise at 20°C for the current condition
    base_d = rms_delta_all(i, idx_20C);
    base_p = rms_psi_all(i, idx_20C);
    base_c = rms_combined_all(i, idx_20C);
    
    % Subtract the baseline across the entire temperature array
    rms_delta_norm(i, :)    = rms_delta_all(i, :) - base_d;
    rms_psi_norm(i, :)      = rms_psi_all(i, :) - base_p;
    rms_combined_norm(i, :) = rms_combined_all(i, :) - base_c;
end

%% === STEP 4: PLOT NORMALIZED NOISE VS. TEMPERATURE ===
% Generate large, presentation-ready figures showing the noise divergence.

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
        
    % Plot Normalized Combined (Thicker line for final emphasis)
    set(0, 'CurrentFigure', h_fig_temp_comb);
    plot(target_temps, rms_combined_norm(i, :), ['-' markers{i}], 'Color', colors(i,:), ...
        'LineWidth', 4, 'MarkerSize', 11, 'MarkerFaceColor', colors(i,:), 'DisplayName', labels{i});
end

% --- FORMAT DELTA FIGURE ---
set(0, 'CurrentFigure', h_fig_temp_delta);
xlabel('Temperature (°C)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Added RMS Error (rel. to 20°C)', 'FontSize', 20, 'FontWeight', 'bold'); 
title(sprintf('NORMALIZED Delta (\\Delta) RMS\n(%d-%d nm)', wl_min, wl_max), 'FontSize', 24, 'FontWeight', 'bold');
xlim([-115 25]); ylim([-1 6]); % Consistent Y-limits based on expected data
yline(0, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off'); 
legend('Location', 'northwest', 'FontSize', 18); 

% --- FORMAT PSI FIGURE ---
set(0, 'CurrentFigure', h_fig_temp_psi);
xlabel('Temperature (°C)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Added RMS Error (rel. to 20°C)', 'FontSize', 20, 'FontWeight', 'bold'); 
title(sprintf('NORMALIZED Psi (\\Psi) RMS\n(%d-%d nm)', wl_min, wl_max), 'FontSize', 24, 'FontWeight', 'bold');
xlim([-115 25]); ylim([-1 6]);
yline(0, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');
legend('Location', 'northwest', 'FontSize', 18); 

% --- FORMAT COMBINED FIGURE ---
set(0, 'CurrentFigure', h_fig_temp_comb);
xlabel('Temperature (°C)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Added Average RMS (rel. to 20°C)', 'FontSize', 20, 'FontWeight', 'bold'); 
title(sprintf('NORMALIZED COMBINED NOISE\n(%d-%d nm)', wl_min, wl_max), 'FontSize', 24, 'FontWeight', 'bold');
xlim([-115 25]); ylim([-1 6]);
yline(0, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');
legend('Location', 'northwest', 'FontSize', 18); 

fprintf('\nAnalysis complete. 3 NORMALIZED POSTER charts generated.\n');