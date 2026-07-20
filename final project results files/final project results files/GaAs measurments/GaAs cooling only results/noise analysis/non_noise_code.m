% =========================================================================
% MASTER SCRIPT: Ellipsometry Noise Analysis (Psi & Delta)
% Evaluates environmental control efficacy (Vacuum)
% Range: 20°C to -70°C 
% =========================================================================
clear all; close all; clc;

%% === STEP 1: USER-DEFINED VARIABLES & SETUP ===
% Define the target raw data Excel file exported from WVASE software.
filename = 'ntype36_cooling_only_raw_data.xlsx';

% Define the temperature array: 20°C down to -70°C in 5-degree increments.
% 'N' stores the total number of temperature steps (27 steps total).
temps_C = 20:-5:-70; 
N       = length(temps_C);

% Define the optical wavelength range for noise evaluation (Visible to Near-IR).
% We exclude UV or Deep-IR if the sensor is noisy in those extreme edges.
wl_noise_min = 200;  
wl_noise_max = 1000; 

% --- PRE-ALLOCATION ---
% We pre-allocate arrays with 'nan' (Not-a-Number) to store our final 
% calculated RMS noise values. Pre-allocation makes MATLAB run faster and 
% prevents memory fragmentation during the loop.
noise_metric_delta = nan(1, N);
noise_metric_psi   = nan(1, N);

% --- COLORMAP SETUP ---
% 'jet' creates a spectrum of colors. 'flipud' reverses it so that 
% hot temperatures (20°C) are Red, and cold temperatures (-110°C) are Blue.
cmap = flipud(jet(N));

%% === STEP 2: INITIALIZE FIGURES ===
% Create the visual layout for the presentation matrices and summaries.
% We use 'units', 'normalized' so the windows scale perfectly on any monitor.

% Figure 1: Global Summary (Overlay of all spectra and RMS trends)
h_fig_main = figure('Color', 'w', 'Name', 'Figure 1: Global Summary', 'units', 'normalized', 'outerposition', [0.05 0.5 0.4 0.45]);
ax1 = subplot(2, 2, 1); hold on; grid on; box on;
xlabel('Wavelength (nm)'); ylabel('Delta Unwrapped (deg)'); title('Delta Spectra Overlay');
ax2 = subplot(2, 2, 2); hold on; grid on; box on;
xlabel('Temperature (C)'); ylabel('RMS Noise'); title('Delta Noise Trend');
ax3 = subplot(2, 2, 3); hold on; grid on; box on;
xlabel('Wavelength (nm)'); ylabel('Psi (deg)'); title('Psi Spectra Overlay');
ax4 = subplot(2, 2, 4); hold on; grid on; box on;
xlabel('Temperature (C)'); ylabel('RMS Noise'); title('Psi Noise Trend');

% Figure 2: Delta Residual Matrix (Grid showing isolated noise at each temp)
h_fig_delta_grid = figure('Color', 'w', 'Name', 'Figure 2: Delta Noise Matrix', 'units', 'normalized', 'outerposition', [0.5 0.05 0.45 0.9]);

% Figure 3: Psi Residual Matrix (Grid showing isolated noise at each temp)
h_fig_psi_grid = figure('Color', 'w', 'Name', 'Figure 3: Psi Noise Matrix', 'units', 'normalized', 'outerposition', [0.05 0.05 0.45 0.4]);

%% === STEP 3: DATA EXTRACTION ===
% We use 'readcell' instead of 'readtable' or 'xlsread' because the raw WVASE 
% Excel file often contains mixed data types (text headers, empty cells, and numbers). 
% 'readcell' safely loads everything into a flexible Cell Array.

fprintf('Loading raw data from %s...\n', filename);
try
    rawDataCell = readcell(filename);
catch
    error('Could not read %s. Please ensure the file is closed in Excel.', filename);
end
[numRows, numCols] = size(rawDataCell);

%% === STEP 4: MAIN PROCESSING & NOISE CALCULATION LOOP ===
% LOOP PURPOSE:
% Iterates through each temperature block in the Excel file. It prevents phase 
% jumps using 'unwrap', and uses a Savitzky-Golay filter to establish a theoretical 
% baseline. The RMS of the difference (residuals) represents the scattering noise.

fprintf('Analyzing %d measurement blocks...\n\n', N);

for i = 1:N
    % --- COLUMN MATH EXPLANATION ---
    % WVASE software exports data in blocks. Each temperature block takes up 
    % exactly 5 columns in the Excel file (e.g., Wavelength, AoI, Psi, Delta, Error).
    % Therefore, block 1 starts at col 1, block 2 starts at col 6, block 3 at 11, etc.
    col_start = (i - 1) * 5 + 1;
    
    % Ensure we don't exceed the actual number of columns in the Excel file
    if col_start + 4 <= numCols
        wl_col    = rawDataCell(:, col_start + 1); % Wavelength column
        psi_col   = rawDataCell(:, col_start + 3); % Psi column
        delta_col = rawDataCell(:, col_start + 4); % Delta column
        
        % Temporary empty arrays to store the clean, numeric data
        wl_raw = []; psi_raw = []; delta_raw = [];
        
        % --- ROW PARSER LOOP (DATA CLEANING) ---
        % This loop goes down row by row. It acts as a strict filter to ignore 
        % text headers (like "nm" or "Psi") and empty cells, extracting ONLY valid numbers.
        for k = 1:numRows
            v_wl = wl_col{k}; 
            v_p  = psi_col{k}; 
            v_d  = delta_col{k};
            
            % Check if all three values in the current row are pure, non-NaN numbers
            if isnumeric(v_wl) && isnumeric(v_p) && isnumeric(v_d) && ...
               ~isnan(v_wl) && ~isnan(v_p) && ~isnan(v_d)
                
                % If valid, append them to our clean arrays
                wl_raw(end+1,1)    = v_wl;
                psi_raw(end+1,1)   = v_p;
                delta_raw(end+1,1) = v_d;
            end
        end
        
        % --- WAVELENGTH FILTERING ---
        % Create a logical mask (True/False) for values within our requested range
        idx = (wl_raw >= wl_noise_min) & (wl_raw <= wl_noise_max);
        
        % Safety check: Only proceed if we have more than 10 valid data points.
        % This prevents the Savitzky-Golay filter from crashing on empty/small arrays.
        if sum(idx) > 10
            % Slice the arrays to include only the requested data points
            w = wl_raw(idx); 
            p = psi_raw(idx); 
            d = delta_raw(idx);
            
            % --- 4A. PHASE CORRECTION (UNWRAP) ---
            % Delta often "wraps" around 360/180 degrees. If we calculate noise 
            % over a raw jump, the RMS will falsely spike. Unwrap fixes this artifact.
            d_unwrapped = unwrap(d * pi / 180) * 180 / pi;
            
            % --- 4B. ISOLATE NOISE (SAVITZKY-GOLAY FILTERING) ---
            % WHY SAVITZKY-GOLAY? Unlike a simple moving average that flattens peaks, 
            % the S-G filter ('sgolay') fits a local polynomial to the data points. 
            % This preserves the natural, physical peaks and valleys of the spectrum 
            % while filtering out the high-frequency random noise.
            % 
            % PARAMETER '45': This is the Window Size. It means the filter uses 45 
            % adjacent data points to calculate the smoothed value of a single point.
            % By subtracting this "ideal" smooth baseline from the raw data, we isolate 
            % the residuals (the actual random scattering noise).
            
            t_d = smoothdata(d_unwrapped, 'sgolay', 45); % Ideal theoretical Delta baseline
            r_d = d_unwrapped - t_d;                     % Residuals (Raw minus Ideal)
            noise_metric_delta(i) = rms(r_d);            % Calculate RMS on the residuals
            
            t_p = smoothdata(p, 'sgolay', 45);           % Ideal theoretical Psi baseline
            r_p = p - t_p;                               % Residuals (Raw minus Ideal)
            noise_metric_psi(i) = rms(r_p);              % Calculate RMS on the residuals
            
            % --- 4C. PLOT SPECTRA ON MAIN SUMMARY (Fig 1) ---
            set(0, 'CurrentFigure', h_fig_main);
            plot(ax1, w, d_unwrapped, 'Color', cmap(i,:)); 
            plot(ax3, w, p, 'Color', cmap(i,:));
            
            % --- 4D. PLOT DELTA RESIDUALS MATRIX (Fig 2) ---
            set(0, 'CurrentFigure', h_fig_delta_grid);
            subplot(5, 6, i); hold on; grid on; box on;
            plot(w, r_d, 'Color', cmap(i,:), 'LineWidth', 0.5);
            yline(0, 'k-', 'LineWidth', 1);
            title(sprintf('%d C (RMS: %.3f)', temps_C(i), noise_metric_delta(i)), 'FontSize', 8);
            if i > 24, xlabel('nm'); end % Add X-labels only to the bottom row to reduce clutter
            xlim([wl_noise_min, wl_noise_max]);
            ylim([-7 7]); 
            
            % --- 4E. PLOT PSI RESIDUALS MATRIX (Fig 3) ---
            set(0, 'CurrentFigure', h_fig_psi_grid);
            subplot(5, 6, i); hold on; grid on; box on;
            plot(w, r_p, 'Color', cmap(i,:), 'LineWidth', 0.5);
            yline(0, 'k-', 'LineWidth', 1);
            title(sprintf('%d C (RMS: %.3f)', temps_C(i), noise_metric_psi(i)), 'FontSize', 8);
            if i > 24, xlabel('nm'); end
            xlim([wl_noise_min, wl_noise_max]);
            ylim([-3 3]);
            
            fprintf('  Processed %4d°C | Delta RMS: %.4f | Psi RMS: %.4f\n', temps_C(i), noise_metric_delta(i), noise_metric_psi(i));
        end
    end
end

%% === STEP 5: FINALIZE SUMMARY BARS (Fig 1) ===
% Populate the final summary bar charts with the calculated RMS values.

set(0, 'CurrentFigure', h_fig_main);

% --- Delta Bar Chart ---
bar(ax2, temps_C, noise_metric_delta, 'FaceColor', 'flat');
% Flip colormap specifically for the bars so -110°C is Blue and 20°C is Red
if ~isempty(ax2.Children), ax2.Children(1).CData = flipud(cmap); end
plot(ax2, temps_C, noise_metric_delta, '-ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');

% --- Psi Bar Chart ---
bar(ax4, temps_C, noise_metric_psi, 'FaceColor', 'flat');
if ~isempty(ax4.Children), ax4.Children(1).CData = flipud(cmap); end
plot(ax4, temps_C, noise_metric_psi, '-ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');

fprintf('\nAnalysis complete. %d measurement blocks mapped successfully.\n', N);