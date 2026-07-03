%% --- ISOLATED BANDGAP THEORY (VARSHNI & JAIN 1992) ---
% Parameters
n_doping = 3.8e18; % Doping concentration [cm^-3]
T_range_C = linspace(25, -115, 200);
T_range_K = T_range_C + 273.15;

% 1. Varshni Model for Intrinsic GaAs
Eg_0 = 1.519; 
alpha_v = 5.405e-4; 
beta_v = 204;
Eg_intrinsic = Eg_0 - (alpha_v .* T_range_K.^2) ./ (T_range_K + beta_v);

% 2. Jain (1992) Model Corrections for heavily doped n-GaAs
m_e = 0.067; m_hh = 0.51; 
m_vc_star = (m_e * m_hh) / (m_e + m_hh) * 9.11e-31; % Reduced effective mass in kg

% Burstein-Moss Shift (Delta E_BM) - Blue Shift
Delta_BM = (( (1.054e-34)^2 / (2 * m_vc_star) ) * (3*pi^2 * n_doping*1e6)^(2/3)) / 1.602e-19; 

% Bandgap Narrowing (Delta E_BGN) - Red Shift
N_norm = n_doping / 1e18;
Delta_BGN = (62 * (N_norm)^(1/3) + 7.4 * (N_norm)^(1/4)) / 1000; 

% Net Bandgap Model
Eg_theory = Eg_intrinsic + Delta_BM - Delta_BGN;

% Plotting
figure('Color', 'w', 'Name', 'Theoretical Bandgap');
hold on; 
plot(T_range_C, Eg_intrinsic, '--k', 'LineWidth', 1.5, 'DisplayName', 'Varshni (Intrinsic)');
plot(T_range_C, Eg_theory, '-r', 'LineWidth', 2, 'DisplayName', 'Jain 1992 (n-GaAs Doped)');
grid on; box on; 
xlabel('Temperature (°C)'); 
ylabel('Band Gap Energy (eV)'); 
title('Theoretical GaAs Bandgap Model');
legend('Location', 'best');
xlim([-115 25]);