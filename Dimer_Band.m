%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [March 2025]
    Description:  [Spectral plot gauge capacitance Dimer]
    --------------------------------------------------------------
%}

clear all;
close all;

% --- Define fixed parameters ---
    gamma = 3;          % Gauge potential
    delta = 0.001;      % Contrast parameter
    ell   = [1, 1];     % Resonator lengths
    s     = [1, 2];     % Resonator spacings
    N     = 300;        % Number of plotting points
    fs    = 18;         % Fontsize of the plot annotation 
    lw    = 2.5;        % Width of the bands

    search_range = 300000;

    % --- Normalise the unit cell to L = 1 ---
    s   = (0.5 / (s(1)   + s(2)))   * s;
    ell = (0.5 / (ell(1) + ell(2))) * ell;

    L = sum(ell) + sum(s);

% --- Generate the spectral bands ---
    % --- Compute the constant decay within the bands ---
    beta_fixed = sum(ell) * gamma / 2 ; 

    % --- initialise plotting range for band/gap functions ---
    alpha_vals   = linspace(-pi, pi, N);
    beta_vals_1  = linspace(-3 + beta_fixed,  3 + beta_fixed, N);
    beta_vals_2_fine  = linspace(- 2 + beta_fixed, 2 + beta_fixed, search_range);
    
    band_functions = zeros(2, N);

    for k = 1:N
        alpha = alpha_vals(k);
        C = quasiperiodic_capacitance_matrix(alpha, beta_fixed, gamma, ell, s);
        ev = eig(C); 
        band_functions(:, k) = sqrt( delta * real(ev) );
    end
    
    % --- Mark the limit of the spectrum ---
    C = quasiperiodic_capacitance_matrix(pi, beta_fixed, gamma, ell, s);
    lower_gaps = real(sqrt( delta * eig(C)));
    C = quasiperiodic_capacitance_matrix(0,  beta_fixed, gamma, ell, s);
    upper_gaps = real(sqrt( delta * eig(C)));

% --- Generate gap bands ---
    % --- Gap band for alpha = 0 ---
    alpha_0 = 0;
    gap_functions_1 = zeros(2, N);

    for k = 1:N
        beta_gap = beta_vals_1(k);
        C = quasiperiodic_capacitance_matrix(alpha_0, beta_gap, gamma, ell, s);
        ev = eig(C);
    
        gap_functions_1(:, k) = NaN(size(gap_functions_1, 1), 1);
    
        for j = 1:min(length(ev), size(gap_functions_1, 1))
            lambda = ev(j);
            if abs(imag(lambda)) < 1e-3 && real(lambda) >= - 1e-2
                gap_functions_1(j, k) = sqrt(delta * real(lambda));
        
            end
        end
    end
    
    % --- Gap band for alpha = pi ---
    alpha_pi  = pi;
    is_real_2 = false(1, search_range);  
    
    for k = 1:search_range
        beta_gap = beta_vals_2_fine(k);
        C = quasiperiodic_capacitance_matrix(alpha_pi, beta_gap, gamma, ell, s);
        ev = eig(C);
    
        for j = 1:length(ev)
            lambda = ev(j);
            if abs(imag(lambda)) < 1e-5 && real(lambda) >= -2e-1
                is_real_2(k) = true;
                break;
            end
        end
    end
    
    % --- extract interval where bands are real ---
    onset_indices_2  = find(diff([false, is_real_2]) == 1);
    offset_indices_2 = find(diff([is_real_2, false]) == -1);
    real_intervals_2 = [beta_vals_2_fine(onset_indices_2); beta_vals_2_fine(offset_indices_2)];
    
    % --- Gap band for alpha = pi ---
    alpha_pi = pi;
    gap_functions_2 = zeros(2, N);
    beta_vals_2  = linspace(real_intervals_2(1, 1), real_intervals_2(2, 1), N);
    
    for k = 1:N
        beta_gap = beta_vals_2(k);
        C = quasiperiodic_capacitance_matrix(alpha_pi, beta_gap, gamma, ell, s);
        ev = eig(C);
        gap_functions_2(:, k) = NaN(size(gap_functions_2, 1), 1);
    
        for j = 1:min(length(ev), size(gap_functions_2, 1))
            lambda = ev(j);
            if abs(imag(lambda)) < 1e-2 && real(lambda) >= - 1e-2
                gap_functions_2(j, k) = sqrt(delta * real(lambda));
            end
        end
    end

% --- Plotting the spectral bands ---
    figure;
    plot(alpha_vals,  real(band_functions(1,:)),  'k-', 'LineWidth', lw); 
    hold on;
    plot(alpha_vals,  real(band_functions(2,:)),  'k-', 'LineWidth', lw); 
    plot(beta_vals_1, real(gap_functions_1(1,:)), 'r-', 'LineWidth', lw); 
    plot(beta_vals_1, real(gap_functions_1(2,:)), 'r-', 'LineWidth', lw); 
    plot(beta_vals_2, real(gap_functions_2(1,:)), 'r-', 'LineWidth', lw);
    plot(beta_vals_2, real(gap_functions_2(2,:)), 'r-', 'LineWidth', lw); 

    yline(upper_gaps(1), 'k--', 'LineWidth', 1);
    yline(lower_gaps(1), 'k--', 'LineWidth', 1);
    yline(upper_gaps(2), 'k--', 'LineWidth', 1);
    yline(lower_gaps(2), 'k--', 'LineWidth', 1); 

    plot([0 0],                   [0 upper_gaps(2)],              'k-', 'LineWidth', 1); 
    plot([pi pi],                 [lower_gaps(2) lower_gaps(1)],  'k-', 'LineWidth', 1);
    plot([-pi -pi],               [lower_gaps(2) lower_gaps(1)],  'k-', 'LineWidth', 1); 
    plot([0 0],                   [upper_gaps(1) 2],              'k-', 'LineWidth', 1); 
    plot([beta_fixed beta_fixed], [upper_gaps(2)  lower_gaps(2)], 'r-', 'LineWidth', 1);
    plot([beta_fixed beta_fixed], [lower_gaps(1)  upper_gaps(1)], 'r-', 'LineWidth', 1);

    xlabel('$\alpha$ and $\beta$ respectively', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\omega^{\alpha, \beta, \gamma}$' , 'Interpreter', 'latex', 'FontSize', fs);
    ylim([0, max(upper_gaps) * 1.2]);
    xticks([-pi/L, 0, pi/L]); 
    xticklabels({'$-\pi/L$', '$0$', '$\pi/L$'});
    set(gca, 'FontSize', fs+4, 'TickLabelInterpreter', 'latex');
    set(gcf, 'Position', [100, 100, 500, 400]); 


%% --- Define the quasiperiodic capacitance matrix ---

function C = quasiperiodic_capacitance_matrix(alpha, beta, gamma, ell, s)
    % Fixed-size version for 2x2 matrix with inputs:
    % alpha, beta, gamma - scalars
    % ell = [ell1, ell2]
    % s   = [s1, s2]

    L = 1; % sum(ell) + sum(s); % total length
    C = zeros(2, 2); % initialize

    % i = j = 1
    C(1,1) = (gamma/s(1)) * (ell(1)/(1 - exp(-gamma*ell(1)))) - (gamma/s(2)) * (ell(1)/(1 - exp(gamma*ell(1))));

    % i = j = 2
    C(2,2) = (gamma/s(2)) * (ell(2)/(1 - exp(-gamma*ell(2)))) - (gamma/s(1)) * (ell(2)/(1 - exp(gamma*ell(2))));

    % i = j + 1 => i = 2, j = 1
    C(2,1) = (gamma/s(1)) * (ell(1)/(1 - exp(gamma*ell(1))));

    % i = j - 1 => i = 1, j = 2
    C(1,2) = -(gamma/s(1)) * (ell(2)/(1 - exp(-gamma*ell(1))));

    % Corner correction
    
    % i = 2, j = 1
    C(2,1) = C(2,1) - exp(+1i*(alpha + 1i*beta)*L) * (gamma / s(2)) * ell(1)/(1-exp(-gamma* ell(1)));

    % i = 1, j = 2
    C(1,2) = C(1,2) + exp(-1i*(alpha + 1i*beta)*L) * (gamma / s(2)) * ell(2)/(1 - exp(gamma*ell(2)));
end
