%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [April 2025]
    Description:  [Defected Dimer system and Green's function]
    --------------------------------------------------------------
%}

clear all;
close all;

% --- Define the system parameters ---
    n     = 200;        % Number of unit cells
    gamma = 3;          % Gauge potential
    delta = 0.001;      % Contrast parameter
    nu    = 1.5;        % Change in wavespeed
    ell   = [1, 1];     % Resonator lengths
    spa   = [1, 2];     % Resonator spacings
    fs    = 18;         % Fontsize in plot annotation
    lw    = 2.5;        % Linewidth of the spectral bands

    search_range = 300000; 
    
    % --- Set the site of the defect ---
    defect_site = floor(n/3);

% --- Normalise the unit cell to L = 1 ---
    spa = (0.5 / (spa(1) + spa(2))) * spa;
    ell = (0.5 / (ell(1) + ell(2))) * ell;

    spa_fix = spa;
    ell_fix = ell;

    s = zeros(1, 2*n);
    l = zeros(1, 2*n);
    
    for i = 1:n
        idx = 2*(i-1) + 1;
        l(idx:idx+1) = [ell(1), ell(1)];
        s(idx:idx+1) = [spa(1), spa(2)];
    end

    N = 2 * n;
    
    % --- Generate the capacitance matrix ---
    capmat = Capacitance(N, s, gamma,  l);

    % --- Define defect Matrix ---
    D = eye(N);
    D(defect_site, defect_site) = 1 + nu;
    Def = D * capmat;
       
    % --- Compute the resonances ---
    [Vect, eigenvalues_matrix] = eig(Def);
    eigenvalues = sort(diag(eigenvalues_matrix));
    eigenvalues = abs(eigenvalues);
    resonances  = sqrt(eigenvalues * delta);

% --- Compute the limit of the spectrum ---

    beta_fixed = sum(ell) * gamma / 2;
    cap0  = quasiperiodic_capacitance_matrix(0,  beta_fixed, gamma, ell, spa);
    cappi = quasiperiodic_capacitance_matrix(pi, beta_fixed, gamma, ell, spa);
    
    lim1 = real(sqrt(delta * eig(cap0)) );
    lim2 = real(sqrt(delta * eig(cappi)));

% --- Retrieve the defect frequncies ---
    def = [];
    if nu > 0
        def(1) = resonances(n);
        def(2) = resonances(2*n);
        resonances(n)   = NaN;
        resonances(2*n) = NaN;
    else
        def(1) = resonances(2);
        def(2) = resonances(n + 1);
        resonances(2)   = NaN;
        resonances(n+1) = NaN;
    end

    disp('----------------------------------------');
    disp(['Defect of size eta= ', num2str(nu)]);
    disp(['Defect Frequency 1: ', num2str(def(1))]);
    disp(['Defect Frequency 2: ', num2str(def(2))]);


% --- Discrete Green's function ---
    cap = delta * Capacitance(N, s, gamma,  l);
    
    % --- Generate a point source ---
    dirac = zeros(1, N);
    dirac(floor(N/2)) = 1;
    dirac = dirac.';

    % --- Measure decay length of Green's function for frequency range ---
    N_freq = 30;
    w = linspace(0.001, 0.2, N_freq);

    decay_l = zeros(1, N_freq);
    decay_r = zeros(1, N_freq);
    I = eye(N);

    for i = 1:N_freq
        Green = (cap-w(i)^2*I) \ I;
        Green = Green * dirac;

        G_l = abs(Green(1:floor(N/2)));
        G_r = abs(Green(floor(N/2):N));

        logG_l = log(G_l);
        logG_r = log(G_r);
        
        x_l = 1:length(logG_l);
        x_r = 1:length(logG_r);

        p_l = polyfit(x_l, logG_l, 1);
        p_r = polyfit(x_r, logG_r, 1);

        decay_l(i) = - 2 * p_l(1); % factor 2 because 2 resonators in unit cell
        decay_r(i) = - 2 * p_r(1);

    end

    % Debugging
    %figure;
    %plot(decay_l, w);
    %hold on;
    %plot(decay_r, w);
    %hold off;

% --- Compute the spectral bands ---
    L = sum(ell_fix) + sum(spa_fix);

% --- Generate the spectral bands ---
    % --- Compute the constant decay within the bands ---
    beta_fixed =   0.5 * sum(ell_fix) * gamma ; 

    % --- initialise plotting range for band/gap functions ---
    alpha_vals   = linspace(-pi, pi, N);
    beta_vals_1  = linspace(-3 + beta_fixed,  3 + beta_fixed, N);
    beta_vals_2_fine  = linspace(- 2 + beta_fixed, 2 + beta_fixed, search_range);
    
    band_functions = zeros(2, N);

    for k = 1:N
        alpha = alpha_vals(k);
        C = quasiperiodic_capacitance_matrix(alpha, beta_fixed, gamma, ell_fix, spa_fix);
        ev = eig(C); 
        band_functions(:, k) = sqrt( delta * real(ev) );
    end
    
    % --- Mark the limit of the spectrum ---
    C = quasiperiodic_capacitance_matrix(pi, beta_fixed, gamma, ell_fix, spa_fix);
    lower_gaps = real(sqrt( delta * eig(C)));

    C = quasiperiodic_capacitance_matrix(0, beta_fixed, gamma, ell_fix, spa_fix);
    upper_gaps = real(sqrt( delta * eig(C)));

% --- Generate gap bands ---
    % --- Gap band for alpha = 0 ---
    alpha_0 = 0;
    gap_functions_1 = zeros(2, N);

    for k = 1:N
        beta_gap = beta_vals_1(k);
        C = quasiperiodic_capacitance_matrix(alpha_0, beta_gap, gamma, ell_fix, spa_fix);
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
    
    % --- Find range for real valued gap functions ---
    for k = 1:search_range
        beta_gap = beta_vals_2_fine(k);
        C = quasiperiodic_capacitance_matrix(alpha_pi, beta_gap, gamma, ell_fix, spa_fix);
        ev = eig(C);
    
        for j = 1:length(ev)
            lambda = ev(j);
            if abs(imag(lambda)) < 1e-5 && real(lambda) >= -2e-1
                is_real_2(k) = true;
                break;
            end
        end
    end

    onset_indices_2  = find(diff([false, is_real_2]) ==  1);
    offset_indices_2 = find(diff([is_real_2, false]) == -1);
    real_intervals_2 = [beta_vals_2_fine(onset_indices_2); beta_vals_2_fine(offset_indices_2)];
  
    % --- Gap band for alpha = pi ---
    alpha_pi = pi;
    gap_functions_2 = zeros(2, N);
    beta_vals_2  = linspace(real_intervals_2(1, 1), real_intervals_2(2, 1), N);
    
    for k = 1:N
        beta_gap = beta_vals_2(k);
        C = quasiperiodic_capacitance_matrix(alpha_pi, beta_gap, gamma, ell_fix, spa_fix);
        ev = eig(C);
        gap_functions_2(:, k) = NaN(size(gap_functions_2, 1), 1);
    
        for j = 1:min(length(ev), size(gap_functions_2, 1))
            lambda = ev(j);
            if abs(imag(lambda)) < 1e-2 && real(lambda) >= - 1e-2
                gap_functions_2(j, k) = sqrt(delta * real(lambda));
            end
        end
    end
    
    % --- Plot complex band structure an overlay decay rates ---
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

    plot(decay_l,  w, 'bx', 'MarkerSize', 8, 'LineWidth', 2);
    plot(decay_r,  w, 'bx', 'MarkerSize', 8, 'LineWidth', 2);

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

%%

% --- Plot the resonances and defect resonances ---
    figure;
    h1 = plot(resonances, 'ko', 'LineWidth', 4, 'MarkerSize', 8); 
    hold on;
    if nu > 0
        % --- Defect in the top spectral Gap ---
        h2 = plot(n,   def(1), 'rx', 'LineWidth', 4, 'MarkerSize', 8);
        plot(2*n,      def(2), 'rx', 'LineWidth', 4, 'MarkerSize', 8);
    else
        % --- Defect in the lower spectral Gap ---
        h2 = plot(2,   def(1), 'rx', 'LineWidth', 4, 'MarkerSize', 8);
        plot(n+1,      def(2), 'rx', 'LineWidth', 4, 'MarkerSize', 8);
    end

    % --- Mark the quasiperiodic spectrum ---
    yline(double(lim1(1)), 'b--', 'LineWidth', 1);
    yline(double(lim2(1)), 'b--', 'LineWidth', 1);
    yline(double(lim1(2)), 'b--', 'LineWidth', 1);
    yline(double(lim2(2)), 'b--', 'LineWidth', 1);

    xlabel('Index',    'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\omega$', 'Interpreter', 'latex', 'FontSize', fs);
    set(gca, 'FontSize', fs+4, 'TickLabelInterpreter', 'latex');
    legend([h1, h2], {'Band Resonance', 'Defect Resonance'}, 'Interpreter', 'latex', 'Location', 'southeast', 'Box', 'on');
    set(gcf, 'Position', [100, 100, 500, 300]);
    ylim([0, 0.19]);


%% --- Define function ---

function capmat = Capacitance(N, s, gamma,  ell)
    
    capmat = zeros(N, N);  
    
    % --- Popolate the matrix ---
    for i = 1:N
        for j = 1:N
            if i == j
            % --- Populate  diagonal ---
                if i == 1
                    capmat(i,j) = (gamma / s(i)) * (ell(i) / (1 - exp(-gamma * ell(i))));
                % Case 2: 1 < i = j < N
                elseif i > 1 && i < N
                    capmat(i,j) = (gamma / s(i))   * (ell(i) / (1 - exp(-gamma * ell(i)))) ...
                                - (gamma / s(i-1)) * (ell(i) / (1 - exp( gamma * ell(i))));
                % Case 3: i = j = N
                else
                    capmat(i,j) = -(gamma / s(i-1)) * (ell(i) / (1 - exp(gamma * ell(i))));
                end

            % --- Populate lower diagonal ---
                elseif i == j - 1
                    % Case 4: 1 <= i = j - 1 <= N - 1
                    capmat(i,j) = -(gamma / s(i)) * (ell(i) / (1 - exp(-gamma * ell(j))));

            % --- Populate upper diagonal ---
                elseif i == j + 1 
                    % Case 5: 2 <= i = j + 1 <= N
                    capmat(i,j) = (gamma / s(j)) * (ell(i) / (1 - exp(gamma * ell(j))));
            end
        end
    end
end


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
