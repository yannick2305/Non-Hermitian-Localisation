%{
    -----------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [Mai 2025]
    Description:  [Spectral plot Defected non-Hermitian Hamiltonian]
    -----------------------------------------------------------------
%}

clear all;
close all;

% --- Parameters ---
    n           = 29;           % Dimension of Hamiltonian
    v           = 1;            % Potential
    gamma       = 0.4;          % Gauge Potential
    defect_site = floor(n/2);   % Position of the defect
    
    fs = 16;   
    lw = 1.5;

% --- Compute the decay lengths for the defects ---
    % --- Set defect range ---
    d_range = linspace(-4, 4, 6);

    % --- Manually add defects at edge of winding region ---
    d_range = [d_range, 0, -4*sinh(gamma)];
    
    % --- Initialise arrays ---
    frequencies           = zeros(size(d_range));
    left_decay_lengths    = zeros(size(d_range));
    right_decay_lengths   = zeros(size(d_range));
    frequencies_2         = zeros(size(d_range));
    left_decay_lengths_2  = zeros(size(d_range));
    right_decay_lengths_2 = zeros(size(d_range));
    
    for idx = 1:length(d_range)
        d = -(2*sinh(gamma) + d_range(idx));
 
        % --- Generate Toeplitz matrix ---
        T = generate_tridiagonal_toeplitz(n, v, gamma);
        T(defect_site, defect_site) = 1 + d;
        
        % --- Eigen decomposition ---
        [V, D] = eig(T);
        eigenvalues = diag(D);
        [sorted_eigenvalues, ind] = sort(real(eigenvalues), 'descend');
        sorted_vectors = V(:, ind);
        
        % --- For positive defect ---
        if d_range(idx) < 0
            % --- Select defect eigenvalue ---
            max_vec = sorted_vectors(:, 1);
            max_val = sorted_eigenvalues(1);

            % --- Split left and right parts ---
            left_part  = abs(max_vec(1:defect_site-1));
            right_part = abs(max_vec(defect_site+1:end));
            
            % --- Fit exponential decay --- 
            x_left    = (defect_site-1:-1:1)';
            x_right   = (1:length(right_part))';
            left_log  = log(left_part);
            right_log = log(right_part);
            
            % --- Fit straight line to log data ---
            p_left  = polyfit(x_left, left_log, 1);
            p_right = polyfit(x_right, right_log, 1);

            lambda_left  = p_left(1);
            lambda_right = p_right(1);
        
            % --- Save the result ---
            frequencies(idx)         = max_val;
            left_decay_lengths(idx)  = -lambda_left;
            right_decay_lengths(idx) =  lambda_right;

            invalid_idx = frequencies == 0;
            frequencies(invalid_idx)         = [];
            left_decay_lengths(invalid_idx)  = [];
            right_decay_lengths(invalid_idx) = [];

        % --- For negative defect ---
        else

            % --- For negative defect ---
            min_vec = sorted_vectors(:, n);
            min_val = sorted_eigenvalues(n);

            % --- Split left and right parts ---
            left_part_2 = abs(min_vec(1:defect_site-1));
            right_part_2 = abs(min_vec(defect_site+1:end));
            
            % --- Fit exponential decay --- 
            x_left_2    = (defect_site-1:-1:1)';
            x_right_2   = (1:length(right_part_2))';
            left_log_2  = log(left_part_2);
            right_log_2 = log(right_part_2);
            
            % --- Fit straight line to log data ---
            p_left_2  = polyfit(x_left_2, left_log_2, 1);
            p_right_2 = polyfit(x_right_2, right_log_2, 1);
      
            lambda_left_2  = p_left_2(1);
            lambda_right_2 = p_right_2(1);
           
            % --- Save the result ---
            frequencies_2(idx) = min_val;
            left_decay_lengths_2(idx)  = -lambda_left_2;
            right_decay_lengths_2(idx) =  lambda_right_2;

            invalid_idx_2 = frequencies_2 == 0;
            frequencies_2(invalid_idx_2)         = [];
            left_decay_lengths_2(invalid_idx_2)  = [];
            right_decay_lengths_2(invalid_idx_2) = [];

        end
    end

% --- Generate the spectral bands ---

    beta_fix     = - gamma;         
    alpha_fix_0  = 0;
    alpha_fix_pi = pi;
    
    % --- Define the range of the bands ---
    alpha = linspace(-pi, pi, 1000);  
    beta  = linspace(-5,   5, 1000);
    
    % --- Compute the band functions ---
    f_alpha    = v - 2 * cos(alpha) .* cosh(gamma + beta_fix);
    f_beta_pi  = v - 2 * cos(alpha_fix_pi) .* cosh(gamma + beta);
    f_beta_0   = v - 2 * cos(alpha_fix_0) .* cosh(gamma + beta);
    
    % --- Compute the limit of the spectrum ---
    Upper_gap = v - 2 * cos(pi) .* cosh(gamma + beta_fix);
    Lower_gap = v - 2 * cos(0) .* cosh(gamma + beta_fix);

    % --- Compute the limit of the winding region ---
    Limit_winding_top    = v - 2 * cos(alpha_fix_pi) .* cosh(gamma + 0); 
    Limit_winding_bottom = v - 2 * cos(alpha_fix_0) .* cosh(gamma + 0); 

    % --- Plot the Band functions ---
    figure;

    xlimit = 7;
    x_coords   = [-xlimit xlimit xlimit -xlimit];
    y_coords   = [Limit_winding_bottom Limit_winding_bottom Limit_winding_top Limit_winding_top];
    y_coords_1 = [ Limit_winding_top Limit_winding_top 10 10];
    y_coords_2 = [ Upper_gap Upper_gap Lower_gap Lower_gap];
    y_coords_3 = [ -4 -4 Limit_winding_bottom Limit_winding_bottom];

    fill(x_coords, y_coords, [0.95 0.95 1]);      % lighter blue
    hold on;
    fill(x_coords, y_coords_1, [0.9 0.95 0.85] ); % lighter green
    fill(x_coords, y_coords_3, [0.9 0.95 0.85] );
    fill(x_coords, y_coords_2, [0.95 0.8 0.9]);   % lighter pink/purple

    plot(alpha, f_alpha, 'k', 'LineWidth', 3);
    
    plot(beta,            f_beta_pi,                 'r-', 'LineWidth', 3);
    plot(beta,            f_beta_0,                  'r-', 'LineWidth', 3);
    plot([0 0],           [Lower_gap-3 Lower_gap],   'k-', 'LineWidth', 1); 
    plot([pi pi],         [Upper_gap 3 * Upper_gap], 'k-', 'LineWidth', 1);
    plot([-pi -pi],       [Upper_gap 3 * Upper_gap], 'k-', 'LineWidth', 1);
    plot([-gamma -gamma], [Lower_gap  Upper_gap],    'r-', 'LineWidth', 1);

    % --- Add the decay lengths ---
    plot(left_decay_lengths(1:end-1),    frequencies(1:end-1),   'bx', 'LineWidth', 4, 'MarkerSize', 8);
    plot(right_decay_lengths(1:end-1),   frequencies(1:end-1),   'bx', 'LineWidth', 4, 'MarkerSize', 8);
    
    plot(left_decay_lengths_2(1:end-1),  frequencies_2(1:end-1), 'bx', 'LineWidth', 4, 'MarkerSize', 8);
    plot(right_decay_lengths_2(1:end-1), frequencies_2(1:end-1), 'bx', 'LineWidth', 4, 'MarkerSize', 8);

    % --- Plot the frequencies at limit of winding region ---
    plot(left_decay_lengths(end),    frequencies(end),   'co', 'LineWidth', 4, 'MarkerSize', 8);
    plot(right_decay_lengths(end),   frequencies(end),   'co', 'LineWidth', 4, 'MarkerSize', 8);

    plot(left_decay_lengths_2(end),  frequencies_2(end), 'co', 'LineWidth', 4, 'MarkerSize', 8);
    plot(right_decay_lengths_2(end), frequencies_2(end), 'co', 'LineWidth', 4, 'MarkerSize', 8);

    grid off;

    % --- LaTeX labels and ticks ---
    xlabel('$\alpha$ and $\beta$ respectively',  'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\lambda^{\alpha, \beta, \gamma}$' , 'Interpreter', 'latex', 'FontSize', fs);
    xlim([-3.5, 3.5]);
    ylim([-3, 5]);
    
    % --- Set LaTeX-formatted ticks ---
    xticks([-pi, 0, pi]); 
    xticklabels({'$-\pi$', '$0$', '$\pi$'});
    set(gca, 'FontSize', fs+4, 'TickLabelInterpreter', 'latex');
    set(gcf, 'Position', [100, 100, 500, 400]);  


%% --- Define the functions ---

function T = generate_tridiagonal_toeplitz(n, v, gamma)
    % Create the diagonals
    main_diag = v * ones(n, 1);             % v on the main diagonal
    below_diag = -exp(-gamma) * ones(n-1, 1); % -e^{-gamma} on the subdiagonal
    above_diag = -exp(gamma) * ones(n-1, 1);  % -e^{gamma} on the superdiagonal

    % Create the tridiagonal Toeplitz matrix
    T = diag(main_diag) + diag(below_diag, -1) + diag(above_diag, 1);
end
