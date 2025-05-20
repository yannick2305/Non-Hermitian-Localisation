%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [March 2025]
    Description:  [Finite Defected system]
    --------------------------------------------------------------
%}

    close all;
    clear all;

% --- Parameters ---
    n     = 300;                % Number of unit cells in the chain
    l1    = 0.5;                % Length of a monomer resonator
    s1    = 0.5;                % Spacing of a monomer resonator
    gamma = 3.2422;             % Gauge potential
    eta   = 1.9;                % Defect parameter (change in wavespeed)
    delta = 0.001;              % Contrast
    fs = 18;                    % Fontsize in the plot
    lw = 3;                     % Linewidth of the modes plot

    % --- Renormalise the lengths ---
    L  = l1 + s1;
    s1 = s1 / L;
    l1 = l1 / L;
    L  = 1;

% --- Find optimal defect size for cloaking ---
    % --- Compute the spectral bands ---
    a = (gamma / s1) * (l1 / (1 - exp(-gamma * l1))) - (gamma / s1) * (l1 / (1 - exp(gamma * l1)));
    b =  gamma / s1  *  l1 / (1 - exp( gamma * l1));
    c = -gamma / s1  *  l1 / (1 - exp(-gamma * l1));
    
    a = delta * a;
    b = delta * b;
    c = delta * c;
    
    % --- Find limit of the winding region ---
    alpha_fixed_2 = pi;
    w_beta_0 = (gamma * l1)/s1 * ( (1 - exp(-1i*(alpha_fixed_2 + 1i* 0)*L)) / (1 - exp(-gamma * l1)) + (exp(1i*(alpha_fixed_2 + 1i* 0)*L)-1) / (1-exp(gamma*l1)) );
    w_beta_0 = sqrt(delta * abs(w_beta_0));

    func = @(eta) sqrt((a * (eta + 1)^2 + sign(eta) * sqrt((eta + 1)^2 * (a^2 * eta^2 + 8 * b * c * eta + 4 * b * c))) / (2 * eta + 1)) - w_beta_0;
    
    eta0 = 1;  % Initial Guess
    options = optimoptions('fsolve', 'Display', 'off');
    eta_solution = fsolve(func, eta0, options);
    
    % Display result
    disp('----------------------------------------');
    fprintf('Cloaking eta:     %.6f\n', eta_solution);

    eta = eta_solution;

% --- Compute the resonances of the defected monomer chain ---

    % --- Initialise resonator chain ---
    l = l1 * ones(1, n);
    s = s1 * ones(1, n);
        
    % --- Generate Capacitance matrix of pristine structure ---
    capmat = Capacitance(n, s, gamma, l);
          
    % --- Define defect Matrix ---
    D = eye(n);
    defect_site = floor(n/2);
    D(defect_site, defect_site) = 1 + eta;

    Def = D * capmat;

    %% --- Truncate eigenvector as Pseudoeigenvector estimate ---

    n_values = 10:4:300;  % Truncation sizes
    first_entries = zeros(size(n_values));  
    
    for idx = 1:length(n_values)
        
        n = n_values(idx);
        l = l1 * ones(1, n);
        s = s1 * ones(1, n);
    
        % --- Generate Capacitance matrix of pristine structure ---
        capmat = Capacitance(n, s, gamma, l);
    
        % --- Define defect Matrix ---
        D = eye(n);
        defect_site = floor(n / 2);
        D(defect_site, defect_site) = 1 + eta;
    
        Def = D * capmat;
    
        % --- Compute eigenvalues and eigenvectors ---
        [V, E] = eig(Def);
        [~, max_idx] = max(abs(diag(E)));  % Defect eigenfrequency
        
        v_max = V(:, max_idx);             % Defect eigenvector
        v_max = v_max / norm(v_max);       % Normalize eigenvector 
        
        first_entries(idx) = v_max(1);     % Store first entry
    end
    
    % --- Plot the result ---
    figure;
    loglog(n_values, abs(real(first_entries)), '-o');
    xlabel('n (matrix size)');
    ylabel('First entry of eigenvector with max eigenfrequency');
    title('First Entry vs. Matrix Size');
    grid on;

    % --- Determine algebraic decay rate ---
    x = n_values(:);                        
    y = abs(real(first_entries(:)));         
    logx = log10(x);
    logy = log10(y);
    
    p = polyfit(logx, logy, 1);
    slope = p(1);
    
    fprintf('Algebraic decay rate: %.4f\n', slope);

%% --- Convergence of the defect eigenfrequency ---

    N   = 32;   % Range of the system sizes
    eta = 0.5;  % Defect size 
    % --- Cosed form of defect resonance in infinite chain ---
    w_def_closed_form = sqrt((a * (eta + 1)^2 + eta/abs(eta) * sqrt((eta + 1)^2 * (a^2 * eta^2 + 8 * b * c * eta + 4 * b * c))) / (2 * eta + 1)); 
    fprintf('--------------------------------------\n');
    disp(['Closed Sol:        ',  num2str(w_def_closed_form)]);
    fprintf('Cloaking eta:      %.6f\n', eta);
    
    non_reciprocity = 0.5 * gamma * l1;
    fprintf('Non-reciprocity:   %.6f\n', non_reciprocity);

    % --- Iterate over chain sizes ---
    sizes = 10:2:N+2; 
    res   = ones(1, length(sizes));
    idx   = 1;
    
    for n = sizes
        % --- Generate defected capacitance matrix ---
        defect_site = floor(n/2);
        l = l1 * ones(1, n);
        s = s1 * ones(1, n);
        capmat = Capacitance(n, s, gamma, l);
        D = eye(n);
        D(defect_site, defect_site) = 1 + eta;
        Def = D * capmat;
        
        % --- Compute the resonances ---
        [~, eigenvalues_matrix] = eig(Def);
        eigenvalues = sort(diag(eigenvalues_matrix));
        eigenvalues = abs(eigenvalues);
        resonances  = sqrt(eigenvalues * delta);
        if eta > 0
            w_def = resonances(n);
        else
            w_def = resonances(2);
        end
        res(idx) = abs(w_def - w_def_closed_form);
        idx = idx + 1;
        
        [V, Dia] = eig(Def);
        eigenvalues = diag(Dia);
        [~, sortIdx] = sort(eigenvalues);
        sortedEigenvectors = V(:, sortIdx);
        %figure;
        %plot(1:n, abs(sortedEigenvectors(:, n)), '-', 'Color', 0.65 * [1, 1, 1], 'LineWidth', lw * 1.5, 'MarkerSize', 15);
    end
    
    % --- Fit a regression line on the semilog plot ---
    log_res = log(res(1:end)); 
    p = polyfit(sizes, log_res, 1); 
    slope = p(1); 
    
    % --- Plot results ---
    figure;
    semilogy(sizes, res(1:end), 'kx', 'LineWidth', 2, 'MarkerSize', 9);  
    hold on;
    
    xlabel('N',                                       'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\bigl|\omega_0 - \omega_0^{(N)}\bigr|$', 'Interpreter', 'latex', 'FontSize', fs);
    set(gca, 'FontSize', fs+4, 'TickLabelInterpreter', 'latex');
    set(gcf, 'Position', [100, 100, 400, 400]); 
    grid on;
    
    disp(['Exponential convergence rate: ', num2str(slope)]);

%% --- Truncated eigenvector as Pseudoeigenvector ---

    N   = 201;
    def = 0.1;              % Diff from tunneling defect size 
    l   = l1 * ones(1, N);
    s   = s1 * ones(1, N);
            
    % Compute large matrix to approximate Toeplitz operator
    C = Capacitance(N, s, gamma, l);
    D = eye(N);
    D(floor(N/2), floor(N/2)) = 1 + eta + def;
    C = C * D;
    
    [v, lambda] = eigs(C, 1, 'largestabs');
    
    n_values  = 11:4:N-2;
    residuals = zeros(size(n_values));
    
    % Loop over n and compute residuals
    for i = 1:length(n_values)
        n = n_values(i);
    
        Cn = Capacitance(n, s, gamma, l);
        D = eye(n);
        D(floor(n/2), floor(n/2)) = 1 + eta + def;
        Cn = Cn * D;
    
        lambda_n = eigs(Cn, 1, 'largestabs');
        
        % Truncate the vector of operator around center
        start_idx = floor((N - n)/2) + 1;
        end_idx   = start_idx + n - 1;
        v_trunc   = v(start_idx:end_idx);  
        
        % Compute residual
        residual_vector = (Cn - lambda_n * eye(n)) * v_trunc;
        residuals(i)    = norm(residual_vector) / norm(v_trunc);
    end

    % --- Plot the residuals ---
    figure;
    semilogy(n_values, residuals, '-o', 'LineWidth', 2);
    %loglog(n_values, residuals, '-o', 'LineWidth', 2);
    xlabel('n');
    ylabel('Residual Norm');
    grid on;


%% --- Function generating the finite Capacitance matrix

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
