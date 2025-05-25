%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [May 2025]
    Description:  [Spectral Convergence and Pseudospectrum]
    --------------------------------------------------------------
%}

    clear all;
    close all;

% --- Define fixed parameters ---
    gamma = 1;          % Gauge potential
    delta = 0.001;      % Contrast parameter
    s1    = 0.5;        % Spacing betweeen the resonators
    l1    = 0.5;        % Length of the resonators
    L     = s1 + l1;    % Length of the unit cell
    Nx    = 1000;       % Number of plotting points in the bands
    fs    = 18;         % Fontsize in plot annotation
    lw    = 3;          % Linewidth of spectral bands
    
    % --- Renormalise the lengths ---
    s1 = s1 / L;
    l1 = l1 / L;
    L  = 1;

% --- Generate the band betaTilde ---
    a = (gamma / s1) * (l1 / (1 - exp(-gamma * l1))) - (gamma / s1) * (l1 / (1 - exp(gamma * l1)));
    b =  gamma / s1  *  l1 / (1 - exp( gamma * l1));
    c = -gamma / s1  *  l1 / (1 - exp(-gamma * l1));

    a = delta * a;
    b = delta * b;
    c = delta * c;

    r = 0.5 * gamma * l1;

    lambda = linspace(0, 0.2^2, 1000);
    betaTilde = acosh(-(a-lambda) ./ (2*b* exp(r)));


%% --- Compute Frequency convergence ---

            % System size and eta range
    N        = 12;
    sizes    = 10:2:N+2;
    etas     = linspace(0.1, 5, 10); 
    slopes   = zeros(size(etas));
    w_closed = zeros(size(etas));
    
    % Loop over eta
    for j = 1:length(etas)
        eta = etas(j);
    
        % --- Closed form defect resonance ---
        w_def_closed_form = sqrt(...
            (a * (eta + 1)^2 + eta/abs(eta) * sqrt((eta + 1)^2 * (a^2 * eta^2 + 8 * b * c * eta + 4 * b * c))) ...
            / (2 * eta + 1) );
    
        w_closed(j) = w_def_closed_form;
    
        res = zeros(1, length(sizes));
        idx = 1;
    
        % Loop over system sizes
        for n = sizes
            defect_site = floor(n/2);
            l = l1 * ones(1, n);
            s = s1 * ones(1, n);
            capmat = Capacitance(n, s, gamma, l);
            D = eye(n);
            D(defect_site, defect_site) = 1 + eta;
            Def = D * capmat;
    
            [~, eigenvalues_matrix] = eig(Def);
            eigenvalues = sort(abs(diag(eigenvalues_matrix)));
            resonances  = sqrt(eigenvalues * delta);
    
            % Pick the right resonance based on eta sign
            if eta > 0
                w_def = resonances(n);
            else
                w_def = resonances(2);
            end
            res(idx) = abs(w_def - w_def_closed_form);
            idx = idx + 1;
        end
    
        % Compute convergence slope
        log_res = log(res);
        p = polyfit(sizes, log_res, 1);
        slopes(j) = -p(1);  
    end


    % --- Plot the convergence rate against beta ---

    figure;
    plot(real(betaTilde), sqrt(lambda),'r-', 'LineWidth', 2);
    hold on;
    plot(slopes, w_closed, 'bx', 'MarkerSize', 8, 'LineWidth', 2);
    
    legend({'$\tilde{\beta}(\lambda)$', '$B(\lambda)$'}, 'Interpreter', 'latex', 'Location', 'best');

    %xlabel('$\tilde{\beta}$ and $B$ respectively', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\lambda$' , 'Interpreter', 'latex', 'FontSize', fs);

    ylim([0.07, 0.16]);

    set(gca, 'FontSize', fs+4, 'TickLabelInterpreter', 'latex');
    set(gcf, 'Position', [100, 100, 400, 250]); 
    grid on;


%% --- Truncate Laurent Operator convergence ---

    etas       = linspace(0.01, 5, 10);  
    N          = 107;
    n_values   = 20:4:N-60;
    
    
    % Preallocate storage
    decay_slopes = zeros(size(etas));
    lambda_ns    = zeros(size(etas));
    
    for j = 1:length(etas)
        eta = etas(j);
        pos = 3;
        defect_pos = floor(N/pos);
    
        % Large matrix to approximate "Operator"
        l = l1 * ones(1, N);
        s = s1 * ones(1, N);
        C = Capacitance(N, s, gamma, l);
        D = eye(N);
        D(defect_pos, defect_pos) = 1 + eta;
        C = C * D;
    
        [v, ~] = eigs(C, 1, 'largestabs'); 
    
        residuals = zeros(size(n_values));
        lambda_n_vals = zeros(size(n_values));
    
        for i = 1:length(n_values)
            n = n_values(i);
            l_n = l1 * ones(1, n);
            s_n = s1 * ones(1, n);
    
            Cn = Capacitance(n, s_n, gamma, l_n);
            Dn = eye(n);
            Dn(floor(n/pos), floor(n/pos)) = 1 + eta;
            Cn = Cn * Dn;
    
            lambda_n = eigs(Cn, 1, 'largestabs');
            lambda_n_vals(i) = lambda_n;
    
            % Truncate eigenvector from operator
            start_idx = floor((N - n)/pos) + 1;
            end_idx   = start_idx + n - 1;
            v_trunc   = v(start_idx:end_idx);
    
            % Residual norm
            residual_vector = (Cn - lambda_n * eye(n)) * v_trunc;
            residuals(i)    = norm(residual_vector);
        end
    
        % Measure exponential decay rate
        p = polyfit(n_values, log(residuals), 1);
        decay_slopes(j) = p(1); 
        lambda_ns(j)    = lambda_n_vals(end); 
    end
    
    figure;
    plot(r-real(betaTilde), sqrt(lambda),'r-', 'LineWidth', 2);
    hold on;
    plot( decay_slopes*pos, sqrt(delta *lambda_ns), 'bx', 'MarkerSize', 8, 'LineWidth', 2);

    legend({'$r-\tilde{\beta}(\lambda)$', '$B(\lambda)$'}, 'Interpreter', 'latex', 'Location', 'best');

    %xlabel('$\tilde{\beta}(\lambda) - r$ and $B$ respectively', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\lambda$' , 'Interpreter', 'latex', 'FontSize', fs);

    ylim([0.07, 0.16]);

    xticks([-2, -1, 0, r]); 
    xticklabels({'$-2$', '$-1$', '$0$', '$r$'});
    set(gca, 'FontSize', fs+4, 'TickLabelInterpreter', 'latex');
    set(gcf, 'Position', [100, 100, 400, 250]); 
    grid on;


%% --- Truncate Toeplitz Operator convergence ---

    % Parameters and constants
    etas       = linspace(0.1, 5, 10); 
    N          = 81;
    n_values   = 30:4:N-40;
    
    
    % Preallocate storage
    decay_slopes = zeros(size(etas));
    lambda_ns    = zeros(size(etas));
    
    for j = 1:length(etas)
        eta = etas(j);
        pos = 2;
        defect_pos = floor(N/pos);
    
        % Large matrix to approximate "Operator"
        l = l1 * ones(1, N);
        s = s1 * ones(1, N);
        C = Capacitance(N, s, gamma, l);
        D = eye(N);
        D(defect_pos, defect_pos) = 1 + eta;
        C = C * D;
    
        [v, ~] = eigs(C, 1, 'largestabs'); 
    
        residuals = zeros(size(n_values));
        lambda_n_vals = zeros(size(n_values));
    
        for i = 1:length(n_values)
            n = n_values(i);
            l_n = l1 * ones(1, n);
            s_n = s1 * ones(1, n);
    
            Cn = Capacitance(n, s_n, gamma, l_n);
            Dn = eye(n);
            Dn(floor(n/pos), floor(n/pos)) = 1 + eta;
            Cn = Cn * Dn;
    
            lambda_n = eigs(Cn, 1, 'largestabs');
            lambda_n_vals(i) = lambda_n;
    
            % Truncate eigenvector from large matrix
            v_trunc   = v(1:n);
    
            % Residual norm
            residual_vector = (Cn - lambda_n * eye(n)) * v_trunc;
            residuals(i)    = norm(residual_vector) ;
        end
    
        % Measure exponential decay rate
        p = polyfit(n_values, log(residuals), 1);
        decay_slopes(j) = p(1); 
        lambda_ns(j)    = lambda_n_vals(end); 
    end

    figure;
    plot(real(betaTilde)-r, sqrt(lambda),'r-', 'LineWidth', 2);
    hold on;
    plot( decay_slopes, sqrt(delta *lambda_ns), 'bx', 'MarkerSize', 8, 'LineWidth', 2);
    legend({'$\tilde{\beta}(\lambda)-r$', '$B(\lambda)$'}, 'Interpreter', 'latex', 'Location', 'best');

    %xlabel('$\tilde{\beta}(\lambda) - r$ and $B$ respectively', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\lambda$' , 'Interpreter', 'latex', 'FontSize', fs);

    ylim([0.07, 0.16]);

    xticks([-r, 0, 1, 2]); 
    xticklabels({'$-r$', '$0$', '$1$', '$2$'});
    set(gca, 'FontSize', fs+4, 'TickLabelInterpreter', 'latex');
    set(gcf, 'Position', [100, 100, 400, 250]); 
    grid on;

%% --- Define the Capacitance Matrix ---

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

