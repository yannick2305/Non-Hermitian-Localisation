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
    n     = 30;                 % Number of unit cells in the chain
    l1    = 0.5;                % Length of a monomer resonator
    s1    = 0.5;                % Spacing of a monomer resonator
    gamma = 3;                  % Gauge potential
    eta   = 1.9;                % Defect parameter (change in wavespeed)
    delta = 0.001;              % Contrast
    defect_site = floor(n/2);   % Site of the defect
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

% --- Compute the resonances of the defected monomer chain ---

    % --- Initialise resonator chain ---
    l = l1 * ones(1, n);
    s = s1 * ones(1, n);
        
    % --- Generate Capacitance matrix of pristine structure ---
    capmat = Capacitance(n, s, gamma, l);
          
    % --- Define defect Matrix ---
    D = eye(n);
    D(defect_site, defect_site) = 1 + eta;

    Def = D * capmat;
       
    % --- Compute the resonances ---
    [Vect, eigenvalues_matrix] = eig(Def);
    eigenvalues = sort(diag(eigenvalues_matrix));
    eigenvalues = abs(eigenvalues);
    resonances  = sqrt(eigenvalues * delta);

    % --- Limit of upper gap ---
    Upper_gap = (gamma * l1)/s1 * ( (1 - exp(-1i*(pi + 1i* (-gamma * l1/2))*L)) / (1 - exp(-gamma * l1)) + (exp(1i*(pi + 1i* (-gamma * l1/2))*L)-1) / (1-exp(gamma*l1)) ); 
    Upper_gap = sqrt(delta * abs(Upper_gap));
    
    % --- Limit of lower gap ---
    Lower_gap = (gamma * l1)/s1 * ( (1 - exp(-1i*(0 + 1i* (-gamma * l1/2))*L)) / (1 - exp(-gamma * l1)) + (exp(1i*(0 + 1i* (-gamma * l1/2))*L)-1) / (1-exp(gamma*l1)) );
    Lower_gap = sqrt(delta * Lower_gap);

    % ---Find defect frequency based on the type of defect ---
    if eta > 0
        % --- Defect in the top spectral Gap ---
        w_def = resonances(n);
        resonances(n) = NaN;

    else
        % --- Defect in the lower spectral Gap ---
        w_def = resonances(2);
        resonances(2) = NaN;
    end   

    % --- Display the defect eigenfrequency ---
    disp(['Defect Frequency: ', num2str(w_def)]);
    
% --- Illustrate the result ---
    figure;
    h1 = plot(resonances, 'ko', 'LineWidth', 4, 'MarkerSize', 8); 
    hold on;

    if eta > 0
        % --- Defect in the top spectral Gap ---
        h2 = plot(n, w_def, 'rx', 'LineWidth', 4, 'MarkerSize', 8);
    else
        % --- Defect in the lower spectral Gap ---
        h2 = plot(2, w_def, 'rx', 'LineWidth', 4, 'MarkerSize', 8);
    end
     
    yline(Lower_gap, 'b--', 'LineWidth', 1);
    yline(Upper_gap, 'b--', 'LineWidth', 1);
    xlabel('Index',    'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\omega$', 'Interpreter', 'latex', 'FontSize', fs);
    set(gca, 'FontSize', fs+4, 'TickLabelInterpreter', 'latex');
    legend([h1, h2], {'Band Resonance', 'Defect Resonance'}, 'Interpreter', 'latex', 'Location', 'southeast', 'Box', 'on');
    grid off;
    ylim([0, 0.159]);
    set(gcf, 'Position', [100, 100, 500, 300]); 


%% --- Plot the eigenmodes ---

    % --- Compute the eigenmodes ---
    [V, Dia] = eig(Def);
    eigenvalues = diag(Dia);
    [~, sortIdx] = sort(eigenvalues);
    sortedEigenvectors = V(:, sortIdx); 

    figure;
    hold on; 

    if eta > 0
        % --- Defect in the top spectral Gap ---
        plot(1:n, abs(sortedEigenvectors(:, n)), '-', 'Color', 0.65 * [1, 1, 1], 'LineWidth', lw * 1.5, 'MarkerSize', 15);
    else
        % --- Defect in the lower spectral Gap ---
        plot(1:n, abs(sortedEigenvectors(:, 2)), '-', 'Color', 0.65 * [1, 1, 1], 'LineWidth', lw * 1.5, 'MarkerSize', 15);
    end

    xlabel('Eigenvector index', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$|u(x)|$',          'Interpreter', 'latex', 'FontSize', fs);
    grid on;
    set(gca, 'YScale', 'log');
    set(gca, 'FontSize', fs+4, 'TickLabelInterpreter', 'latex');
    set(gcf, 'Position', [100, 100, 400, 400]); 
    hold off;


%% --- Convergence to the infinite system ---

    N = 52; % Range of the system sizes

    % --- Cosed form of defect resonance in infinite chain ---
    w_def_closed_form = sqrt((a * (eta + 1)^2 + eta/abs(eta) * sqrt((eta + 1)^2 * (a^2 * eta^2 + 8 * b * c * eta + 4 * b * c))) / (2 * eta + 1)); 
    disp(['Closed Sol:       ', num2str(w_def_closed_form)]);
    
    % --- Iterate over chain sizes ---
    sizes = 3:2:N+2; % System sizes
    res = ones(1, floor(N/2) );
    idx = 1;
    
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
    end
    
    % --- Fit a regression line on the semilog plot ---
    log_res = log(res(1:end)); 
    p = polyfit(sizes, log_res, 1); 
    slope = p(1); 
    
    % --- Plot results ---
    figure;
    semilogy(sizes, res(1:end), 'kx', 'LineWidth', 2, 'MarkerSize', 9);  
    hold on;
    
    % --- Plot regression line ---
    %semilogy(sizes, exp(polyval(p, sizes)), 'r--', 'LineWidth', 2);
    
    xlabel('N', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\bigl|\omega_0 - \omega_0^{(N)}\bigr|$', 'Interpreter', 'latex', 'FontSize', fs);
    set(gca, 'FontSize', fs+4, 'TickLabelInterpreter', 'latex');
    set(gcf, 'Position', [100, 100, 400, 400]); 
    grid on;
    
    disp(['Exponential convergence rate: ', num2str(slope)]);



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
