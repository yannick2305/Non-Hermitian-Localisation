%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [March 2025]
    Description:  [Cloaking in monomer chain]
    --------------------------------------------------------------
%}

close all;
clear all;

% --- Parameters ---
    n  = 61; % Use odd number  % Number of unit cells in the chain
    l1 = 0.5;                   % Length of a monomer resonator
    s1 = 0.5;                   % Spacing of a monomer resonator
    gamma = 3;                  % Gauge potential
    defect_site_1 = floor(n/3); % Site of the defect
    defect_site_2 = n-floor(n/3);
    delta = 0.001;              % Contrast
    fs = 18;                    % Fontsize in the plot
    lw = 3;                     % Linewidth of the modes plot
    
    % --- Renormalise the lengths ---
    L = l1 + s1;
    s1 = s1 / L;
    l1 = l1 / L;
    L  = 1;

% --- Compute the resonances of the chain ---
    % --- Initialise resonator chain ---
    l = l1 * ones(1, n);
    s = s1 * ones(1, n);

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

    % Define the function whose root we want to find
    func = @(eta) sqrt((a * (eta + 1)^2 + sign(eta) * sqrt((eta + 1)^2 * (a^2 * eta^2 + 8 * b * c * eta + 4 * b * c))) / (2 * eta + 1)) - w_beta_0;
    
    eta0 = 1;  % Initial Guess
    options = optimoptions('fsolve', 'Display', 'off');
    eta_solution = fsolve(func, eta0, options);
    
    % Display result
    disp('----------------------------------------');
    fprintf('Computed eta:     %.6f\n', eta_solution);

    nu = eta_solution;
        
% --- Generate Capacitance matrix of pristine structure ---
    capmat = Capacitance(n, s, gamma, l);
          
    % --- Define defect Matrix ---
    D = eye(n);
    D(defect_site_1, defect_site_1) = 1 + nu;
    D(defect_site_2, defect_site_2) = 1 + nu;

    Def = D * capmat;
    

%% --- Plot the eigenmodes ---

    % --- Compute the eigenmodes ---
    [V, Dia] = eig(Def);
    eigenvalues = diag(Dia);
    [~, sortIdx] = sort(eigenvalues);
    sortedEigenvectors = V(:, sortIdx); 

    % --- Plot the eigenvector entries ---
    figure;

    plot(1:n, abs(sortedEigenvectors(:, n)), '-', 'Color', 0.65 * [1, 1, 1], 'LineWidth', lw * 1.5, 'MarkerSize', 15);

    xlabel('Eigenvector index', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$|u(x)|$',          'Interpreter', 'latex', 'FontSize', fs);
    grid on;
    set(gca, 'YScale', 'log');
    set(gca, 'FontSize', fs+4, 'TickLabelInterpreter', 'latex');
    set(gcf, 'Position', [100, 100, 400, 400]); 
    hold off;


%% --- Plot the resonances ---

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
    
    % --- Defect in the top spectral Gap ---
    w_def = zeros(1,2);
    w_def(1) = resonances(n);
    w_def(2) = resonances(n-1);
    resonances(n)   = NaN;
    resonances(n-1) = NaN;

    % --- Display the defect eigenfrequency ---
    disp(['Defect Frequency: ', num2str(w_def(1))]);

    figure;
    h1 = plot(resonances, 'ko', 'LineWidth', 4, 'MarkerSize', 8); 
    hold on;

    % --- Defect in the top spectral Gap ---
    h2 = plot(n-1:n, w_def, 'rx', 'LineWidth', 4, 'MarkerSize', 8);

    yline(Lower_gap, 'b--', 'LineWidth', 1);
    yline(Upper_gap, 'b--', 'LineWidth', 1);
    xlabel('Index',    'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\omega$', 'Interpreter', 'latex', 'FontSize', fs);
    set(gca, 'FontSize', fs+4, 'TickLabelInterpreter', 'latex');
    legend([h1, h2], {'Band Resonance', 'Defect Resonance'}, ...
    'Interpreter', 'latex', ...
    'Location', 'southeast', ...
    'Box', 'on');
    grid off;
    ylim([0, w_def(1) * 1.2]);
    set(gcf, 'Position', [100, 100, 500, 300]); 


%% --- Function generating the finite Capacitance matrix

function capmat = Capacitance(N, s, gamma_va,  ell)

    gamma = zeros(1,N);
    
    for i = 1:N
        if i <= floor(N/2)
            gamma(i) = gamma_va;
        else
            gamma(i) = -gamma_va;
        end
    end

    
    capmat = zeros(N, N);  
    
    % --- Popolate the matrix ---
    for i = 1:N
        for j = 1:N
            if i == j
            % --- Populate  diagonal ---
                % Case 1: 1 = i = j
                if i == 1
                    capmat(i,j) = (gamma(i) / s(i)) * (ell(i) / (1 - exp(-gamma(i) * ell(i))));
                % Case 2: 1 < i = j < N
                elseif i > 1 && i < N
                    capmat(i,j) = (gamma(i) / s(i))   * (ell(i) / (1 - exp(-gamma(i) * ell(i)))) ...
                                - (gamma(i) / s(i-1)) * (ell(i) / (1 - exp( gamma(i) * ell(i))));
                % Case 3: i = j = N
                else
                    capmat(i,j) = -(gamma(i) / s(i-1)) * (ell(i) / (1 - exp(gamma(i) * ell(i))));
                end

            % --- Populate lower diagonal ---
                elseif i == j - 1
                    % Case 4: 1 <= i = j - 1 <= N - 1
                    capmat(i,j) = -(gamma(j) / s(i)) * (ell(i) / (1 - exp(-gamma(j) * ell(j))));

            % --- Populate upper diagonal ---
                elseif i == j + 1 
                    % Case 5: 2 <= i = j + 1 <= N
                    capmat(i,j) = (gamma(j) / s(j)) * (ell(i) / (1 - exp(gamma(j) * ell(j))));
            end
        end
    end
end

