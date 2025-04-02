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
    n  = 30;                    % Number of unit cells in the chain
    l1 = 0.5;                   % Length of a monomer resonator
    s1 = 0.5;                   % Spacing of a monomer resonator
    gamma = 1;                  % Gauge potential
    nu = 1.5;                   % Defect parameter (change in wavespeed)
    defect_site = floor(n/2);   % Site of the defect
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
        
    % --- Generate Capacitance matrix of pristine structure ---
    capmat = Capacitance(n, s, gamma, l);
          
    % --- Define defect Matrix ---
    D = eye(n);
    D(defect_site, defect_site) = 1 + nu;
    Def = D * capmat;
       
    % --- Compute the resonances ---
    [~, eigenvalues_matrix] = eig(Def);
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
    if nu > 0
        % --- Defect in the top spectral Gap ---
        w_def = resonances(n);
    else
        % --- Defect in the lower spectral Gap ---
        w_def = resonances(2);
    end   

    % --- Display the defect eigenfrequency ---
    disp('----------------------------------------');
    disp(['Defect Frequency: ', num2str(w_def)]);
    
% --- Illustrate the result ---
    figure;
    plot(resonances, 'ko', 'LineWidth', 4, 'MarkerSize', 8); 
    hold on;
    
    if nu > 0
        % --- Defect in the top spectral Gap ---
        plot(n, resonances(n), 'ro', 'LineWidth', 4, 'MarkerSize', 8);
    else
        % --- Defect in the lower spectral Gap ---
        plot(2, resonances(2), 'ro', 'LineWidth', 4, 'MarkerSize', 8);
    end
     
    yline(Lower_gap, 'b--', 'LineWidth', 1);
    yline(Upper_gap, 'b--', 'LineWidth', 1);
    xlabel('Index',    'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\omega$', 'Interpreter', 'latex', 'FontSize', fs);
    set(gca, 'FontSize', fs+4, 'TickLabelInterpreter', 'latex');
    grid off;
    set(gcf, 'Position', [100, 100, 500, 500]); 


%% --- Plot the eigenmodes ---

    % --- Compute the eigenmodes ---
    [V, Dia] = eig(Def);
    eigenvalues = diag(Dia);
    [~, sortIdx] = sort(eigenvalues);
    sortedEigenvectors = V(:, sortIdx); 

    figure;
    hold on; 

    if nu > 0
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
    hold off;


%% --- Convergence to the infinite system ---

    N = 10; % Range of the system sizes
    
    % --- Precomputed value of the infinite chain ---
    w_def_int = 0.1124085417663; % nu = 1.5, gamma = 1, s1 = l1 = 0.5
    
    % --- Iterate over chain sizes ---
    res = ones(1,N);
    sizes = 3:N; % System sizes
    
    for n = sizes
        % --- Generate defected capacitance matrix ---
        defect_site = floor(n/2);
        l = l1 * ones(1, n);
        s = s1 * ones(1, n);
        capmat = Capacitance(n, s, gamma, l);
        D = eye(n);
        D(defect_site, defect_site) = 1 + nu;
        Def = D * capmat;
        
        % --- Compute the resonances ---
        [~, eigenvalues_matrix] = eig(Def);
        eigenvalues = sort(diag(eigenvalues_matrix));
        eigenvalues = abs(eigenvalues);
        resonances  = sqrt(eigenvalues * delta);
            
        w_def = resonances(n);
        res(n) = abs(w_def - w_def_int);
    end
    
    % --- Fit a regression line on the semilog plot ---
    log_res = log(res(3:end)); % Ignore initial ones
    p = polyfit(sizes, log_res, 1); % Linear fit: log(res) = p(1) * size + p(2)
    slope = p(1); % Extract slope
    
    % --- Plot results ---
    figure;
    semilogy(sizes, res(3:end), 'bo', 'LineWidth', 2, 'MarkerSize', 8);  
    hold on;
    
    % --- Plot regression line ---
    semilogy(sizes, exp(polyval(p, sizes)), 'r--', 'LineWidth', 2);
    
    xlabel('Size of the system', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('Frequency error', 'Interpreter', 'latex', 'FontSize', fs);
    set(gca, 'FontSize', fs+4, 'TickLabelInterpreter', 'latex');
    grid on;
    legend('Difference', 'Regression Line', 'Interpreter', 'latex', 'FontSize', fs);
    
    % Display slope
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
