%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [March 2025]
    Description:  [Defected Dimer system]
    --------------------------------------------------------------
%}


clear all;
close all;

n = 400; 
gamma = 2;                  % Gauge potential
delta = 0.001;              % Contrast parameter
nu = 1.5;                   % Change in wavespeed
ell   = [0.25, 0.25];       % Resonator lengths
spa     = [1, 2];           % Resonator spacings
defect_site = floor(n/3);

fs = 18;

spa = (0.5 / (spa(1) + spa(2))) * spa;

s = [];
l = [];

for i = 1:n
    l = [l, ell(1), ell(1)]; 
    s = [s, spa(1), spa(2)];
end

N = length(l);

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

    beta_fixed = -sum(ell) * gamma / 2;
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
    disp(['Defect Frequency 1: ', num2str(def(1))]);
    disp(['Defect Frequency 2: ', num2str(def(2))]);



% --- Discrete Green's function ---
    N = 400;
    cap = delta * Capacitance(N, s, gamma,  l);
    
    dirac = zeros(1, N);
    dirac(floor(N/2)) = 1;
    dirac = dirac.';

    N_freq = 30;

    w = linspace(0.001, 0.15, N_freq);

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

        decay_l(i) = p_l(1);
        decay_r(i) = p_r(1);

    end

    figure;
    plot(decay_l, w);
    hold on;
    plot(decay_r, w);
    hold off;


% --- Plot the spectrum ---
    figure;
    h1 = plot(resonances, 'ko', 'LineWidth', 4, 'MarkerSize', 8); 
    hold on;
    if nu > 0
        % --- Defect in the top spectral Gap ---
        h2 = plot(n,   def(1), 'rx', 'LineWidth', 4, 'MarkerSize', 8);
        h2 = plot(2*n, def(2), 'rx', 'LineWidth', 4, 'MarkerSize', 8);
    else
        % --- Defect in the lower spectral Gap ---
        h2 = plot(2,   def(1), 'rx', 'LineWidth', 4, 'MarkerSize', 8);
        h2 = plot(n+1, def(2), 'rx', 'LineWidth', 4, 'MarkerSize', 8);
    end
    yline(double(lim1(1)), 'b--', 'LineWidth', 1);
    yline(double(lim2(1)), 'b--', 'LineWidth', 1);
    yline(double(lim1(2)), 'b--', 'LineWidth', 1);
    yline(double(lim2(2)), 'b--', 'LineWidth', 1);

    xlabel('Index',    'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\omega$', 'Interpreter', 'latex', 'FontSize', fs);
    set(gca, 'FontSize', fs+4, 'TickLabelInterpreter', 'latex');
    legend([h1, h2], {'Band Resonance', 'Defect Resonance'}, ...
        'Interpreter', 'latex', ...
    'Location', 'southeast', ...
    'Box', 'on');
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
    C(2,1) = C(2,1) - exp(-1i*(alpha + 1i*beta)*L) * (gamma / s(2)) * ell(1)/(1-exp(-gamma* ell(1)));

    % i = 1, j = 2
    C(1,2) = C(1,2) + exp( 1i*(alpha + 1i*beta)*L) * (gamma / s(2)) * ell(2)/(1 - exp(gamma*ell(2)));
end
