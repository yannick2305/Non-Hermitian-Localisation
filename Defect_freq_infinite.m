%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [March 2025]
    Description:  [Compute defect frequency for a single defect]
    --------------------------------------------------------------
%}

clear all;
close all;

% --- Define fixed parameters ---
    gamma = 1;          % Gauge potential
    delta = 0.001;      % Contrast parameter
    
    nu = 1.5;           % Defect parameter (change in wavespeed)
    s1 = 0.5;           % Resonator spacing
    l1 = 0.5;           % Resonator length
    L  = s1 + l1;       % Unit cell length
    fs = 18;            % Fontsize in the plot
    Nw = 1000;          % Number of frequency steps

    % --- Parmeters for Adaptive Gaussian integration ---
    torr = 1e-15;   % Accuracy (absolute)
    N = 16;         % Order of the scheme

% --- Compute the spectral bands ---
    a = (gamma / s1) * (l1 / (1 - exp(-gamma * l1))) - (gamma / s1) * (l1 / (1 - exp(gamma * l1)));
    b =  gamma / s1  *  l1 / (1 - exp( gamma * l1));
    c = -gamma / s1  *  l1 / (1 - exp(-gamma * l1));
    
    a = delta * a;
    b = delta * b;
    c = delta * c;

    % --- Limit of lower gap ---
    Lower_gap = (gamma * l1)/s1 * ( (1 - exp(-1i*(0 + 1i* (-gamma*l1/2))*L)) / (1 - exp(-gamma * l1)) + (exp(1i*(0 + 1i* (-gamma*l1/2))*L)-1) / (1-exp(gamma*l1)) );
    Lower_gap = sqrt(delta * abs(Lower_gap));

    % --- Limit of upper gap ---
    Upper_gap = real((gamma * l1)/s1 * ( (1 - exp(-1i*(pi + 1i* -gamma*l1/2)*L)) / (1 - exp(-gamma * l1)) + (exp(1i*(pi + 1i* -gamma*l1/2)*L)-1) / (1-exp(gamma*l1)) )); 
    Upper_gap = sqrt(delta * abs(Upper_gap));

% --- Compute a solution ---
    % --- Set parameters for equation to solve: I(w) = C ---
    C = 2 * pi / nu;

    % --- Define the Eigenvalues of Capacitance ---
    g = @(alpha) (c * sqrt(b/c) * exp(-1j * alpha) + a + b * sqrt(c/b) * exp(1j * alpha));
    
    % --- plot the spectral band (debugging) ---
    %alpha = linspace(-pi, pi, 100);
    %plot( alpha, sqrt(g(alpha)) );

    % --- Upper Gap Range ---
    w = linspace(Upper_gap * 1.0001, 0.12, Nw); % Start just above the spectrum otherwise integral explodes

    % --- Lower Gap Range ---
    %w = linspace(0.0001, Lower_gap * 0.9999, Nw);

    % --- Evaluate the integral for the frequency range --- 
    integrals = arrayfun(@(w_value) adapgauss(@(alpha) (g(alpha)) ./ (w_value^2 - g(alpha)), -pi, pi, torr, N) - C, w);

    % --- Find the defect eigengrequency ---
    [~, idx] = min(abs(integrals));
    w_def = w(idx);
    disp('----------------------------------------');
    disp(['Defect Frequency: ', num2str(w_def)]);
    
% --- Illustrate the result ---
    figure;
    plot(w, integrals, 'LineWidth', 2);
    hold on;
    plot(w(idx), integrals(idx), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    xlabel('$w$',            'Interpreter', 'latex', 'FontSize', fs);
    ylabel('Integral Value', 'Interpreter', 'latex', 'FontSize', fs);
    set(gca, 'FontSize', fs+4, 'TickLabelInterpreter', 'latex');
    grid on;


%% --- Refined search ---

    % --- Refined search Range ---
    w_r = linspace(w(idx) - 0.00001, w(idx) + 0.00001, 5000);
    
    % --- Evaluate the integral for the refined frequency range --- 
    integrals_r = arrayfun(@(w_value) adapgauss(@(alpha) (g(alpha)) ./ (w_value^2 - g(alpha)), -pi, pi, torr, N) - C, w_r);
    
    % --- Find the defect eigengrequency ---
    [~, idx_r] = min(abs(integrals_r));
    w_def_r = w_r(idx_r);
    disp('----------------------------------------');
    disp(['Defect Frequency (refined): ', num2str(w_def_r)]);
    
% --- Illustrate the result ---
    figure;
    plot(w_r, integrals_r, 'LineWidth', 2);
    hold on;
    plot(w_r(idx_r), integrals_r(idx_r), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    xlabel('$w$',            'Interpreter', 'latex', 'FontSize', fs);
    ylabel('Integral Value', 'Interpreter', 'latex', 'FontSize', fs);
    set(gca, 'FontSize', fs, 'TickLabelInterpreter', 'latex');
    grid on;
