%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [April 2025]
    Description:  [Spectral plot monomer chain]
    --------------------------------------------------------------
%}

    clear all;
    close all;

% --- Define fixed parameters ---
    gamma = 3;          % Gauge potential
    delta = 0.001;      % Contrast parameter
    s1    = 0.5;        % Spacing betweeen the resonators
    l1    = 0.5;        % Length of the resonators
    L     = s1 + l1;    % Length of the unit cell
    Nx    = 200;        % Number of plotting points in the bands
    fs    = 18;         % Fontsize in plot annotation
    lw    = 2.5;        % Linewidth of spectral bands
    
    % --- Renormalise the lengths ---
    s1 = s1 / L;
    l1 = l1 / L;
    L  = 1;

% --- Lower Gap function ---
    alpha_fixed = 0; 
    beta_lower_gap = linspace(0, gamma*l1, Nx);
    w_beta_lower_gap = (gamma * l1)/s1 * ( (1 - exp(+1i*(alpha_fixed + 1i* beta_lower_gap)*L)) / (1 - exp(-gamma * l1)) + (exp(-1i*(alpha_fixed + 1i* beta_lower_gap)*L)-1) / (1-exp(gamma*l1)) );
    w_beta_lower_gap = sqrt(delta *  abs(w_beta_lower_gap));

    % --- Limit of lower gap ---
    Lower_gap = (gamma * l1)/s1 * ( (1 - exp(-1i*(0 + 1i* (-gamma*l1/2))*L)) / (1 - exp(-gamma * l1)) + (exp(1i*(0 + 1i* (-gamma*l1/2))*L)-1) / (1-exp(gamma*l1)) );
    Lower_gap = sqrt(delta * abs(Lower_gap));

% --- Upper Gap function ---
    alpha_fixed_2 = pi;
    beta_upper_gap = linspace(-3, 3, Nx);
    w_beta_upper_gap = (gamma * l1)/s1 * ( (1 - exp(+1i*(alpha_fixed_2 + 1i* beta_upper_gap)*L)) / (1 - exp(-gamma * l1)) + (exp(-1i*(alpha_fixed_2 + 1i* beta_upper_gap)*L)-1) / (1-exp(gamma*l1)) );
    w_beta_upper_gap = sqrt(delta * abs(w_beta_upper_gap));

% --- Band function ---
    beta_fixed = (gamma * l1/2);
    alpha = linspace(-pi, pi, Nx);
    w_alpha = (gamma * l1)/s1 * ( (1 - exp(+1i*(alpha + 1i* beta_fixed)*L)) / (1 - exp(-gamma * l1)) + (exp(-1i*(alpha + 1i* beta_fixed)*L)-1) / (1-exp(gamma*l1)) );
    w_alpha = sqrt(delta * abs(w_alpha));

    % --- Limit of upper gap ---
    Upper_gap = real((gamma * l1)/s1 * ( (1 - exp(+1i*(pi + 1i* beta_fixed)*L)) / (1 - exp(-gamma * l1)) + (exp(-1i*(pi + 1i* beta_fixed)*L)-1) / (1-exp(gamma*l1)) )); 
    Upper_gap = sqrt(delta * abs(Upper_gap));
    
% --- Symbol band function ---
    a = (gamma / s1) * (l1 / (1 - exp(-gamma * l1))) - (gamma / s1) * (l1 / (1 - exp(gamma * l1)));
    b =  gamma / s1  *  l1 / (1 - exp( gamma * l1));
    c = -gamma / s1  *  l1 / (1 - exp(-gamma * l1));

    % --- Alternative definition of band function ---
    f_z = - 2 * sqrt(b*c) * cos(alpha) + a;
    f_z = sqrt(delta * abs(f_z));
   
% --- Compute the decay of Green's function ---

    a = delta * a;
    b = delta * b;
    c = delta * c;

    % --- Decay to the left of the defect ---
    omega_1 = linspace(0, Upper_gap, 12 );

    r = 0.5 * log(b/c);
    d  = (a - omega_1.^2)/ (2 * sqrt(b*c));
    dl = real( -(r - acosh(d)) );
    dr = real( -(r + acosh(d)) );

    % --- Decay to the right of the defect ---
    omega_2 = linspace(Upper_gap, 1.5 * Upper_gap, 9);

    d   = (a - omega_2.^2)/ (2 * sqrt(b*c));
    dl2 = real( -(r - acosh(d)) );
    dr2 = real( -(r + acosh(d)) );

% --- Plot the spectral bands ---
    figure;
    plot(beta_lower_gap, real(w_beta_lower_gap), 'r', 'LineWidth', lw);
    hold on;
    plot(alpha,          real(f_z),              'k', 'LineWidth', lw);
    plot(beta_upper_gap, real(w_beta_upper_gap), 'r', 'LineWidth', lw);

    % --- Plot the decay of the Green's function ---
    %{
    plot(dl,  omega_1, 'bx', 'MarkerSize', 8, 'LineWidth', 2);
    plot(dr,  omega_1, 'bx', 'MarkerSize', 8, 'LineWidth', 2);
    plot(dl2, omega_2, 'bx', 'MarkerSize', 8, 'LineWidth', 2);
    plot(dr2, omega_2, 'bx', 'MarkerSize', 8, 'LineWidth', 2);
    %}

    % --- Plot fixed alpha and beta in the bands ---
    plot([0 0],                   [0 Lower_gap],             'k-', 'LineWidth', 1); 
    plot([pi pi],                 [Upper_gap 3 * Upper_gap], 'k-', 'LineWidth', 1);
    plot([-pi -pi],               [Upper_gap 3 * Upper_gap], 'k-', 'LineWidth', 1);
    plot([gamma*l1/2 gamma*l1/2], [Lower_gap  Upper_gap],    'r-', 'LineWidth', 1);

    % --- Mark limit of Spectrum/Gap ---
    yline(Upper_gap, 'k--', 'LineWidth', 1);
    yline(Lower_gap, 'k--', 'LineWidth', 1);
    ylim([0, Upper_gap * 1.3]);

    % --- LaTeX labels and ticks ---
    xlabel('$\alpha$ and $\beta$ respectively', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\omega^{\alpha, \beta, \gamma}$' , 'Interpreter', 'latex', 'FontSize', fs);
    
    % --- Set LaTeX-formatted ticks ---
    xticks([-pi/L, 0, pi/L]); 
    xticklabels({'$-\pi/L$', '$0$', '$\pi/L$'});
    set(gca, 'FontSize', fs+4, 'TickLabelInterpreter', 'latex');
    set(gcf, 'Position', [100, 100, 500, 400]); 

    grid off;
    hold off;
