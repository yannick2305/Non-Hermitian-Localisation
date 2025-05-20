%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [April 2025]
    Description:  [Localisation regions of monomer Band]
    --------------------------------------------------------------
%}

    clear all;
    close all;

% --- Define fixed parameters ---
    gamma = 3;          % Gauge potential
    delta = 0.001;      % Contrast parameter
    
    s1 = 0.5;           % Spacing betweeen the resonators
    l1 = 0.5;           % Length of the resonators
    L  = s1 + l1;       % Length of the unit cell
    Nx = 100;           % Number of plotting points in the bands
    fs = 14;            % Fontsize in plot annotation
    lw = 3.5;           % Linewidth of the spectral bands
    
    % --- Renormalise the unit cell ---
    s1 = s1 / L;
    l1 = l1 / L;
    L  = 1;

% --- Lower Gap function ---
    alpha_fixed = 0; 
    beta_lower_gap = linspace(0, gamma*l1, Nx);
    w_beta_lower_gap = (gamma * l1)/s1 * ( (1 - exp(+1i*(alpha_fixed + 1i* beta_lower_gap)*L)) / (1 - exp(-gamma * l1)) + (exp(-1i*(alpha_fixed + 1i* beta_lower_gap)*L)-1) / (1-exp(gamma*l1)) );
    w_beta_lower_gap = sqrt(delta *  abs(w_beta_lower_gap));

    % --- Limit of lower gap ---
    Lower_gap = (gamma * l1)/s1 * ( (1 - exp(+1i*(0 + 1i* (gamma*l1/2))*L)) / (1 - exp(-gamma * l1)) + (exp(-1i*(0 + 1i* (gamma*l1/2))*L)-1) / (1-exp(gamma*l1)) );
    Lower_gap = sqrt(delta * abs(Lower_gap));

% --- Upper Gap function ---
    alpha_fixed_2 = pi;
    beta_upper_gap = linspace(-3, 3, Nx);
    w_beta_upper_gap = (gamma * l1)/s1 * ( (1 - exp(+1i*(alpha_fixed_2 + 1i* beta_upper_gap)*L)) / (1 - exp(-gamma * l1)) + (exp(-1i*(alpha_fixed_2 + 1i* beta_upper_gap)*L)-1) / (1-exp(gamma*l1)) );
    w_beta_upper_gap = sqrt(delta * abs(w_beta_upper_gap));

% --- Band function ---
    beta_fixed =  gamma * l1/2;
    alpha = linspace(-pi, pi, Nx);
    w_alpha = (gamma * l1)/s1 * ( (1 - exp(+1i*(alpha + 1i* beta_fixed)*L)) / (1 - exp(-gamma * l1)) + (exp(-1i*(alpha + 1i* beta_fixed)*L)-1) / (1-exp(gamma*l1)) );
    w_alpha = sqrt(delta * abs(w_alpha));

    % --- Limit of upper gap ---
    Upper_gap = real((gamma * l1)/s1 * ( (1 - exp(+1i*(pi + 1i* beta_fixed)*L)) / (1 - exp(-gamma * l1)) + (exp(-1i*(pi + 1i* beta_fixed)*L)-1) / (1-exp(gamma*l1)) )); 
    Upper_gap = sqrt(delta * abs(Upper_gap));
    
    w_beta_0 = (gamma * l1)/s1 * ( (1 - exp(+1i*(alpha_fixed_2 + 1i* 0)*L)) / (1 - exp(-gamma * l1)) + (exp(-1i*(alpha_fixed_2 + 1i* 0)*L)-1) / (1-exp(gamma*l1)) );
    w_beta_0 = sqrt(delta * abs(w_beta_0));

    disp('----------------------------------------');
    disp(['Limit winding region : ', num2str(w_beta_0)]);

% --- Plot the spectral bands ---
    figure;
    xlimit = 7;
    x_coords = [-xlimit xlimit xlimit -xlimit];
    y_coords = [-0 0 w_beta_0 w_beta_0];
    y_coords_1 = [ w_beta_0 w_beta_0 4 4];
    y_coords_2 = [ Upper_gap Upper_gap Lower_gap Lower_gap];


    fill(x_coords, y_coords, [0.95 0.95 1]);      % lighter blue
    hold on;
    fill(x_coords, y_coords_1, [0.9 0.95 0.85] );  % lighter green
    fill(x_coords, y_coords_2, [0.95 0.8 0.9]);   % lighter pink/purple
 
    plot(beta_lower_gap, real(w_beta_lower_gap), 'r', 'LineWidth', lw);
  
    plot(alpha,          real(w_alpha),          'k', 'LineWidth', lw);
    plot(beta_upper_gap, real(w_beta_upper_gap), 'r', 'LineWidth', lw);

    % --- Plot fixed alpha and beta in the bands ---
    plot([0 0],                     [0 Lower_gap],             'k-', 'LineWidth', 1); 
    plot([pi pi],                   [Upper_gap 3 * Upper_gap], 'k-', 'LineWidth', 1);
    plot([-pi -pi],                 [Upper_gap 3 * Upper_gap], 'k-', 'LineWidth', 1);
    plot([gamma*l1/2, gamma*l1/2],  [Lower_gap  Upper_gap],    'r-', 'LineWidth', 1);

    xlim([-3.5, xlimit]);
    ylim([0, Upper_gap * 1.5]);

    % --- LaTeX labels and ticks ---
    xlabel('$\alpha$ and $\beta$ respectively', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\omega^{\alpha, \beta, \gamma}$' , 'Interpreter', 'latex', 'FontSize', fs);
    xticks([-pi/L, 0, pi/L]); 
    xticklabels({'$-\pi/L$', '$0$', '$\pi/L$'});
    set(gca, 'FontSize', fs+4, 'TickLabelInterpreter', 'latex');
    set(gcf, 'Position', [100, 100, 900, 400]); 

    grid off;

    % --- Inserts to illustrate eigenmode behavior ---
    sizeIn = 0.12;

    % --- Collapsed symbol region ---
        v = logspace(0, -14, 20); 
        idx = 1:length(v);       
        
        % --- Create inset axes ---
        inset_pos = [0.73 0.45 sizeIn sizeIn]; 
        axInset = axes('Position', inset_pos);
        box on;
    
        semilogy(idx, v, 'r-', 'LineWidth', 2);

        %xlabel('Index',     'Interpreter', 'latex', 'FontSize', 10);
        %ylabel('$|u(x)|$' , 'Interpreter', 'latex', 'FontSize', 10);
       
  % --- Winding region ---
        v1 = logspace(0, -2, 10); 
        v2 = logspace(-2.5, -10, 10);
        v = [v1, v2];
        idx = 1:length(v);       
        
        inset_pos = [0.73 0.18 sizeIn sizeIn]; 
        axInset = axes('Position', inset_pos);
        box on;
 
        semilogy(idx, v, 'b-', 'LineWidth', 2);

        %xlabel('Index',     'Interpreter', 'latex', 'FontSize', 10);
        %ylabel('$|u(x)|$' , 'Interpreter', 'latex', 'FontSize', 10);

 % --- Bukl localisation region ---
        v1 = logspace(-2, 0, 10);  
        v2 = logspace(-0.8, -9, 10);
        v = [v1, v2];
    
        idx = 1:length(v);      
        
        inset_pos = [0.73 0.755 sizeIn sizeIn];  
        axInset = axes('Position', inset_pos);
        box on;
        
        semilogy(idx, v, 'g-', 'LineWidth', 2);

        %xlabel('Index',     'Interpreter', 'latex', 'FontSize', 10);
        %ylabel('$|u(x)|$' , 'Interpreter', 'latex', 'FontSize', 10);
        
    hold off;
    exportgraphics(gcf, 'Spectral_plot_Regions_new.pdf');
