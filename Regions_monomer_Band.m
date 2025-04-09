%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [March 2025]
    Description:  [regions of monomer Band]
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
    Nx = 10000;           % Number of plotting points in the bands
    fs = 16;            % Fontsize in plot annotation
    lw = 3.5;
    
    % --- Renormalise the lengths ---
    s1 = s1 / L;
    l1 = l1 / L;
    L  = 1;

% --- Lower Gap function ---
    alpha_fixed = 0; 
    beta_lower_gap = linspace(-gamma*l1 , 0, Nx);
    w_beta_lower_gap = (gamma * l1)/s1 * ( (1 - exp(-1i*(alpha_fixed + 1i* beta_lower_gap)*L)) / (1 - exp(-gamma * l1)) + (exp(1i*(alpha_fixed + 1i* beta_lower_gap)*L)-1) / (1-exp(gamma*l1)) );
    w_beta_lower_gap = sqrt(delta *  abs(w_beta_lower_gap));

    % --- Limit of lower gap ---
    Lower_gap = (gamma * l1)/s1 * ( (1 - exp(-1i*(0 + 1i* (-gamma*l1/2))*L)) / (1 - exp(-gamma * l1)) + (exp(1i*(0 + 1i* (-gamma*l1/2))*L)-1) / (1-exp(gamma*l1)) );
    Lower_gap = sqrt(delta * abs(Lower_gap));

% --- Upper Gap function ---
    alpha_fixed_2 = pi;
    beta_upper_gap = linspace(-3, 3, Nx);
    w_beta_upper_gap = (gamma * l1)/s1 * ( (1 - exp(-1i*(alpha_fixed_2 + 1i* beta_upper_gap)*L)) / (1 - exp(-gamma * l1)) + (exp(1i*(alpha_fixed_2 + 1i* beta_upper_gap)*L)-1) / (1-exp(gamma*l1)) );
    w_beta_upper_gap = sqrt(delta * abs(w_beta_upper_gap));

% --- Band function ---
    beta_fixed = - (gamma * l1/2);
    alpha = linspace(-pi, pi, Nx);
    w_alpha = (gamma * l1)/s1 * ( (1 - exp(-1i*(alpha + 1i* beta_fixed)*L)) / (1 - exp(-gamma * l1)) + (exp(1i*(alpha + 1i* beta_fixed)*L)-1) / (1-exp(gamma*l1)) );
    w_alpha = sqrt(delta * abs(w_alpha));

    % --- Limit of upper gap ---
    Upper_gap = real((gamma * l1)/s1 * ( (1 - exp(-1i*(pi + 1i* beta_fixed)*L)) / (1 - exp(-gamma * l1)) + (exp(1i*(pi + 1i* beta_fixed)*L)-1) / (1-exp(gamma*l1)) )); 
    Upper_gap = sqrt(delta * abs(Upper_gap));
    
% --- Symbol band function ---
    a = (gamma / s1) * (l1 / (1 - exp(-gamma * l1))) - (gamma / s1) * (l1 / (1 - exp(gamma * l1)));
    b =  gamma / s1  *  l1 / (1 - exp( gamma * l1));
    c = -gamma / s1  *  l1 / (1 - exp(-gamma * l1));

    %f_z = c * sqrt(b/c) * exp(-1j * alpha) + a + b * sqrt(c/b) * exp(1j * alpha);
    f_z = -2 * sqrt(b*c) * cos(alpha) + a;
    f_z = sqrt(delta * abs(f_z));
   
% --- Compute the decay of Green's function ---

    a = delta * a;
    b = delta * b;
    c = delta * c;

    % --- Decay to the left of the defect ---
    omega_1 = linspace(0, Lower_gap, 3 );

    d  = (a - omega_1.^2)/ (2 * sqrt(b*c));
    dl = real(- (acosh(d) + log(sqrt(b*c)) - log(abs(b))));
    dr = real(   acosh(d) + log(sqrt(b*c)) - log(abs(c)));

    % --- Decay to the right of the defect ---
    omega_2 = linspace(Upper_gap, 1.5 * Upper_gap, 9);

    d   = (a - omega_2.^2)/ (2 * sqrt(b*c));
    dl2 = real(-(acosh(d) + log(sqrt(b*c)) - log(abs(b)) ));
    dr2 = real(  acosh(d) + log(sqrt(b*c)) - log(abs(c)));



    w_beta_0 = (gamma * l1)/s1 * ( (1 - exp(-1i*(alpha_fixed_2 + 1i* 0)*L)) / (1 - exp(-gamma * l1)) + (exp(1i*(alpha_fixed_2 + 1i* 0)*L)-1) / (1-exp(gamma*l1)) );
    w_beta_0 = sqrt(delta * abs(w_beta_0));

    format long
    disp('w_beta_0 =');
    disp(w_beta_0);


% --- Plot the spectral bands ---
    figure;
    xlimit = 6;
    x_coords = [-xlimit xlimit xlimit -xlimit];
    y_coords = [-0 0 w_beta_0 w_beta_0];
    y_coords_1 = [ w_beta_0 w_beta_0 4 4];
    y_coords_2 = [ Upper_gap Upper_gap Lower_gap Lower_gap];


    fill(x_coords, y_coords, [0.95 0.95 1]);      % lighter blue
    hold on;
    fill(x_coords, y_coords_1, [0.9 0.95 0.85] );  % lighter green
    fill(x_coords, y_coords_2, [0.95 0.8 0.9]);   % lighter pink/purple

    plot(beta_lower_gap, real(w_beta_lower_gap), 'r', 'LineWidth', lw);
    
    plot(alpha,          real(w_alpha),          'k',  'LineWidth', lw);
    plot(beta_upper_gap, real(w_beta_upper_gap), 'r',  'LineWidth', lw);

    % --- Plot the decay of the Green's function ---
    %{
    plot(dl,  omega_1, 'bx', 'MarkerSize', 8, 'LineWidth', 2);
    plot(dr,  omega_1, 'bx', 'MarkerSize', 8, 'LineWidth', 2);
    plot(dl2, omega_2, 'bx', 'MarkerSize', 8, 'LineWidth', 2);
    plot(dr2, omega_2, 'bx', 'MarkerSize', 8, 'LineWidth', 2);
    %}

    % --- Plot fixed alpha and beta in the bands ---
    plot([0 0],                     [0 Lower_gap],             'k-', 'LineWidth', 1); 
    plot([pi pi],                   [Upper_gap 3 * Upper_gap], 'k-', 'LineWidth', 1);
    plot([-pi -pi],                 [Upper_gap 3 * Upper_gap], 'k-', 'LineWidth', 1);
    plot([-gamma*l1/2 -gamma*l1/2], [Lower_gap  Upper_gap],    'r-', 'LineWidth', 1);

    % --- Plot the symbol bands ---
    %plot(alpha,          real(f_z),          'r:',  'LineWidth', 2);

    % --- Verify that imaginary part vanishes (debugging) ---
    %plot(beta_lower_gap, imag(w_beta_lower_gap), 'r--', 'LineWidth', 2);
    %plot(beta_upper_gap, imag(w_beta_upper_gap), 'r--', 'LineWidth', 2);
    %plot(alpha,          imag(w_alpha),          'k--', 'LineWidth', 2);

    % --- Mark limit of Spectrum/Gap ---
    %yline(Upper_gap, 'k--', 'LineWidth', 1);
    %yline(Lower_gap, 'k--', 'LineWidth', 1);

    xlim([-3.5, xlimit]);
    ylim([0, Upper_gap * 1.4]);

    %{

 yline(0.097193636630524, 'b--', 'LineWidth', 1);
    % Define alpha values over [-pi, pi]
    alpha_wind = linspace(-pi, pi, 10000);
    
    % Define unit circle points: z = e^{i*alpha}
    z = exp(1i * alpha_wind);
    
    % Define your complex function f(z)
    % Example: f(z) = (z^2 + 1)/(z - 0.5)
    f = @(z) c*z + a + b ./ z;   % <-- Change this to your desired function
    
    % Evaluate f(z)
    fz = sqrt(f(z));
    
    % Plot the path of f(z) in the complex plane
    plot(100 *imag(fz), real(fz), 'b', 'LineWidth', 2);

    %}
    % --- LaTeX labels and ticks ---
    xlabel('$\alpha$ and $\beta$ respectively', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\omega^{\alpha, \beta, \gamma}$' , 'Interpreter', 'latex', 'FontSize', fs);
    
    % --- Set LaTeX-formatted ticks ---
    xticks([-pi/L, 0, pi/L]); 
    xticklabels({'$-\pi/L$', '$0$', '$\pi/L$'});
    set(gca, 'FontSize', fs+4, 'TickLabelInterpreter', 'latex');
    set(gcf, 'Position', [100, 100, 900, 450]);

    
   

    grid off;

    sizeIn = 0.12;
% --- Inserts ---

    % Sample vector for inset plot
    v = logspace(0, -14, 20);  % Just an example: decreasing vector
    idx = 1:length(v);        % Indices
    
    % Create inset axes
    inset_pos = [0.73 0.45 sizeIn sizeIn];  % [left bottom width height], normalized
    axInset = axes('Position', inset_pos);
    box on;
    
    % Inset semilogy plot
    semilogy(idx, v, 'r-', 'LineWidth', 2);
    %title('Inset: semilogy');
    %xlabel('Index', 'Interpreter', 'latex', 'FontSize', 10);
    %ylabel('$|u(x)|$' , 'Interpreter', 'latex', 'FontSize', 10);
       
  % --- 2 ---
    % Sample vector for inset plot
    v1 = logspace(0, -2, 10);  % Just an example: decreasing vector
    v2 = logspace(-2.5, -10, 10);
    v = [v1, v2];
    idx = 1:length(v);        % Indices
    
    % Create inset axes
    inset_pos = [0.73 0.18 sizeIn sizeIn];  % [left bottom width height], normalized
    axInset = axes('Position', inset_pos);
    box on;
    
    % Inset semilogy plot
    semilogy(idx, v, 'b-', 'LineWidth', 2);
    %title('Inset: semilogy');
    %xlabel('Index', 'Interpreter', 'latex', 'FontSize', 10);
    %ylabel('$|u(x)|$' , 'Interpreter', 'latex', 'FontSize', 10);

 % --- 3 ---
  % Sample vector for inset plot
    v1 = logspace(-2, 0, 10);  % Just an example: decreasing vector
    v2 = logspace(-0.8, -9, 10);
    v = [v1, v2];

    idx = 1:length(v);        % Indices
    
    % Create inset axes
    inset_pos = [0.73 0.78 sizeIn sizeIn];  % [left bottom width height], normalized
    axInset = axes('Position', inset_pos);
    box on;
    
    % Inset semilogy plot
    semilogy(idx, v, 'g-', 'LineWidth', 2);
    %title('Inset: semilogy');
    %xlabel('Index', 'Interpreter', 'latex', 'FontSize', 10);
    %ylabel('$|u(x)|$' , 'Interpreter', 'latex', 'FontSize', 10);


    hold off;
    exportgraphics(gcf, 'Spectral_plot_inlay.pdf');


