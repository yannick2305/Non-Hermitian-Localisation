%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [April 2025]
    Description:  [Complex band structure surface plot monomer]
    --------------------------------------------------------------
%}

clear all;
close all;

% --- Define fixed parameters ---
    gamma  = 3;         % Gauge potential
    delta  = 0.001;     % Contrast parameter
    s1     = 0.5;       % Spacing betweeen the resonators
    l1     = 0.5;       % Length of the resonators
    L      = s1 + l1;   % Length of the unit cell
    lw     = 6;         % Linewidth of the spectral band
    fs     = 26;        % Fontsize plot annotation
    N_mesh = 50;        % Number of discretisation points for the surface
    N_line = 100;       % Number of discretisation points for spectral band

    % --- Initialise mesh ---
    alpha = linspace(-pi, pi, N_mesh);                  
    beta  = linspace(0,  1.5, N_mesh);
    [Alpha, Beta] = meshgrid(alpha, beta);

    z = Alpha - 1i * Beta;

    % --- Compute the spectral bands ---
    term1 = (1 - exp(-1i * z * L)) ./ (1 - exp(-gamma * l1));
    term2 = (exp(1i * z * L) - 1)  ./ (1 - exp( gamma * l1));
    w_beta_upper_gap = (gamma * l1) / s1 * (term1 + term2);
    w_beta_upper_gap = sqrt(w_beta_upper_gap * delta);
   
    real_part = real(w_beta_upper_gap);
    imag_part = imag(w_beta_upper_gap);

% --- Generate the surface plot ---
    figure;
    z_offset = 3e-3;
    
    mesh(Alpha, Beta, real_part - z_offset, imag_part);  
    xlabel('$\alpha$',                    'Interpreter', 'latex', 'FontSize', fs+4);
    ylabel('$\beta$',                     'Interpreter', 'latex', 'FontSize', fs+4);
    zlabel('$Re(\omega^{\alpha,\beta})$', 'Interpreter', 'latex', 'FontSize', fs+5);

    colorbar; 
    shading interp;  
    colormap jet;   
    cb = colorbar; 
    cb.Label.String = '$Im(\omega^{\alpha,\beta})$';
    cb.TickLabelInterpreter = 'latex'; 
    cb.FontSize = fs;
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize = fs; 
    cb.Label.Rotation = 270; 
    cb.Label.VerticalAlignment = 'bottom'; 
    set(gcf, 'Position', [100, 100, 1100, 600]);

    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = fs+5;

% --- Real spectral band ---
    beta_fixed = gamma * l1 /2;
    alpha = linspace(-pi, pi,  N_line);
    
    z_line = alpha - 1i * beta_fixed;
    term1_line  = (1 - exp(-1i * z_line * L)) ./ (1 - exp(-gamma * l1));
    term2_line  = (exp(1i * z_line * L) - 1)  ./ (1 - exp( gamma * l1));
    w_beta_line = (gamma * l1) / s1 * (term1_line + term2_line);
    w_beta_line = sqrt(delta * w_beta_line);

    hold on;
   
    plot3(alpha, beta_fixed * ones(size(alpha)), real(w_beta_line) , 'k-', 'LineWidth', lw);

% --- Plot gap band below ---

    % --- Band at alpha = pi ---
    beta = linspace(0, gamma/2, N_line);
    alpha_fixed = 0;

    z_line_2 = alpha_fixed - 1i * beta;
    
    term1_line  = (1 - exp(-1i * z_line_2 * L)) ./ (1 - exp(-gamma * l1));
    term2_line  = (exp(1i * z_line_2 * L) - 1)  ./ (1 - exp( gamma * l1));
    w_beta_line = (gamma * l1) / s1 * (term1_line + term2_line);
    w_beta_line = sqrt(delta * w_beta_line);

    plot3(alpha_fixed * ones(size(alpha)), beta,  real(w_beta_line), 'k-', 'LineWidth', lw);

% --- Plot gap band above ---
    add = 0.5;
    % --- Band at alpha = pi ---
    beta = linspace(- add, gamma/2 + add,  N_line);
    alpha_fixed = pi;

    z_line_2 = alpha_fixed - 1i * beta;
    
    term1_line = (1 - exp(-1i * z_line_2 * L)) ./ (1 - exp(-gamma * l1));
    term2_line = (exp(1i * z_line_2 * L) - 1) ./ (1 - exp(gamma * l1));
    w_beta_line = (gamma * l1) / s1 * (term1_line + term2_line);
    w_beta_line = sqrt(delta * w_beta_line);

    plot3(alpha_fixed * ones(size(alpha)), beta,  real(w_beta_line), 'k-', 'LineWidth', lw);
    
    % --- Band at alpha = -pi ---
    beta = linspace(- add, gamma/2 + add,  N_line);
    alpha_fixed = -pi;

    z_line_2 = alpha_fixed - 1i * beta;
    
    term1_line  = (1 - exp(-1i * z_line_2 * L)) ./ (1 - exp(-gamma * l1));
    term2_line  = (exp(1i * z_line_2 * L) - 1)  ./ (1 - exp( gamma * l1));
    w_beta_line = (gamma * l1) / s1 * (term1_line + term2_line);
    w_beta_line = sqrt(delta * w_beta_line);

    plot3(alpha_fixed * ones(size(alpha)), beta,  real(w_beta_line), 'k-', 'LineWidth', lw);

    view([-13, 51]);
    
    axis tight;

    %exportgraphics(gcf, 'Monomer_surface.pdf', 'ContentType', 'vector');
    exportgraphics(gcf, 'Monomer_surface_s.pdf', 'Resolution', 300);

