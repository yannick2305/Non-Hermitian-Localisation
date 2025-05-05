%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [Mai 2025]
    Description:  [Spectral plot non-Hermitian Hamiltonian]
    --------------------------------------------------------------
%}

% --- Define the Parameters ---
    v     = 1;             
    gamma = 0.4;      
    fs    = 18;

    beta_fix     = - gamma;         
    alpha_fix_0  = 0;
    alpha_fix_pi = pi;
    
 % --- Define the range of the bands ---
    alpha = linspace(-pi, pi, 1000);  
    beta  = linspace(-5,   5, 1000);
    
% --- Compute the band functions ---
    f_alpha    = v - 2 * cos(alpha) .* cosh(gamma + beta_fix);
    f_beta_pi  = v - 2 * cos(alpha_fix_pi) .* cosh(gamma + beta);
    f_beta_0   = v - 2 * cos(alpha_fix_0) .* cosh(gamma + beta);
    
    % --- Compute the limit of the spectrum ---
    Upper_gap = v - 2 * cos(pi) .* cosh(gamma + beta_fix);
    Lower_gap = v - 2 * cos(0) .* cosh(gamma + beta_fix);

% --- Plot the Band functions ---
    figure;
    plot(alpha, f_alpha, 'k', 'LineWidth', 2);
    hold on;
    plot(beta,            f_beta_pi,                 'r-', 'LineWidth', 2);
    plot(beta,            f_beta_0,                  'r-', 'LineWidth', 2);
    plot([0 0],           [Lower_gap-3 Lower_gap],   'k-', 'LineWidth', 1); 
    plot([pi pi],         [Upper_gap 3 * Upper_gap], 'k-', 'LineWidth', 1);
    plot([-pi -pi],       [Upper_gap 3 * Upper_gap], 'k-', 'LineWidth', 1);
    plot([-gamma -gamma], [Lower_gap  Upper_gap],    'r-', 'LineWidth', 1);

    yline(Upper_gap, 'k--', 'LineWidth', 1);
    yline(Lower_gap, 'k--', 'LineWidth', 1);
    grid off;
    ylim([Lower_gap - 1, Upper_gap + 1]);
 % --- LaTeX labels and ticks ---
    xlabel('$\alpha$ and $\beta$ respectively', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\lambda^{\alpha, \beta, \gamma}$' , 'Interpreter', 'latex', 'FontSize', fs);
    
    % --- Set LaTeX-formatted ticks ---
    xticks([-pi, 0, pi]); 
    xticklabels({'$-\pi$', '$0$', '$\pi$'});
    set(gca, 'FontSize', fs+4, 'TickLabelInterpreter', 'latex');
    set(gcf, 'Position', [100, 100, 500, 400]); 
