%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [December 2024]
    Description:  [1D Disorderes system for Skin effect]
    --------------------------------------------------------------
%}

close all;
clear all;

% --- Parameters ---
n = 10;         % Number of unit cells in the chain
l0 = 0.3;       % Length of a monomer resonator
s0 = 0.7;       % Spacing of a monomer resonator
len = l0;       % total length of resonator
gamma = 3;      % Gauge potential

delta = 0.001;  % Contrast


%% --- Generate random resonator chain ---

    % Remark: While it is possible to generate a completely random chain by
    %         uncommenting the code below, we refrain from doing so as the 
    %         decay due the intefrace is an a posteriori estimate, i.e. one
    %         has to know the dimer positions in order to estimate the decay
    %         lenght. A pregenerated resonaor chain can be found below.
    
    %{
    % --- Initialize arrays ---
    l = [];
    s = [];
    
    % --- Generate th disordered chain of resonators ---
        for i = 1:n
            if rand() < 0.9 
                % Monomer
                l = [l, l0]; 
                s = [s, s0]; 
            else
                % Dimer
                l = [l, l0/2, l0/2]; 
                s = [s, s0/2, s0/2];
            end
        end
    
    N = length(l);
    
    % Display results
    fprintf('Resonator lengths (l):\n');
    disp(l);
    
    fprintf('Spacings (a):\n');
    disp(s);
    %}

%% --- Generate deterministic resonator chain ---
    
    % --- Initialize arrays ---
    l = [];
    s = [];
    
    % --- Generate deterministic chain ---
    Def_pos = [floor((1*n)/2)]; %, floor((2*n)/4), floor((3*n)/4)]; %floor((2*n)/3); 
    for i = 1:n
        if i == Def_pos(1)   % || i == Def_pos(2) || i == Def_pos(3)
            % Dimer
            l = [l, l0/3, (l0 - l0/3)]; 
            s = [s, s0/3, (s0 - s0/3)]; 
        else
            % Monomer
            l = [l, l0]; 
            s = [s, s0]; 
        end
    end
    
    % --- Compute the decay length of the interface mode ---
        N = length(l);
        % --- Get the defect eigenfrequency ---
        capmat = Capacitance(N, s, gamma, l);
        rowSums = sum(capmat, 2);
        [~, eigenvalues_matrix] = eig(capmat);
        eigenvalues = sort(diag(eigenvalues_matrix));
    
        w_def = sqrt(delta * eigenvalues(N));
        decay_interface =  acosh((w_def^2 * s0) / (2 * delta) - 1);
    


%%  --- Plot the eigenvector entries and illustrate exponential decay ---

    plot_all_eigenvectors(N, s, gamma, l, len, decay_interface , Def_pos);


%% --- Plot the geometry of the resonator chain ---

    PlotGeometry(N, l, s);


%% --- Plot the eigenvalues ---

     PlotEigenvalues(N, s, gamma, l)


%% --- Defining functions ---


function plot_all_eigenvectors(N, s, gamma, l, len, decay, pos)

    lw = 3;     % Thickness of the resonators
    fs = 24;    % Fontsize of the annotations

    % --- Generate the eigenvectors of the Capacitance matrix ---
        capmat = Capacitance(N, s, gamma, l);
        [V, D] = eig(capmat);
        eigenvalues = diag(D);
        [~, sortIdx] = sort(eigenvalues);
        sortedEigenvectors = V(:, sortIdx);    
    
    % --- Plot the eigenvector entries ---
        figure;
        hold on;    

        %{
        for i = 1:(N-2)
            plot(1:N, (abs(sortedEigenvectors(:, i))), '-', 'Color', 0.65 * [1, 1, 1], 'LineWidth', lw/2.5, 'MarkerSize', 15);
        end
        
        x = 0:N;
        y = exp(- gamma * 0.5 * sum(len) * (1 / length(len)) * x);

        % Remark:   We are plotting the mode evaluated at the resonators,
        %           therefore the decay rate has to be devided by the number of
        %           resonators in the unit cell, i.e. length(len).

        plot(x, (y), 'r', 'LineWidth', lw); 
        %}
   
        % --- Uncomment to plot the defect ---
            %
            for i = N-3:(N)
                plot(1:N, (abs(sortedEigenvectors(:, i))), '-', 'Color', 0.65 * [1, 1, 1], 'LineWidth', lw * 1.5, 'MarkerSize', 15);
            end
            % --- Add the predicted exponential decay rate ---
            x = 0:N; 
            constant = 1;
            %y =  constant * exp(- (gamma * 0.5 * len - decay) * (x - pos(1) - 1)); 
            %plot(x, (y), 'r', 'LineWidth', lw/2); 
            %y2 =  constant * exp(- (gamma * 0.5 * len + decay) * (x - pos(1) -1)); 
            %plot(x, (y2), 'r', 'LineWidth', lw/2); 
            
            %x2 = 0 : pos(2) + 2;
            %y3 =  constant * exp(- (gamma * 0.5 * len - decay) * (x2 - pos(2) - 1) - 9); 
            %plot(x2, (y3), 'r', 'LineWidth', lw/2); 
    
            %x3 = pos(2) +2 : N;
            %y4 =  constant * exp(- (gamma * 0.5 * len + decay) * (x3 - pos(2) -1) - 6); 
            %plot(x3, (y4), 'r', 'LineWidth', lw/2); 
            %}

    % --- Formatting the plot ---
        xlabel('Eigenvector index', 'Interpreter', 'latex', 'FontSize', fs);
        ylabel('$|u(x)|$', 'Interpreter', 'latex', 'FontSize', fs);
        ylim([-25, 1]);
        xlim([0, floor(N/10)*10]);
        grid on;
        set(gca, 'YScale', 'log');
        set(gca, 'FontSize', fs-4);
        hold off;
end

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

function PlotGeometry(N, l, s)

    x_start = 0;  % Initialize starting position of the first segment            

    figure;
    hold on;
    for i = 1:N
        x_end = x_start + l(i);
        line([x_start, x_end], [0, 0], 'Color', 'k', 'LineWidth', 4);
        if i < N
            x_start = x_end + s(i);
        end
    end

    % --- Formatting the plot ---
    set(gcf, 'Position', [100, 100, 600, 100]);  
    xlabel('x-position', 'Interpreter', 'latex', 'FontSize', 20);
    set(gca, 'FontSize', 18-4);
    set(gca, 'YTick', []); 
    ylim([-1,1]);
    axis equal;  
    grid off;  
    hold off;
end

function PlotEigenvalues(N, s, gamma, l)
    
    % --- Generate the Capacitance matrix ---
    capmat = Capacitance(N, s, gamma, l);
    
    % --- Compute the eigenvalues of the capacitance matrix ---
    [~, eigenvalues_matrix] = eig(capmat);
    eigenvalues = sort(diag(eigenvalues_matrix));  
    
    % --- Plot the eigenvalues ---
    figure;
    plot(eigenvalues, 'bo', 'LineWidth', 2, 'MarkerSize', 8);  
    xlabel('Index', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Eigenvalue', 'Interpreter', 'latex', 'FontSize', 14);
    grid on; 
end
