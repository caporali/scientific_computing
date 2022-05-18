% preset: clear windows, load files and set font for plots
clear; clc;
load('../input/data_GR_Morrison.mat');
addpath('./methods');
addpath('./functions/input');
addpath('./functions/output');
addpath('./functions/side');
set(0,'DefaultAxesTitleFontWeight','normal');

% set
set = 1;
alpha = 0.85;
tol = 1e-14;
    % All
    if set == 1 
        W = sparse(w_All);
        ex = sparse(expr_data);	
    % Down
    elseif set == 2
        W = sparse(w_Down);
        ex = sparse(expr_dataDown);	
    % Up
    elseif set == 3
        W = sparse(w_Up);
        ex = sparse(expr_dataUp);	
    % Random
    elseif set == 4
        n = 300000;
        W = random_GR_matrix(n,7);
        ex = ones(n, 1)/n;
        W = sparse(W);
        ex = sparse(ex);
        
        % save the matrix and the vector
        % save('../input/random_300000.mat','W','ex');
    % Saved Random (load the matrix and the vector)
    elseif set == 5
        W = load('../input/random_300000.mat').W;
        ex = load('../input/random_300000.mat').ex;        
    end

% show input
fprintf("--- \ninput \n\n");

figure;
spy(W);
title('$$W$$','interpreter','latex');
xlabel(['$$W \in \{0, 1\}^{',int2str(size(W,1)),' \times ',int2str(size(W,2)),'}$$ / n. of non-zero elements = ',int2str(full(sum(sum(W))))],'interpreter','latex');
ax = gca;
ax.XAxis.FontSize = 8;
ax.YAxis.FontSize = 8;
ax.XLabel.FontSize = 10;
ax.Title.FontSize = 12;

% text output
fprintf("	size(W) =  %d x %d \n", size(W,1), size(W,2));
fprintf("	n. of non-zero elements =  %d \n", full(sum(sum(W))));

% main
    % remove dangling nodes
    W = remove_dangling(W);

    xreal = check_GR(W,ex,alpha);

    figure;

    % 1
    x = cg_GR(W,ex,alpha,tol);
    norm_diff(x,xreal);
    sin_angle(x, xreal);

    % 2
    x = pcg_GR(W,ex,alpha,tol);
    norm_diff(x,xreal);
    sin_angle(x, xreal);

    % 3
    x = jacobi_GR(W,ex,alpha,tol);
    norm_diff(x,xreal);
    sin_angle(x, xreal);

    % 4
    x = modified_arnoldi_GR(W,ex,alpha,3,tol,1000);
    norm_diff(x,xreal);
    sin_angle(x, xreal);

    % 5
    x = gauss_seidel_GR(W,ex,alpha,tol);
    norm_diff(x,xreal);
    sin_angle(x, xreal);

    % 6
    x = power_GR(W,ex,alpha,tol);
    norm_diff(x,xreal);
    sin_angle(x, xreal);

    % 7
    x = richardson_GR(W,ex,alpha,tol);
    norm_diff(x,xreal);
    sin_angle(x, xreal);

% plot sets
legend("cg\_GR", "pcg\_GR", "jacobi\_GR", "modified\_arnoldi\_GR", "gauss\_seidel\_GR", "power\_GR", "richardson\_GR");

% eigenvalues    
plot_eigenvalues(W,alpha);