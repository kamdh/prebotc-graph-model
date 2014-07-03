function myPcolor(X,Y,C)
    h = imagesc(X(:,1), Y(1,:), C');
    axis xy;
    % ylabels = linspace(min(Y(1,:)), max(Y(1,:)), 6);
    % xlabels = linspace(min(X(:,1)), max(X(:,1)), 6);
    % set(gca, 'ytick', ylabels)
    % set(gca, 'xtick', xlabels)
