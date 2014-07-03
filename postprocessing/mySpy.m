function mySpy(X,Y,C)
    h = imagesc(X, Y, C');
    axis xy;
    set(gca, 'ytick', Y(1,:))
    set(gca, 'xtick', X(:,1))
