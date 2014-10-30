function plotPartic(fnPost, plotDir, partic, varargin)
    if nargin == 3
        opThresh = 0.2;
    else
        opThresh = varargin{1};
    end

    B = load(fnPost);
    [Y,I]=sort(B.vertex_types);
    trans = 12;
    binWidth = double( max(diff(B.bins)) ) / 1000;
    thebins = double( B.bins(trans:end) )/1000;
    butterIntBin = B.butter_int_bin(trans:end);
    binCt = B.spike_mat_bin(:,trans:end);
    ops = B.ops;
    opAngle = angle(ops);
    opAbs = abs(ops);
    peaklocs = B.pop_burst_peak;
    maxindex = length(B.bins)-trans;
    peaklocs(peaklocs < trans | peaklocs > maxindex) = [];
    peaklocs = peaklocs - trans + 2; % +1 for trans, +1 for
                                     % python 0-indexing
    
    figure
    set(gcf, 'position', [1118, 727, 2675, 500])
    %% subplot 1
    h1 = subplot(2,1,1);
    axis(h1)
    plot(thebins, butterIntBin, 'k-')
    hold on
    plot(thebins(peaklocs), butterIntBin(peaklocs), 'ko')
    ylabel('\int preBot (Hz/neuron)', 'fontsize', 20)
    tmp = get(gca, 'ylim');
    axis tight
    set(h1, 'ylim', [0, tmp(2)]);
    %% subplot 2
    h2 = subplot(2,1,2);
    axis(h2)
    imagesc(thebins, 1:size(binCt,1), binCt(I,:));
    xlabel('t (s)', 'fontsize', 20)
    colormap(flipud(colormap('gray')))
    set(h2, 'xticklabel', get(h1, 'xticklabel'))
    %spy(binCt(I,:), 'k', 1)
    %xlabel(['time bin (' num2str(binWidth*1000) ' ms)'])
    ylabel('sorted neurons', 'fontsize', 20)
    %% finalize and print
    set(h2, 'xtick', get(h1,'xtick'))
    set(h2, 'xticklabel', get(h1,'xticklabel'))
    %plt = [plotDir, '/ts_raster.eps']
    %print('-deps', plt)
    plt = [plotDir, '/ts_combined_' partic '.eps'];
    print('-deps', plt)
    %% order parameter plot
    figure
    plot(opAngle, opAbs, 'k.')
    xlabel('angle (radians)')
    ylabel('magnitude')
    hline(opThresh, 'k--')
    ylim([0, 1])
    xlim([-pi, pi])
    title('neuron order parameters')
    plt = [plotDir, '/op_' partic '.eps'];
    print('-deps', plt)
    