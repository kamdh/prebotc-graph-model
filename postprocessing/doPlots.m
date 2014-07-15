clear all
close all


set(0,'DefaultFigurePaperPositionMode','auto')
set(0,'DefaultAxesFontSize', 18)
set(0,'DefaultAxesFontName','Helvetica')
set(0,'DefaultAxesLineWidth',1);

dopartics = 1;
docombined = 0;
partics = {'er_n300_k1.0_deg_pI0.00_rep1',...
           'er_n300_k1.0_deg_pI0.20_rep1',...
           'er_n300_k1.0_deg_pI0.50_rep1',...
           'er_n300_k1.0_deg_pI0.80_rep1',...
           'er_n300_k1.0_deg_pI1.00_rep1'};
% partics = {
%     'er_n300_k3.0_deg_pI0.00_rep1', ...
%     'er_n300_k3.0_deg_pI0.10_rep1', ...
%     'er_n300_k3.0_deg_pI0.50_rep1', ...
%     'er_n300_k3.0_deg_pI0.80_rep1', ...
%     'er_n300_k3.0_deg_pI1.00_rep1', ...
%     'er_n300_k2.0_deg_pI0.00_rep1', ...
%     'er_n300_k2.0_deg_pI0.10_rep1', ...
%     'er_n300_k2.0_deg_pI0.50_rep1', ...
%     'er_n300_k2.0_deg_pI0.80_rep1', ...
%     'er_n300_k2.0_deg_pI1.00_rep1', ...
%     'er_n300_k6.0_deg_pI0.00_rep1', ...
%     'er_n300_k6.0_deg_pI0.10_rep1', ...
%     'er_n300_k6.0_deg_pI0.50_rep1', ...
%     'er_n300_k6.0_deg_pI0.80_rep1', ...
%     'er_n300_k6.0_deg_pI1.00_rep1' 
%           };
projName = 'random_fine_2';


%% start processing
dataDir = [getenv('HOME') '/work/prebotc/data/', projName]
fn = [dataDir, '/post/collected.mat'];
plotDir = [getenv('HOME') '/work/prebotc/data/', projName, ...
           '/plots'];

if docombined
    load(fn)

    figure
    myPcolor(X,Y, chiArray);
    title('\chi', 'fontsize', 32)
    xlabel('\langle k \rangle', 'fontsize', 24)
    ylabel('p_I','fontsize', 24)
    caxis([0,1])
    colorbar
    colormap('gray')
    plt = [plotDir, '/chi.eps']
    print('-deps', plt)

    figure
    myPcolor(X,Y, dutyCycle)
    title('duty cycle', 'fontsize', 32)
    xlabel('\langle k \rangle', 'fontsize', 24)
    ylabel('p_I','fontsize', 24)
    colorbar
    colormap('gray')
    plt = [plotDir, '/duty_cycle.eps']
    print('-deps', plt)

    figure
    myPcolor(X,Y, fMax)
    title('peak frequency (1/s)','fontsize', 32)
    xlabel('\langle k \rangle', 'fontsize', 24)
    ylabel('p_I','fontsize', 24)
    colorbar
    colormap('gray')
    plt = [plotDir, '/peak_freq.eps']
    print('-deps', plt)

    figure
    myPcolor(X,Y, lag)
    title('dominant period (s)','fontsize', 32)
    xlabel('\langle k \rangle', 'fontsize', 24)
    ylabel('p_I','fontsize', 24)
    colorbar
    colormap('gray')
    plt = [plotDir, '/lag.eps']
    print('-deps', plt)

    figure
    myPcolor(X,Y, muB / 1000)
    title('mean burst duration (s)','fontsize', 32)
    xlabel('\langle k \rangle','fontsize', 24)
    ylabel('p_I','fontsize', 24)
    colorbar
    colormap('gray')
    plt = [plotDir, '/mean_burst.eps']
    print('-deps', plt)

    figure
    myPcolor(X,Y, muIBI / 1000)
    title('mean IBI (s)','fontsize', 32)
    xlabel('\langle k \rangle','fontsize', 24)
    ylabel('p_I','fontsize', 24)
    colorbar
    colormap('gray')
    plt = [plotDir, '/mean_IBI.eps']
    print('-deps', plt)

    figure
    myPcolor(X,Y, cvIBI)
    title('c.v. of IBIs','fontsize', 32)
    xlabel('\langle k \rangle','fontsize', 24)
    ylabel('p_I','fontsize', 24)
    colorbar
    colormap('gray')
    plt = [plotDir, '/cv_IBIs.eps']
    print('-deps', plt)

    figure
    myPcolor(X,Y, cvB)
    title('c.v. of burst durations','fontsize', 32)
    xlabel('\langle k \rangle','fontsize', 24)
    ylabel('p_I','fontsize', 24)
    colorbar
    colormap('gray')
    plt = [plotDir, '/cv_bursts.eps']
    print('-deps', plt)


    figure
    myPcolor(X,Y, op_angle_mean)
    title('mean OP phase','fontsize', 32)
    xlabel('\langle k \rangle','fontsize', 24)
    ylabel('p_I','fontsize', 24)
    colorbar
    colormap('gray')
    plt = [plotDir, '/op_angle_mean.eps']
    print('-deps', plt)
    figure
    myPcolor(X,Y, op_angle_std)
    title('std dev OP phase','fontsize', 32)
    xlabel('\langle k \rangle','fontsize', 24)
    ylabel('p_I','fontsize', 24)
    colorbar
    colormap('gray')
    plt = [plotDir, '/op_angle_std.eps']
    print('-deps', plt)


    close all
end

if dopartics
    %%%%%%%%% plot some trajectories
    idx = 1;
    for partic=partics
        partic = char(partic);
        fn1 = [dataDir, '/output/', partic, '.mat'];
        fn2 = [dataDir, '/post/', partic, '_post.mat'];

        A = load(fn1, 'vTypes');
        vTypes = A.vTypes;
        clear A
        [Y,I] = sort(vTypes);

        
        B = load(fn2);
        trans = 12;
        binWidth = double( max(diff(B.bins)) ) / 1000;
        thebins = double( B.bins(trans:end-trans) )/1000;
        butterIntBin = B.butter_int_bin(trans:end-trans);
        binCt = B.spike_mat_bin(:,trans:end-trans);
        %peaklocs = B.b
        
        figure
        set(gcf, 'position', [1118         727        2675         500])
        %% subplot 1
        h1 = subplot(2,1,1);
        axis(h1)
        plot(thebins, butterIntBin, 'k-')
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
        plt = [plotDir, '/ts_combined_' partic '.eps']
        print('-deps', plt)
    end
end