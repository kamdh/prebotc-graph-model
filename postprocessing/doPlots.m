clear all
close all


set(0,'DefaultFigurePaperPositionMode','auto')
set(0,'DefaultAxesFontSize', 18)
set(0,'DefaultAxesFontName','Helvetica')
set(0,'DefaultAxesLineWidth',1);

partics = {
    'er_n300_k3.0_deg_pI0.00_rep1', ...
    'er_n300_k3.0_deg_pI0.20_rep1', ...
    'er_n300_k3.0_deg_pI0.50_rep1', ...
    'er_n300_k3.0_deg_pI0.80_rep1', ...
    'er_n300_k3.0_deg_pI1.00_rep1', ...
    'er_n300_k1.0_deg_pI0.00_rep1', ...
    'er_n300_k1.0_deg_pI0.20_rep1', ...
    'er_n300_k1.0_deg_pI0.50_rep1', ...
    'er_n300_k1.0_deg_pI0.80_rep1', ...
    'er_n300_k1.0_deg_pI1.00_rep1', ...
    'er_n300_k6.0_deg_pI0.00_rep1', ...
    'er_n300_k6.0_deg_pI0.20_rep1', ...
    'er_n300_k6.0_deg_pI0.50_rep1', ...
    'er_n300_k6.0_deg_pI0.80_rep1', ...
    'er_n300_k6.0_deg_pI1.00_rep1' 
          };
projName = 'random_fine_2';


%% start processing
dataDir = [getenv('HOME') '/work/prebotc/data/', projName]
fn = [dataDir, '/post/collected.mat'];
plotDir = [getenv('HOME') '/work/prebotc/data/', projName, ...
           '/plots'];
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
    butterIntBin = B.butterIntBin(trans:end-trans);
    binCt = B.spikeMatBin(:,trans:end-trans);

    figure
    h1= subplot(2,1,1);
    plot(thebins, butterIntBin, 'k-')
    ylabel('\int preBot (Hz/neuron)', 'fontsize', 20)
    tmp = get(gca, 'ylim');
    %set(h1, 'ylim', [0, 35]);
    set(h1, 'ylim', [0, tmp(2)]);
    %xlabel('t (s)')
    %set(gcf, 'position', [1118         727        2675         321])
    %plt = [plotDir, '/ts_int.eps']
    %print('-deps', plt)
    h2 = subplot(2,1,2);
    imagesc(binCt(I,:));
    xlabel('t (s)', 'fontsize', 20)
    colormap(flipud(colormap('gray')))
    set(h2, 'xticklabel', get(h1, 'xticklabel'))
    %spy(binCt(I,:), 'k', 1)
    %xlabel(['time bin (' num2str(binWidth*1000) ' ms)'])
    ylabel('sorted neurons', 'fontsize', 20)
    set(gcf, 'position', [1118         727        2675         500])
    %plt = [plotDir, '/ts_raster.eps']
    %print('-deps', plt)
    plt = [plotDir, '/ts_combined_' partic '.eps']
    print('-deps', plt)
end