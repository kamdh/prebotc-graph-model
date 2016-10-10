clear all
close all 

set(0,'DefaultFigurePaperPositionMode','auto')

%set(0,'DefaultFigurePosition',[25,25,500,500]);
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultAxesFontName','Helvetica')
set(0,'DefaultAxesLineWidth',1);

 
fn = [getenv('HOME'), '/prebotc-graph-model/data/2d_sweep/retry_collection_update'];
fn_messy = [getenv('HOME'), '/prebotc-graph-model/data/2d_sweep/post/er_n600_intra2.0_inter1.5_deg_pI0.50_rep5_gE2.5_gI2.5_post.mat'];
fn_clean = [getenv('HOME'), '/prebotc-graph-model/data/2d_sweep/post/er_n600_intra1.0_inter4.0_deg_pI0.50_rep7_gE2.5_gI2.5_post.mat'];
plot_fn = [getenv('HOME'), '/prebotc-graph-model/data/2d_sweep/'];
load(fn)
sorted_data = sortrows(phi_avrg_std,[1,2,7]);
mat_data = cell2mat(sorted_data(:,[3:6]));
mean_angle = mat_data(:,[1]);
variance_angle = mat_data(:,[2]);
chi_1 = mat_data(:,[3]);
chi_2 = mat_data(:,[4]);
avg_mean = reshape(mean(reshape(mean_angle,8,[])),[9,9]);
avg_var = reshape(mean(reshape(variance_angle,8,[])),[9,9]);
avg_chi1 = reshape(mean(reshape(chi_1,8,[])),[9,9]);
avg_chi2 = reshape(mean(reshape(chi_2,8,[])),[9,9]);
avg_chi = (avg_chi1 + avg_chi2) / 2;
X = 0:.5:4;
Y = 0:.5:4;

figure 
myPcolor(X,Y,transpose(avg_mean),'Spectral','clim',[.2,.8])
set(gca,'FontSize',20);
%titlestr = sprintf('Phase difference       ');
%title(titlestr,'fontsize',24,'FontWeight','normal');
xlabel('k_{intra}','Fontsize',28);
ylabel('k_{inter}','Fontsize',28);
set(gca, 'XTick', [0,1,2,3,4,5,6,7,8]);
set(gca, 'XTickLabel',[0,.5,1,1.5,2,2.5,3,3.5,4]);
plt = [plot_fn,'mean_average.eps'];
print('-depsc',plt);

figure
myPcolor(X,Y,transpose(avg_var),'YlOrRd','clim',[0,1.0])
title('Phase order        ','fontsize',24,'FontWeight','normal');
xlabel('k_{intra}','Fontsize',28);
ylabel('k_{inter}','Fontsize',28);
set(gca, 'XTick', [0,1,2,3,4,5,6,7,8]);
set(gca, 'XTickLabel',[0,.5,1,1.5,2,2.5,3,3.5,4]);
plt = [plot_fn,'mean_variance.eps'];
print('-depsc',plt);

figure
myPcolor(X,Y,transpose(avg_chi),'YlOrRd','clim',[0,1.0])
%title('Synchrony      ','fontsize',24,'FontWeight','normal');
xlabel('k_{intra}','Fontsize',28);
ylabel('k_{inter}','Fontsize',28);
set(gca, 'XTick', [0,1,2,3,4,5,6,7,8]);
set(gca, 'XTickLabel',[0,.5,1,1.5,2,2.5,3,3.5,4]);
plt = [plot_fn,'mean_chi.eps'];
print('-depsc',plt);

load(fn_clean);
figure
set(gcf,'position',[1118,727,1605,500]);
h1 = subplot(2,1,1)
axis(h1)
plot(bins/1000,butter_int_bin,'--k');
hold on
plot(bins/1000,butter_int_bin2,'k');
ylabel('x^{\rm int} (Hz/cell)','fontsize',36)
axis tight
ylim([0,20])
axis([0,12,0,20])
tmp = get(gca,'ylim');
set(h1,'ylim',[0,tmp(2)],'fontsize',34);
ann_x = [73,6];
ann_y = [0.15,0.3]*tmp(2);
drawnow
thebins = double(bins) / 1000;
binCt = spike_mat_bin;
h2 = subplot(2,1,2);
axis(h2)
imagesc(thebins, 1:size(binCt,1), binCt);
axis([0,12,0,600])
colormap(flipud(colormap('gray')))
set(h2,'fontsize',34);
xlabel('time(s)','fontsize',36)
ylabel('neurons','fontsize',36)
plt = [plot_fn,'clean.eps'];
print('-depsc',plt);

B = load(fn_messy);
figure
%1605 / 12 is position per 1 unit x axis
set(gcf,'position',[1118,727,1605,500]);
h1 = subplot(2,1,1)
axis(h1)
peaks1 = B.pop_burst_peak1;
peaks2 = B.pop_burst_peak2;
bins1 = bins/1000;
plot(bins/1000,butter_int_bin,'--k');
hold on
plot(bins/1000,butter_int_bin2,'k');
ylabel('x^{\rm int} (Hz/cell)','fontsize',36)
axis tight
ylim([0,20])
axis([0,12,0,20])
tmp = get(gca,'ylim');
set(h1,'ylim',[0,tmp(2)],'fontsize',34);
axis([0,12,0,20])
ann_x = [73,6];
ann_y = [0.15,0.3]*tmp(2);
drawnow
thebins = double(bins) / 1000;
binCt = spike_mat_bin;
h2 = subplot(2,1,2);
imagesc(thebins, 1:size(binCt,1), binCt);
axis([0,20,0,600])
colormap(flipud(colormap('gray')))
set(h2,'fontsize',34);
xlabel('time(s)','fontsize',36)
ylabel('neurons','fontsize',36)
plt = [plot_fn,'messy.eps'];
print('-depsc',plt);


