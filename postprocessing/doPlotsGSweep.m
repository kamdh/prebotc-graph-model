clear all
close all
% for hyak
set(0,'DefaultFigureVisible','off');
% others
set(0,'DefaultFigurePaperPositionMode','auto')
set(0,'DefaultAxesFontSize', 18)
set(0,'DefaultAxesFontName','Helvetica')
set(0,'DefaultAxesLineWidth',1);

chi_threshold=0.25;
num_neurons=300;
expir_threshold=15;

dopartics = 0;
docombined = 1;
partics = {'er_n300_k1.0_deg_pI0.00_rep1',...
           'er_n300_k1.0_deg_pI0.10_rep1',...
           'er_n300_k1.0_deg_pI0.20_rep1',...
           'er_n300_k1.0_deg_pI0.50_rep1',...
           'er_n300_k1.0_deg_pI1.00_rep1', ...
           'er_n300_k3.0_deg_pI0.00_rep1', ...
           'er_n300_k3.0_deg_pI0.10_rep1', ...
           'er_n300_k3.0_deg_pI0.20_rep1', ...
           'er_n300_k3.0_deg_pI0.50_rep1', ...
           'er_n300_k3.0_deg_pI1.00_rep1', ...
           'er_n300_k2.0_deg_pI0.00_rep1', ...
           'er_n300_k2.0_deg_pI0.10_rep1', ...
           'er_n300_k2.0_deg_pI0.20_rep1', ...
           'er_n300_k2.0_deg_pI0.50_rep1', ...
           'er_n300_k2.0_deg_pI1.00_rep1', ...
           'er_n300_k6.0_deg_pI0.00_rep1', ...
           'er_n300_k6.0_deg_pI0.10_rep1', ...
           'er_n300_k6.0_deg_pI0.20_rep1', ...
           'er_n300_k6.0_deg_pI0.50_rep1', ...
           'er_n300_k6.0_deg_pI1.00_rep1' 
          };
projName = 'random_fine_g_sweep';
opThresh = 0.2;

%% start processing
dataDir = [getenv('HOME') '/work/prebotc/data/', projName]
fn = [dataDir, '/post/collected.mat'];
plotDir = [getenv('HOME') '/work/prebotc/data/', projName, ...
           '/plots'];

x_axis_label='connectivity k_{avg}';
y_axis_label='inhibitory fraction p_I';

load(fn)
numgE = length(gEs);
numgI = length(gIs);
if docombined
    for gEidx = 1:numgE
        for gIidx = 1:numgI
            gE = gEs(gEidx);
            gI = gIs(gIidx);
            
            pltGStr = sprintf('gE_%1.1f_gI_%1.1f', gE, gI);

            figure
            myPcolor(X,Y, chiArray(:,:,gEidx, gIidx));
            titlestr=sprintf('synchrony \\chi\ng_E = %1.1f, g_I = %1.1f', ...
                             gE,gI);
            title(titlestr, 'fontsize', 32)
            xlabel(x_axis_label, 'fontsize', 24)
            ylabel(y_axis_label,'fontsize', 24)
            %axis([0,1])
            % colorbar
            %colormap('gray')
            plt = [plotDir, '/', pltGStr, '_chi.eps']
            print('-depsc', plt)

            chi_mask=chiArray(:,:,gEidx,gIidx) > chi_threshold;
            
            % figure
            % myPcolor(X,Y, dutyCycle(:,:,gEidx, gIidx))
            % titlestr = sprintf('duty cycle\ng_E = %1.1f, g_I = %1.1f', gE,gI);
            % title(titlestr, 'fontsize', 32)
            % xlabel(x_axis_label, 'fontsize', 24)
            % ylabel(y_axis_label,'fontsize', 24)
            % % colorbar
            % %colormap('gray')
            % plt = [plotDir, '/', pltGStr, '_duty_cycle.eps']
            % print('-depsc', plt)

            figure
            myPcolor(X,Y, fMax(:,:,gEidx, gIidx))
            titlestr=sprintf('peak frequency (Hz)\ng_E = %1.1f, g_I = %1.1f',...
                             gE,gI);
            title(titlestr, 'fontsize', 32)
            title('peak frequency (1/s)','fontsize', 32)
            xlabel(x_axis_label, 'fontsize', 24)
            ylabel(y_axis_label,'fontsize', 24)
            % colorbar
            %colormap('gray')
            plt = [plotDir, '/', pltGStr, '_peak_freq.eps']
            print('-depsc', plt)

            figure
            myPcolor(X,Y, lag(:,:,gEidx, gIidx))
            titlestr=sprintf('dominant period (s)\ng_E = %1.1f, g_I = %1.1f',...
                             gE,gI);
            title(titlestr, 'fontsize', 32)
            xlabel(x_axis_label, 'fontsize', 24)
            ylabel(y_axis_label,'fontsize', 24)
            % colorbar
            %colormap('gray')
            plt = [plotDir, '/', pltGStr, '_lag.eps']
            print('-depsc', plt)

            % figure
            % myPcolor(X,Y, muB(:,:,gEidx, gIidx) / 1000)
            % titlestr=sprintf(['mean burst duration (s)\ng_E = %1.1f, ' ...
            %                   'g_I = %1.1f'], ...
            %                  gE,gI);
            % title(titlestr, 'fontsize', 32)
            % title('mean burst duration (s)','fontsize', 32)
            % xlabel(x_axis_label,'fontsize', 24)
            % ylabel(y_axis_label,'fontsize', 24)
            % % colorbar
            % %colormap('gray')
            % plt = [plotDir, '/', pltGStr, '_mean_burst.eps']
            % print('-depsc', plt)

            % figure
            % myPcolor(X,Y, muIBI(:,:,gEidx, gIidx) / 1000)
            % titlestr=sprintf('mean IBI (s)\ng_E = %1.1f, g_I = %1.1f',...
            %                  gE,gI);
            % title(titlestr, 'fontsize', 32)
            % xlabel(x_axis_label,'fontsize', 24)
            % ylabel(y_axis_label,'fontsize', 24)
            % % colorbar
            % %colormap('gray')
            % plt = [plotDir, '/', pltGStr, '_mean_IBI.eps']
            % print('-depsc', plt)

            % figure
            % myPcolor(X,Y, cvIBI(:,:,gEidx, gIidx))
            % titlestr=sprintf('CV of IBIs\ng_E = %1.1f, g_I = %1.1f',...
            %                  gE,gI);
            % title(titlestr, 'fontsize', 32)
            % xlabel(x_axis_label,'fontsize', 24)
            % ylabel(y_axis_label,'fontsize', 24)
            % % colorbar
            % %colormap('gray')
            % plt = [plotDir, '/', pltGStr, '_cv_IBIs.eps']
            % print('-depsc', plt)

            % figure
            % myPcolor(X,Y, cvB(:,:,gEidx, gIidx))
            % titlestr=sprintf('CV burst lengths\ng_E = %1.1f, g_I = %1.1f',...
            %                  gE,gI);
            % title(titlestr, 'fontsize', 32)
            % xlabel(x_axis_label,'fontsize', 24)
            % ylabel(y_axis_label,'fontsize', 24)
            % % colorbar
            % %colormap('gray')
            % plt = [plotDir, '/', pltGStr, '_cv_bursts.eps']
            % print('-depsc', plt)

            figure
            myPcolor(X,Y, op_angle_mean(:,:,gEidx, gIidx))
            titlestr=sprintf('mean OP phase\ng_E = %1.1f, g_I = %1.1f',...
                             gE,gI);
            title(titlestr, 'fontsize', 32)
            xlabel(x_axis_label,'fontsize', 24)
            ylabel(y_axis_label,'fontsize', 24)
            % colorbar
            %colormap('gray')
            plt = [plotDir, '/', pltGStr, '_op_angle_mean.eps']
            print('-depsc', plt)

            figure
            myPcolor(X,Y, op_angle_std(:,:,gEidx, gIidx))
            titlestr=sprintf(['std. of OP phase\ng_E = ' ...
                              '%1.1f, g_I = %1.1f'],...
                             gE,gI);
            title(titlestr, 'fontsize', 32)
            xlabel(x_axis_label,'fontsize', 24)
            ylabel(y_axis_label,'fontsize', 24)
            % colorbar
            %colormap('gray')
            plt = [plotDir, '/', pltGStr, '_op_angle_std.eps']
            print('-depsc', plt)

            figure
            tmp=num_expir(:,:,gEidx, gIidx);
            tmp(~chi_mask)=nan;
            myPcolor(X,Y, tmp)
            titlestr=sprintf(['# expiratory\ng_E = ' ...
                              '%1.1f, g_I = %1.1f'],...
                             gE,gI);
            title(titlestr, 'fontsize', 32)
            xlabel(x_axis_label,'fontsize', 24)
            ylabel(y_axis_label,'fontsize', 24)
            % colorbar
            %colormap('gray')
            plt = [plotDir, '/', pltGStr, '_num_expir.eps']
            print('-depsc', plt)
            
            expir_mask=num_expir(:,:,gEidx,gIidx) > ...
                expir_threshold;
            intersection_mask = chi_mask & expir_mask;

            % figure
            % myPcolor(X,Y, intersection_mask, flipud(gray(2)))
            % titlestr=sprintf(['synchrony & expiration\ng_E = ' ...
            %                   '%1.1f, g_I = %1.1f'],...
            %                  gE,gI);
            % title(titlestr, 'fontsize', 32)
            % xlabel(x_axis_label,'fontsize', 24)
            % ylabel(y_axis_label,'fontsize', 24)
            % colorbar off
            % %colormap('gray')
            % plt = [plotDir, '/', pltGStr, '_expir_times_chi.eps']
            % print('-depsc', plt)
            close all
        end % gI
    end % gE
end

if dopartics
    %%%%%%%%% plot some trajectories
    idx = 1;
    for partic=partics
        partic = char(partic);
        for gEidx = 1:numgE
            for gIidx = 1:numgI
                gE = gEs(gEidx);
                gI = gIs(gIidx);
                pltGStr = sprintf('gE_%1.1f_gI_%1.1f_', gE, gI);
                fileGStr = sprintf('gE%1.1f_gI%1.1f', gE, gI);
                fn1 = [dataDir, '/output/', partic, '_', fileGStr, '.mat'];
                fn2 = [dataDir, '/post/', partic, '_', fileGStr, '_post.mat'];
                
                %% old
                % A = load(fn2, 'vTypes');
                % vTypes = A.vTypes;
                % clear A
                % [Y,I] = sort(vTypes);

                B = load(fn2);
                vTypes = B.vertex_types;
                [Y,I] = sort(vTypes);
                trans = 12;
                binWidth = double( max(diff(B.bins)) ) / 1000;
                thebins = double( B.bins(trans:end) )/1000;
                butterIntBin = B.butter_int_bin(trans:end);
                binCt = B.spike_mat_bin(:,trans:end);
                ops = B.ops;
                opAngle = angle(ops);
                opAbs = abs(ops);
                peaklocs = B.pop_burst_peak;
                maxindex = length(B.bins)-trans
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
                %print('-depsc', plt)
                plt = [plotDir, '/', pltGStr, 'ts_combined_' partic '.eps']
                print('-depsc', plt)
                %% order parameter plot
                figure
                plot(opAngle/pi, opAbs, 'k.')
                xlabel('angle (\pi radians)')
                ylabel('magnitude')
                hline(opThresh, 'k--')
                ylim([0, 1])
                xlim([-1, 1])
                title('neuron order parameters')
                plt = [plotDir, '/', pltGStr, 'op_' partic '.eps']
                print('-depsc', plt)

                close all
            end
        end
    end
end
