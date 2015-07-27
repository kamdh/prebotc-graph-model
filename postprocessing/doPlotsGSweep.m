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
fontsz=24;

dopartics = 1;
docombined = 1;
do_k_vs_pI=1;
do_gE_vs_gI=1;
%partics = {'er_n300_k6.0_deg_pI0.20_rep1'};
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
           'er_n300_k6.0_deg_pI0.30_rep1', ...
           'er_n300_k6.0_deg_pI0.40_rep1', ...
           'er_n300_k6.0_deg_pI0.50_rep1', ...
           'er_n300_k6.0_deg_pI0.60_rep1', ...
           'er_n300_k6.0_deg_pI0.70_rep1', ...
           'er_n300_k6.0_deg_pI0.80_rep1', ...
           'er_n300_k6.0_deg_pI0.90_rep1', ...
           'er_n300_k6.0_deg_pI1.00_rep1' 
          };
projName = 'random_extend';
opThresh = 0.2;

%% start processing
dataDir = [getenv('HOME') '/work/prebotc/data/', projName]
fn = [dataDir, '/post/collected.mat'];
plotDir = [getenv('HOME') '/work/prebotc/data/', projName, ...
           '/plots'];


load(fn)
numgE = length(gEs);
numgI = length(gIs);
num_k = length(ks);
num_pI= length(pIs);
if docombined
    if do_k_vs_pI
    for gEidx = 1:numgE
        for gIidx = 1:numgI
            x_axis_label='connectivity k_{avg}';
            y_axis_label='inhibitory fraction p_I';
            gE = gEs(gEidx);
            gI = gIs(gIidx);
            
            pltGStr = sprintf('gE_%1.1f_gI_%1.1f', gE, gI);

            figure
            myPcolor(X,Y, chiArray(:,:,gEidx, gIidx), 'clim', [0,1]);
            titlestr=sprintf('Synchrony \\chi\ng_E = %1.1f, g_I = %1.1f', ...
                             gE,gI);
            title(titlestr, 'fontsize', fontsz)
            xlabel(x_axis_label, 'fontsize', fontsz)
            ylabel(y_axis_label,'fontsize', fontsz)
            %axis([0,1])
            % colorbar
            %colormap('gray')
            plt = [plotDir, '/', pltGStr, '_chi.eps']
            print('-depsc', plt)

            chi_mask=chiArray(:,:,gEidx,gIidx) > chi_threshold;
            
            
            figure
            myPcolor(X,Y, chiArray_std(:,:,gEidx, gIidx));
            titlestr=sprintf('SD of \\chi\ng_E = %1.1f, g_I = %1.1f', ...
                             gE,gI);
            title(titlestr, 'fontsize', fontsz)
            xlabel(x_axis_label, 'fontsize', fontsz)
            ylabel(y_axis_label,'fontsize', fontsz)
            %axis([0,1])
            % colorbar
            %colormap('gray')
            plt = [plotDir, '/', pltGStr, '_chi_std.eps']
            print('-depsc', plt)
            

            figure
            myPcolor(X,Y, amplitude_irregularity(:,:,gEidx, gIidx));
            titlestr=sprintf(['Amplitude irregularity\ng_E = %1.1f, ' ...
                              'g_I = %1.1f'],gE,gI);
            title(titlestr, 'fontsize', fontsz)
            xlabel(x_axis_label, 'fontsize', fontsz)
            ylabel(y_axis_label,'fontsize', fontsz)
            %axis([0,1])
            % colorbar
            %colormap('gray')
            plt = [plotDir, '/', pltGStr, '_amp_irreg.eps']
            print('-depsc', plt)

            figure
            myPcolor(X,Y, ibi_irregularity(:,:,gEidx, gIidx));
            titlestr=sprintf(['IBI irregularity\ng_E = %1.1f, ' ...
                              'g_I = %1.1f'],gE,gI);
            title(titlestr, 'fontsize', fontsz)
            xlabel(x_axis_label, 'fontsize', fontsz)
            ylabel(y_axis_label,'fontsize', fontsz)
            %axis([0,1])
            % colorbar
            %colormap('gray')
            plt = [plotDir, '/', pltGStr, '_ibi_irreg.eps']
            print('-depsc', plt)
            
            % figure
            % myPcolor(X,Y, dutyCycle(:,:,gEidx, gIidx))
            % titlestr = sprintf('duty cycle\ng_E = %1.1f, g_I = %1.1f', gE,gI);
            % title(titlestr, 'fontsize', fontsz)
            % xlabel(x_axis_label, 'fontsize', fontsz)
            % ylabel(y_axis_label,'fontsize', fontsz)
            % % colorbar
            % %colormap('gray')
            % plt = [plotDir, '/', pltGStr, '_duty_cycle.eps']
            % print('-depsc', plt)

            figure
            myPcolor(X,Y, fMax(:,:,gEidx, gIidx))
            titlestr=sprintf('Peak frequency (Hz)\ng_E = %1.1f, g_I = %1.1f',...
                             gE,gI);
            title(titlestr, 'fontsize', fontsz)
            title('peak frequency (1/s)','fontsize', fontsz)
            xlabel(x_axis_label, 'fontsize', fontsz)
            ylabel(y_axis_label,'fontsize', fontsz)
            % colorbar
            %colormap('gray')
            plt = [plotDir, '/', pltGStr, '_peak_freq.eps']
            print('-depsc', plt)

            figure
            myPcolor(X,Y, lag(:,:,gEidx, gIidx))
            titlestr=sprintf('Period (s)\ng_E = %1.1f, g_I = %1.1f',...
                             gE,gI);
            title(titlestr, 'fontsize', fontsz)
            xlabel(x_axis_label, 'fontsize', fontsz)
            ylabel(y_axis_label,'fontsize', fontsz)
            % colorbar
            %colormap('gray')
            plt = [plotDir, '/', pltGStr, '_lag.eps']
            print('-depsc', plt)

            % figure
            % myPcolor(X,Y, muB(:,:,gEidx, gIidx) / 1000)
            % titlestr=sprintf(['mean burst duration (s)\ng_E = %1.1f, ' ...
            %                   'g_I = %1.1f'], ...
            %                  gE,gI);
            % title(titlestr, 'fontsize', fontsz)
            % title('mean burst duration (s)','fontsize', fontsz)
            % xlabel(x_axis_label,'fontsize', fontsz)
            % ylabel(y_axis_label,'fontsize', fontsz)
            % % colorbar
            % %colormap('gray')
            % plt = [plotDir, '/', pltGStr, '_mean_burst.eps']
            % print('-depsc', plt)

            % figure
            % myPcolor(X,Y, muIBI(:,:,gEidx, gIidx) / 1000)
            % titlestr=sprintf('mean IBI (s)\ng_E = %1.1f, g_I = %1.1f',...
            %                  gE,gI);
            % title(titlestr, 'fontsize', fontsz)
            % xlabel(x_axis_label,'fontsize', fontsz)
            % ylabel(y_axis_label,'fontsize', fontsz)
            % % colorbar
            % %colormap('gray')
            % plt = [plotDir, '/', pltGStr, '_mean_IBI.eps']
            % print('-depsc', plt)

            % figure
            % myPcolor(X,Y, cvIBI(:,:,gEidx, gIidx))
            % titlestr=sprintf('CV of IBIs\ng_E = %1.1f, g_I = %1.1f',...
            %                  gE,gI);
            % title(titlestr, 'fontsize', fontsz)
            % xlabel(x_axis_label,'fontsize', fontsz)
            % ylabel(y_axis_label,'fontsize', fontsz)
            % % colorbar
            % %colormap('gray')
            % plt = [plotDir, '/', pltGStr, '_cv_IBIs.eps']
            % print('-depsc', plt)

            % figure
            % myPcolor(X,Y, cvB(:,:,gEidx, gIidx))
            % titlestr=sprintf('CV burst lengths\ng_E = %1.1f, g_I = %1.1f',...
            %                  gE,gI);
            % title(titlestr, 'fontsize', fontsz)
            % xlabel(x_axis_label,'fontsize', fontsz)
            % ylabel(y_axis_label,'fontsize', fontsz)
            % % colorbar
            % %colormap('gray')
            % plt = [plotDir, '/', pltGStr, '_cv_bursts.eps']
            % print('-depsc', plt)

            figure
            myPcolor(X,Y, op_angle_mean(:,:,gEidx, gIidx))
            titlestr=sprintf('Mean phase arg(z)\ng_E = %1.1f, g_I = %1.1f',...
                             gE,gI);
            title(titlestr, 'fontsize', fontsz)
            xlabel(x_axis_label,'fontsize', fontsz)
            ylabel(y_axis_label,'fontsize', fontsz)
            % colorbar
            %colormap('gray')
            plt = [plotDir, '/', pltGStr, '_op_angle_mean.eps']
            print('-depsc', plt)

            figure
            myPcolor(X,Y, op_angle_std(:,:,gEidx, gIidx))
            titlestr=sprintf(['SD of phase arg(z)\ng_E = ' ...
                              '%1.1f, g_I = %1.1f'],...
                             gE,gI);
            title(titlestr, 'fontsize', fontsz)
            xlabel(x_axis_label,'fontsize', fontsz)
            ylabel(y_axis_label,'fontsize', fontsz)
            % colorbar
            %colormap('gray')
            plt = [plotDir, '/', pltGStr, '_op_angle_std.eps']
            print('-depsc', plt)

            figure
            tmp=num_expir(:,:,gEidx, gIidx)./num_neurons*100;
            tmp(~chi_mask)=nan;
            myPcolor(X,Y, tmp)
            titlestr=sprintf(['Percent expiratory\ng_E = ' ...
                              '%1.1f, g_I = %1.1f'],...
                             gE,gI);
            title(titlestr, 'fontsize', fontsz)
            xlabel(x_axis_label,'fontsize', fontsz)
            ylabel(y_axis_label,'fontsize', fontsz)
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
            % title(titlestr, 'fontsize', fontsz)
            % xlabel(x_axis_label,'fontsize', fontsz)
            % ylabel(y_axis_label,'fontsize', fontsz)
            % colorbar off
            % %colormap('gray')
            % plt = [plotDir, '/', pltGStr, '_expir_times_chi.eps']
            % print('-depsc', plt)
            close all
        end % gI
    end % gE
    end
    
    if do_gE_vs_gI
    %% plot gE vs gI
    for k_idx=1:num_k
        for pI_idx=1:num_pI
            k=ks(k_idx);
            pI=pIs(pI_idx);
            plt_str = sprintf('k_%1.1f_pI_%1.1f', k,pI);
            x_axis_label='g_E (nS)';
            y_axis_label='g_I (nS)';

            figure
            myPcolor(Xg,Yg, squeeze(chiArray(k_idx, pI_idx, :,:)),...
                     'clim', [0,1])
            titlestr=sprintf(['Synchrony \\chi\nk_{avg} = %1.1f, p_I ' ...
                              '= %1.1f'], k,pI);
            title(titlestr, 'fontsize', fontsz)
            xlabel(x_axis_label, 'fontsize', fontsz)
            ylabel(y_axis_label,'fontsize', fontsz)
            %axis([0,1])
            % colorbar
            %colormap('gray')
            plt = [plotDir, '/', plt_str, '_chi.eps']
            print('-depsc', plt)

            chi_mask=squeeze(chiArray(k_idx, pI_idx, :,:) > chi_threshold);

            figure
            myPcolor(Xg,Yg, squeeze(chiArray_std(k_idx, pI_idx, :,:)))
            titlestr=sprintf(['SD of \\chi\nk_{avg} = %1.1f, p_I ' ...
                              '= %1.1f'], k,pI);
            title(titlestr, 'fontsize', fontsz)
            xlabel(x_axis_label, 'fontsize', fontsz)
            ylabel(y_axis_label,'fontsize', fontsz)
            %axis([0,1])
            % colorbar
            %colormap('gray')
            plt = [plotDir, '/', plt_str, '_chi_std.eps']
            print('-depsc', plt)
        
            figure
            myPcolor(Xg,Yg,squeeze(amplitude_irregularity(k_idx,pI_idx,:,:)), ...
                     'clim', [0,1]);
            titlestr=sprintf(['Amplitude irregularity\nk_{avg} = %1.1f, ' ...
                              'p_I = %1.1f'],k,pI);
            title(titlestr, 'fontsize', fontsz)
            xlabel(x_axis_label, 'fontsize', fontsz)
            ylabel(y_axis_label,'fontsize', fontsz)
            %axis([0,1])
            % colorbar
            %colormap('gray')
            plt = [plotDir, '/', plt_str, '_amp_irreg.eps']
            print('-depsc', plt)

            figure
            myPcolor(Xg,Yg,squeeze(ibi_irregularity(k_idx,pI_idx,:,:)), ...
                     'clim', [0,1]);
            titlestr=sprintf(['IBI irregularity\nk_{avg} = %1.1f, ' ...
                              'p_I = %1.1f'],k,pI);
            title(titlestr, 'fontsize', fontsz)
            xlabel(x_axis_label, 'fontsize', fontsz)
            ylabel(y_axis_label,'fontsize', fontsz)
            %axis([0,1])
            % colorbar
            %colormap('gray')
            plt = [plotDir, '/', plt_str, '_ibi_irreg.eps']
            print('-depsc', plt)
        end
    end
    end
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
                chival=B.chi;
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
                
                
                set(0,'DefaultAxesFontSize', fontsz)
                
                figure
                set(gcf, 'position', [1118, 727, 2675, 500])
                %% subplot 1
                h1 = subplot(2,1,1);
                axis(h1)
                plot(thebins, butterIntBin, 'k-')
                hold on
                plot(thebins(peaklocs), butterIntBin(peaklocs), 'ko')
                %ylabel('\int preBot (Hz/neuron)', 'fontsize', 20)
                ylabel('x^{\rm int} (Hz/cell)','fontsize',28)
                axis tight
                tmp = get(gca, 'ylim');
                tmp(2)=1.05*tmp(2); % add 5 pct padding to ylim
                set(h1, 'ylim', [0, tmp(2)]);
                str=sprintf('\\chi = %1.3f',chival);
                ann_x=[73,6];
                ann_y=[0.15,0.3]*tmp(2);
                drawnow
                rectangle('position', ...
                          [ann_x(1)-.15,ann_y(1)-0.02*tmp(2),...
                           ann_x(2),ann_y(2)],...
                          'facecolor', 'w')
                text(ann_x(1),ann_y(1), str, ...
                     'fontsize', 28, 'verticalalignment', 'bottom');
                %% subplot 2
                h2 = subplot(2,1,2);
                axis(h2)
                imagesc(thebins, 1:size(binCt,1), binCt(I,:));
                xlabel('time (s)', 'fontsize', 28)
                colormap(flipud(colormap('gray')))
                % set(h2, 'xtick', get(h1, 'xtick'))
                % set(h2, 'xticklabel', get(h1, 'xticklabel'))
                %spy(binCt(I,:), 'k', 1)
                %xlabel(['time bin (' num2str(binWidth*1000) ' ms)'])
                ylabel('neurons', 'fontsize', 28)
                %% finalize and print
                set(h2, 'xtick', 10:10:70)
                set(h2, 'xticklabel', get(h1,'xticklabel'))
                % ann_x=[73,6];
                % ann_y=[0.3,0.2]*max(buttrIntBin);
                % [ann_nf_x,ann_nf_y]=ds2nfu(h1,ann_x,ann_y);
                % annotation('textbox', ...
                %            [ann_nf_x(1),ann_nf_y(1),...
                %             ann_nf_x(2),ann_nf_y(2)], ...
                %            'String', str, ...
                %            'fontsize', fontsz, ...
                %            'backgroundcolor', 'w',...
                %            'fitboxtotext', 'off');
                %plt = [plotDir, '/ts_raster.eps']
                %print('-depsc', plt)
                plt = [plotDir, '/', pltGStr, 'ts_combined_' partic '.eps']
                print('-deps', plt)
                %pause                
                
                %% order parameter plot
                set(0,'DefaultAxesFontSize', 18)
                figure
                plot(opAngle/pi, opAbs, 'k.')
                xlabel('angle (\pi radians)', 'fontsize', fontsz)
                ylabel('magnitude', 'fontsize', fontsz)
                hline(opThresh, 'k--')
                %y=get(gca,'ylim');
                hold on;
                plot([-.5,-.5],[opThresh,1],'k--');
                plot([.5,.5],[opThresh,1],'k--');
                ylim([0, 1])
                xlim([-1, 1])
                title('Phase-locking variable','fontsize',fontsz)
                plt = [plotDir, '/', pltGStr, 'op_' partic '.eps']
                print('-deps', plt)

                
                close all
                
            end
        end
    end
end
