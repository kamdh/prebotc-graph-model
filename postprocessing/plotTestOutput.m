clear all
close all

load('../model/output/testoutput.mat')


T = dt*(1e-3)*(1:size(Y,2));

figure

subplot(2,2,1);
plot(T, Y(1,:), 'k-', 'linewidth', .5)
ylim([-60 10])
title('CS')

subplot(2,2,2);
%plot(Y(7+1,:))
plot(T, Y(2,:), 'k-', 'linewidth', .5)
ylim([-60 10])
title('CI')

subplot(2,2,3);
%plot(Y(2*7+1,:))
plot(T, Y(3,:), 'k-', 'linewidth', .5)
ylim([-60 10])
title('TS')

subplot(2,2,4);
%plot(Y(3*7+1,:))
plot(T, Y(4,:), 'k-', 'linewidth', .5)
ylim([-60 10])
title('Q')

set(gcf, 'position', [0 0 1400 1200])
print('-deps', 'output/voltages_no_network.eps')