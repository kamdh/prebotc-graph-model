clear all
close all

load('test.mat')

figure
subplot(2,2,1);
plot(Y(1,:))
title('CS')
subplot(2,2,2);
plot(Y(1*7+1,:))
title('CI')
subplot(2,2,3);
plot(Y(2*7+1,:))
title('TS')
subplot(2,2,4);
plot(Y(3*7+1,:))
title('Sil')
