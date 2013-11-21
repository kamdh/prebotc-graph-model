% plotOutputFromMat
% 
% makes some plots of my model output

clear all
close all

%% parameters
printit=0;
fn = 'avg_synapses_dt_1e-3.mat';
outfn = 'postprocessing/run_both_synapses';
firstStep = 10000;
finalStep = 49000;
minISIstep = 1;
binWidth = 15e-3; % in s
numNeuron = 30;
numEqnsPerNeuron = 7;
vThresh = -0.02;
dt = 1e-3;

%% run analysis
A = load(fn);
%A_py = A_py.voltages';
A = A.Y';

% compare early stages
T = firstStep:finalStep;
T = T*dt;

A = A(firstStep:finalStep,:);
%idx = 1:1:30*7;
%A = A(:,idx);

%% plot all variables
% figure
% plot(A, 'linewidth', 1);
% title('python model')
% ylim([-1,2])

%% plot just voltages
Vidx = 1:numEqnsPerNeuron:numEqnsPerNeuron*numNeuron;
V = A(:, Vidx);
% figure
% plot(T, V, '-');
% ylim([-.08, 0.03])
% title('python')

%% count spikes and bin
spikeT = spikeTimes(V, T, vThresh, minISIstep);
[binCt, bins] = binSpikes(spikeT, T, binWidth, dt);
% sum over all neurons
pop = sum(binCt, 1);
%% plot the raw bin counts
figure
bar(bins(1:end-1), pop, 'hist')
title('python')
%% bin counts smoothed with a Butterworth filter
filtOut = filterSpikes(pop, binWidth);
figure
bar(bins(1:end-1), filtOut, 'hist')
title('py')
axis tight
