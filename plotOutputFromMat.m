% plotOutputFromMat
% 
% makes some plots of my model output

clear all
close all

%% parameters
printit=0;
%fn = 'avg_synapses_dt_1e-3.mat';
fn = 'testER_1_d.mat';
outfn = 'postprocessing/run_both_synapses';
firstStep = 25000;
finalStep = 150000;
minISIstep = 1;
binWidth = 10e-3; % in s
numNeuron = 300;
numEqnsPerNeuron = 7;
vThresh = -0.0;
dt = 1e-4;
save_full = 0;

%% run analysis
V = load(fn);
%A_py = A_py.voltages';
%A = A.Y';

vTypes=V.vTypes;

% compare early stages
T = firstStep:finalStep;
T = T*dt;

%A = A(firstStep:finalStep,:);
%idx = 1:1:30*7;
%A = A(:,idx);

%% plot all variables
% figure
% plot(A, 'linewidth', 1);
% title('python model')
% ylim([-1,2])

%% plot just voltages
if save_full
    Vidx = 1:numEqnsPerNeuron:numEqnsPerNeuron*numNeuron;
else
    Vidx = 1:numNeuron;
end
V = V.Y(Vidx,firstStep:finalStep)';
%clear A
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
% figure
% bar(bins(1:end-1), pop, 'hist')
% title('raw spike counts')
%% bin counts smoothed with a Butterworth filter
filtOut = filterSpikes(pop, binWidth);
figure
bar(bins(1:end-1), filtOut/(binWidth*numNeuron), 'hist')
title('smoothed firing rates')
ylabel('avg firing rate per neuron (Hz)')
xlabel('time')
axis tight

%% psd
[Pxx,F]=periodogram(zscore(filtOut), [], 512, 1/binWidth);
[Y,I] = max(Pxx);
peakFreq = F(I);
fprintf('peak frequency: %1.2fHz (%1.2fs burst period)\n', peakFreq, ...
        1/peakFreq);