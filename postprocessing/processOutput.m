fn = '../model/output/er_300_2_deg_0.00_b.mat';
units_ms = 1;


load(fn)
firstStep = 1;
finalStep = length(Y);
if units_ms
    minISIstep = 10;
    binWidth = 10; % in ms
else
    minISIstep = 1
    binWidth = 10e-3; % in s
end
numNeuron = size(vTypes,1);
numEqnsPerNeuron = 5;
vThresh = -0.0;
save_full = 0;

%% run analysis
V = load(fn);
%A_py = A_py.voltages';
%A = A.Y';

%vTypes=V.vTypes;

% compare early stages
T = firstStep:finalStep;
T = T*double(dt);

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
V = Y(Vidx,firstStep:finalStep)';
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
if units_ms
    filtOut = filterSpikes(pop, binWidth*10^-3);
else
    filtOut = filterSpikes(pop, binWidth);
end

figure
if units_ms
    bar(bins(1:end-1), filtOut/(binWidth*10^-3*numNeuron), 'hist')
else
    bar(bins(1:end-1), filtOut/(binWidth*numNeuron), 'hist')
end
title('smoothed firing rates')
ylabel('avg firing rate per neuron (Hz)')
xlabel('time')
axis tight

%% psd
if units_ms
    [Pxx,F]=periodogram(zscore(filtOut), [], 512, 10^3/binWidth);
else
    [Pxx,F]=periodogram(zscore(filtOut), [], 512, 1/binWidth);
end
[Y,I] = max(Pxx);
peakFreq = F(I);
fprintf('peak frequency %1.2fHz; %1.2fs burst period\n', peakFreq, ...
        1/peakFreq);

[Y, I] = sort(vTypes);
binCtSort = binCt(I,:);
figure;
spy(binCtSort)
