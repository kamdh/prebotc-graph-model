% compareOutputs
% debugging script to compare Tatiana's C++ code with my Python

clear all
close all

%% parameters
printit=0;

% fn1 = 'both_synapses.mat';
% fn2 = '../dashevskiy/both_synapses.txt';
% fn1 = 'orig_synapse.mat';
% fn2 = '../dashevskiy/orig_synapse.txt';
% fn1 = 'test_avg.mat';
% fn2 = '../dashevskiy/test_avg.txt';
% fn1 = 'avg_synapses.mat';
%fn2 = '../dashevskiy/avg_synapses_dt_1e-4.txt';
fn1 = 'avg_synapses_dt_5e-5.mat';
fn2 = '../dashevskiy/avg_synapses_dt_5e-5.txt';
outfn = 'postprocessing/run_both_synapses';
finalStep = 20000;
minISIstep = 2;
binWidth = 30e-3; % in s
numNeuron = 30;
numEqnsPerNeuron = 7;
vThresh = 0.0;
dt = 5e-5;

%% run analysis
A_py = load(fn1);
%A_py = A_py.voltages';
A_py = A_py.Y';
A_cpp = load(fn2);
A_cpp = A_cpp(2:end, 2:end);


% compare early stages
T = 1:finalStep;
T = T*dt;

A_cpp = A_cpp(1:finalStep,:);
idx = 1:30*10;
idx2 = [8:10:30*10, 9:10:30*10, 10:10:30*10];
idx3 = [9:10:30*10, 10:10:30*10];
idx(idx2) = [];
A_cpp = [ A_cpp(:,idx), A_cpp(:, idx3) ];

A_py = A_py(1:finalStep,:);
%idx = 1:1:30*7;
%A_py = A_py(:,idx);


%% plot all variables
subplot(2,1,1);
plot(A_py, 'linewidth', 1);
title('python model')
ylim([-1,2])
subplot(2,1,2);
plot(A_cpp, 'linewidth', 1)
title('C++ model')
ylim([-1,2])

%% plot relative errors
err_rel = abs((A_py - A_cpp)./A_cpp);
figure
imagesc(err_rel)
colorbar
ylabel('time')
xlabel('state variable')
title('relative error of python model')

%% plot absolute errors
err_abs = abs((A_py - A_cpp));
figure
imagesc(err_abs)
colorbar
ylabel('time')
xlabel('state variable')
title('absolute error of python model')
disp(['max absolute error: ' num2str(max(err_abs(:)))])
disp(['sum of abs error * dt: ' num2str(dt*sum(err_abs(:)))])

if printit
    figure(1)
    print('-depsc', [outfn '_states.eps'])
    figure(2)
    print('-depsc', [outfn '_err_rel.eps'])
    figure(3)
    print('-depsc', [outfn '_err_abs.eps'])
end

%% plot just voltages
Vidx = 1:numEqnsPerNeuron:numEqnsPerNeuron*numNeuron;
V_py = A_py(:, Vidx);
V_cpp = A_cpp(:, Vidx);
figure
subplot(2,1,1)
plot(T, V_py, '-');
ylim([-.08, 0.03])
title('python')
subplot(2,1,2)
plot(T, V_cpp, '-')
ylim([-.08, 0.03])
title('C++')

%% count spikes and bin
spikeT_py = spikeFilt(V_py, T, vThresh, minISIstep);
spikeT_cpp = spikeFilt(V_cpp, T, vThresh, minISIstep);
[binCt_py, bins] = binSpikes(spikeT_py, T, binWidth, dt);
[binCt_cpp, bins] = binSpikes(spikeT_cpp, T, binWidth, dt);
% sum over all neurons
pop_py = sum(binCt_py, 1);
pop_cpp = sum(binCt_cpp, 1);
%% plot the raw bin counts
figure
subplot(2,1,1)
bar(bins(1:end-1), pop_py, 'hist')
title('python')
subplot(2,1,2)
bar(bins(1:end-1), pop_cpp, 'hist')
title('C++')
%% bin counts smoothed with a Butterworth filter
filt_py = filterSpikes(pop_py, binWidth);
filt_cpp = filterSpikes(pop_cpp, binWidth);
figure
subplot(2,1,1)
bar(bins(1:end-1), filt_py, 'hist')
title('py')
axis tight
subplot(2,1,2)
bar(bins(1:end-1), filt_cpp, 'hist')
title('cpp')
axis tight