% function: spikeTimes

function spikeTime=spikeTimes(V, T, thresh, minISIstep, varargin)
    if ismatrix(V) && ~isvector(V)
        % If we get a matrix, assume voltages of diff neurons
        % (num time steps) x (num neuron)
        [numTime, numNeur] = size(V);
        if numNeur > numTime
            error(['number of neurons > number of time steps; do you ' ...
                   'need to transpose your voltage?'])
        end
        spikeTime = cell(numNeur, 1);
        for neur = 1:numNeur
            % recursion!
            spikeTime{neur} = spikeTimes(V(:, neur), T, thresh, ...
                                         minISIstep, neur);
        end
    else
        % We got a vector, so find the peaks
        s = warning('error', ['signal:findpeaks:' ...
                            'largeMinPeakHeight']);
        warning('error', 'signal:findpeaks:largeMinPeakHeight');
        try
            [spikeV, spikeLoc] = findpeaks(V, ...
                                           'minpeakheight', thresh, ...
                                           'minpeakdistance', ...
                                           minISIstep);
        catch
            fprintf('no spikes found for neuron id %d\n', ...
                    varargin{1});
            spikeLoc = [];
        end
        warning(s);
        % Return the spike times
        spikeTime = T(spikeLoc);
    end
        