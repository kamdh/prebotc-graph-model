% function: spikeTimes

function spikeTime=spikeTimes(V, T, thresh, minISIstep)
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
            spikeTime{neur} = spikeTimes(V(:, neur), T, thresh, minISIstep);
        end
    else
        % We got a vector, so find the peaks
        [spikeV, spikeLoc] = findpeaks(V, ...
                                       'minpeakheight', thresh, ...
                                       'minpeakdistance', minISIstep);
        % Return the spike times
        spikeTime = T(spikeLoc);
    end
        