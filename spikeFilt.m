% function: spikeFilt

function spikeTime=spikeFilt(V, T, thresh, minISIstep)
    if ismatrix(V) && ~isvector(V)
        % If we get a matrix, assume voltages of diff neurons
        % Short dimension indexes neurons
        S = size(V);
        [numNeur, neurIdx] = min(S); % n
        [tmp, timeIdx] = max(S);
        V = reshape(V, [S(neurIdx), S(timeIdx)]);
        spikeTime = cell(numNeur, 1);
        for neur = 1:numNeur
            % recursion!
            spikeTime{neur} = spikeFilt(V(neur, :), T, thresh, minISIstep);
        end
    else
        % We got a vector, so find the peaks
        [spikeV, spikeLoc] = findpeaks(V, ...
                                       'minpeakheight', thresh, ...
                                       'minpeakdistance', minISIstep);
        % Return the spike times
        spikeTime = T(spikeLoc);
    end
    