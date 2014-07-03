function [binCts, bins] = binSpikes(spikeT, T, binWidth, dt)
    if iscell(spikeT)
        % If we get a cell array, assume it's a collection of
        % spikeTs
        numNeur = length(spikeT);
        numBins = ceil(length(T)*dt/binWidth)-1;
        binCts = zeros(numNeur, numBins);
        for neur = 1:numNeur
            % recursion!
            [binCts(neur, :), bins] = binSpikes( spikeT{neur}, T, ...
                                                 binWidth, dt );
        end
    else
        binWidthN = round(binWidth/dt);
        binEdges = 1:binWidthN:length(T);
        binEdgeT = T(binEdges);
        binCts = zeros(1, length(binEdgeT)-1);
        for bin = 1:length(binEdgeT)-1
            binCts(bin) = sum( spikeT >= binEdgeT(bin) & ...
                               spikeT < binEdgeT(bin+1) );
        end
        bins = binEdgeT;
    end
