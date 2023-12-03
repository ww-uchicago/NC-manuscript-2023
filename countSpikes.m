function out = countSpikes(SortedDSData,varargin)
%computes spike counts and means for the repetitions of each direction

%preallocate spikeCount matrix
SortedDSData.spikeCounts = zeros(length(SortedDSData.directions),...
    length(SortedDSData.spikeTimes{1}));

for d = 1:length(SortedDSData.directions)
    for r = 1:length(SortedDSData.spikeTimes{d})
        counts(r) = length(SortedDSData.spikeTimes{d}{r});
    end
    SortedDSData.spikeCounts(d,:) = counts;

end
SortedDSData.meanSpikeCounts = mean(SortedDSData.spikeCounts,2);
out = SortedDSData;
return

