function varargout = bn_convertPLXtoEEG(filename)
% bn_convertPLXtoEEG

% Edited 2015-04-01 by Gabriele Arnulfo <gabriele.arnulfo@gmail.com>


dataStruct = readPLXFileC( filename,'all');
channelNames = [{dataStruct.ContinuousChannels.Name}];

% Extract LFPs channels 
lfpsChannelsMask = ~cellfun(@isempty,regexp(channelNames, '^LfpChannel\d+'));

if(sum(lfpsChannelsMask) == 0)
		error('No LFPs recorded');
end

minDataLength = min(dataStruct.ContSampleCounts(lfpsChannelsMask));

data = [{dataStruct.ContinuousChannels(lfpsChannelsMask).Values}];
data = cellfun(@(x) x(1:minDataLength),data,'UniformOutput',false);
data = double([data{:}]);


data = (data .* dataStruct.ContMaxMagnitudeMV)./ ( (0.5*2^dataStruct.BitsPerContSample).*dataStruct.ContinuousChannels(end).ADGain );

[p,f,e] = fileparts(filename);
labels  = regexp(filename, '_[L|R](\d+\w+)q?_','match');
labels  = [labels{:}];
labels  = labels(3:end-1);

hdr.label = regexp(labels,'[A-Z]','match');
hdr.nChans = numel(hdr.label); 
hdr.Fs = dataStruct.ContinuousChannels(end).ADFrequency;
hdr.chanunit = repmat({'mV'},[hdr.nChans,1]);

% write brainvision support ch x T matrices
write_brainvision_eeg(fullfile(p,cleanFilenameString(f)),hdr, data');

end

function str = cleanFilenameString(str)
		str = regexprep(str,'[\s+|[.]]','_');
end


