function eventStruct = bn_readExternalEventFile(filename)
%BN_READEXTERNALEVENTFILE reads and external event file with comma separated values where 
%	number of rows represents single events and number of colmuns representes event types
%	EVENTSTRUCT = BN_READEXTERNALEVENTFILE(FILENAME) Long description
%

% Edited 2015-05-13 by Gabriele Arnulfo <gabriele.arnulfo@gmail.com>

eventsRaw = importdata(filename); %,'delimiter',',');

[nEvents nEventTypes] = size(eventsRaw.data);

% the first column represents the event indices 
% and should not be counted as an event class
nEventTypes = nEventTypes;

for evTypeIdx = 1:nEventTypes
		eventStruct(evTypeIdx).type			= repmat({'custom'},[1 nEvents]);
		eventStruct(evTypeIdx).label 		= eventsRaw.colheaders{evTypeIdx};
		eventStruct(evTypeIdx).samples	= nan(1,nEvents);
		eventStruct(evTypeIdx).times 		= eventsRaw.data(:,evTypeIdx);
		eventStruct(evTypeIdx).epochs 	= ones(1,nEvents);
		eventStruct(evTypeIdx).length 	= ones(1,nEvents);
		eventStruct(evTypeIdx).chan_num = zeros(1,nEvents);
 
end


end
