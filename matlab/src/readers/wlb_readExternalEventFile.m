function eventStruct = wlb_readExternalEventFile(filename, trialIdx)
%BN_READEXTERNALEVENTFILE reads and external event file with comma separated values where 
%	number of rows represents single events and number of colmuns representes event types
%	EVENTSTRUCT = BN_READEXTERNALEVENTFILE(FILENAME) Long description
%

% Edited 2015-05-13 by Gabriele Arnulfo <gabriele.arnulfo@gmail.com>

		eventsRaw = importdata(filename); 

		[nEvents, nEventTypes] = size(eventsRaw.data);

		% the first column represents the event indices (i.e. trial number)
		% and should not be counted as an event class
		% nEventTypes = nEventTypes-1;
		if ~isempty(trialIdx)	
				nEvents = numel(trialIdx); 
		else
				trialIdx = 1:nEvents;
		end

		% init array of structure
		eventStruct = struct('label',eventsRaw.colheaders);

		for evTypeIdx = 1:nEventTypes

				eventStruct(evTypeIdx).type			= repmat({'custom'},[1 nEvents]);
				eventStruct(evTypeIdx).label 		= eventsRaw.colheaders{evTypeIdx};
				eventStruct(evTypeIdx).samples	= nan(1,nEvents);
				eventStruct(evTypeIdx).times 		= eventsRaw.data(trialIdx,evTypeIdx);
				eventStruct(evTypeIdx).epochs 	= ones(1,nEvents);
				eventStruct(evTypeIdx).length 	= ones(1,nEvents);
				eventStruct(evTypeIdx).chan_num = zeros(1,nEvents);
		 
		end

end
