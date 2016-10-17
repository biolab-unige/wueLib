%subjects = {'wue04' 'wue06' 'wue07' 'wue09' 'wue10' 'wue11'};
subjects = {'wue11'};

dataFolder = '/media/lgabri/My Passport/gait-2016';

for subjIdx = 1:numel(subjects)
	% read the vmrk file only
	vmrkFilenames = dir(fullfile(dataFolder,subjects{subjIdx},'*vmrk'));
	% read the events from file
	eveFilename = (fullfile(dataFolder,subjects{subjIdx},strcat(subjects{subjIdx},'_walking_off_events.csv')));
	eventsData = importdata(eveFilename);

	% unwrap event data structure
	[eventGroupNames,~,origOrder] = unique(eventsData.colheaders);
	nEventGroups = numel(eventGroupNames);
	eventGroups = struct('label',repmat([],nEventGroups,1),...
												'type',repmat('Gait',nEventGroups,1),...
												'epochs',repmat([],nEventGroups,1),...
												'samples',repmat([],nEventGroups,1),...
												'length',repmat([],nEventGroups,1),...
												'chan_num',repmat([],nEventGroups,1));

	for eventGroupIdx = 1:nEventGroups
			eventGroups(eventGroupIdx).label = eventGroupNames{eventGroupIdx};

			% we store times (ms) in samples for sake of simple coding
			% % WARN % those below are milliseconds and not samples
			eventGroups(eventGroupIdx).samples = eventsData.data(:,strcmp(eventsData.colheaders,...
																															eventGroupNames(eventGroupIdx)));

	end


	% check whether we have as many trials (ie files) as rows in table (ie trials)
	if size(eventsData.data,1) < numel(vmrkFilenames)
			error(' Not enought events/files');
	end

	for trialIdx = 1:numel(vmrkFilenames)
			% build paths
			vmrkCurrFname = fullfile(dataFolder,subjects{subjIdx},vmrkFilenames(trialIdx).name);
			vhdrCurrFname = strrep(vmrkCurrFname,'vmrk','vhdr');
 
			% read hdr information
			eegHdr = read_brainvision_vhdr(vhdrCurrFname);

			% read events from old vmrk file
			eegEvents = read_brainvision_vmrk(vmrkCurrFname);

			% extract trial-specific events from eventGroups
			newEvents = struct(eegEvents);
			for eventGroupIdx = 1:nEventGroups
					newEvents(eventGroupIdx).label = eventGroups(eventGroupIdx).label;
					newEvents(eventGroupIdx).samples = round(eventGroups(eventGroupIdx).samples(trialIdx,:).*eegHdr.Fs);
					nEventInTrial = numel(eventGroups(eventGroupIdx).samples(trialIdx,:));
					newEvents(eventGroupIdx).type = repmat({'Gait'},nEventInTrial,1);
					newEvents(eventGroupIdx).epochs = ones(nEventInTrial,1);
					newEvents(eventGroupIdx).length = ones(nEventInTrial,1);
					newEvents(eventGroupIdx).chan_num = zeros(nEventInTrial,1);
			end

			% append to original events
			eegEvents = [eegEvents newEvents];

			% write vmrk back to disk
			[writeOutPath,~] = fileparts(vmrkCurrFname);
			write_brainvision_vmrk(writeOutPath,eegHdr,eegEvents);
	end
		
	
end
