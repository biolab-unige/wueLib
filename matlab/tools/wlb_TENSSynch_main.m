function status = bn_TENSSynch_main(varargin)
%BN_TENSSYNCH_MAIN Main interface for synchrnization wue data
%	FLAG = BN_TENSSYNCH_MAIN(VARARGIN) 


% Edited 2014-09-22 by Gabriele Arnulfo <gabriele.arnulfo@gmail.com>

%%% DEFINE SUPPORTED FILE EXTENSIONS %%%%
%ext_hdeeg 				= '*.eeg';
%ext_pcs 					= '*.xml';
%ext_emg 					= '*.txt';
%event_timeoffset 	= 0;

%%% DEFINE INPUT STRUCTURE %%%
p = inputParser;
p.addParamValue('path_pcs',[],@ischar);
p.addParamValue('path_emg',[],@ischar);
p.addParamValue('path_hdeeg',[],@ischar);
p.addParamValue('path_events',[],@ischar);
p.addParamValue('events_label',[],@ischar);
p.addParamValue('outdir',[],@ischar);

p.parse(varargin{:});
recordingModalities = cell(1,3);

if ~isempty(p.Results.path_hdeeg) 
		path_hdeeg = p.Results.path_hdeeg; 
		recordingModalities{1} = 'eeg';
end
if ~isempty(p.Results.path_emg) 
		path_emg = p.Results.path_emg; 
		recordingModalities{2} = 'emg';
end
if ~isempty(p.Results.path_pcs) 
		path_pcs = p.Results.path_pcs; 
		recordingModalities{3} = 'pcs';	
end

if ~isempty(p.Results.path_events)
		path_events = p.Results.path_events;
end

if(isempty(p.Results.outdir))
		% search for a non-empty directory between recordingModalities
		fieldNames = regexp(fieldnames(p.Results),'^.*path.*','match');
		fieldNames(cellfun(@isempty,fieldNames)) = [];
		fieldNames = [fieldNames{:}];

		nonEmptyField = zeros(size(fieldNames));
		for field = 1:numel(fieldNames)
			nonEmptyField(field) = ~isempty(p.Results.(fieldNames{field}));
		end; clear field;
		[outdir, ~] = fileparts(p.Results.(fieldNames{find(nonEmptyField,1,'first')}));
else
    if(~exist(p.Results.outdir,'dir'))
        mkdir(p.Results.outdir);
    end
		outdir = p.Results.outdir;
end

%if(~isempty(p.Results.events_path))
%    event_timeoffset = 1;
%    event_offset_dir = p.Results.events_path;
%    event_marker_label = p.Results.events_labels;
%end

recordingModalities(cellfun(@isempty,recordingModalities)) = [];
recordingModalities = strjoin(recordingModalities,'_');

switch recordingModalities
		case 'eeg_emg'
				status = bn_EEGEMGSynch(path_hdeeg,path_emg);
		case 'eeg_emg_pcs'
				status = bn_EEGEMGPCSSynch(path_hdeeg,path_emg,path_pcs);
		case 'emg_pcs'
				status = bn_EMGPCSSynch(path_emg,path_pcs,outdir,path_events);
		case 'eeg_pcs'
				status = bn_EEGPCSSynch(path_hdeeg,path_pcs);
		otherwise
				status = -1;
				error('Unsupported modality');
end

end
