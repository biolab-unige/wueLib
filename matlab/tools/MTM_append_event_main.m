function eegFname = MTM_append_event_main(varargin)
%BNTENSSYNCHMAIN Main interface for synchrnization wue data
%	FLAG = BNTENSSYNCHMAIN(VARARGIN) 
%	USAGE:
%		

% Edited 2014-09-22 by Gabriele Arnulfo <gabriele.arnulfo@gmail.com>

%%% DEFINE INPUT STRUCTURE %%%
		p = inputParser;
		p.addParamValue('pathHdeeg','',@ischar);
		p.addParamValue('pathEventFile','',@ischar);
		p.addParamValue('fNameFilters',{''},@iscell);
		p.addParamValue('outdir','',@ischar);
		p.addParamValue('ecg_channel','',@ischar);

		p.parse(varargin{:});
        
		fnameFilters = p.Results.fNameFilters;
        
		if(isempty(p.Results.outdir))
				% search for a non-empty directory between recordingModalities
				fieldNames = regexp(fieldnames(p.Results),'^.*path.*','match');
				fieldNames(cellfun(@isempty,fieldNames)) = [];
				fieldNames = [fieldNames{:}];

				nonEmptyField = zeros(size(fieldNames));
				for field = 1:numel(fieldNames)
					nonEmptyField(field) = ~isempty(p.Results.(fieldNames{field}));
				end; clear field;
				outdir = p.Results.(fieldNames{find(nonEmptyField,1,'first')});
		else
				if(~exist(p.Results.outdir,'dir'))
					 mkdir(p.Results.outdir);
				end
				outdir = p.Results.outdir;
		end

		eegFname = MTM_append_event(p.Results.pathHdeeg,p.Results.pathEventFile,outdir,fnameFilters,p.Results.ecg_channel);
end