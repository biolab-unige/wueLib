function status = wlb_TENSSynch_main(varargin)
%BNTENSSYNCHMAIN Main interface for synchrnization wue data
%	FLAG = BNTENSSYNCHMAIN(VARARGIN) 
%	USAGE:
%		
global globDebug;

% Edited 2014-09-22 by Gabriele Arnulfo <gabriele.arnulfo@gmail.com>

%%% DEFINE INPUT STRUCTURE %%%
		p = inputParser;
		p.addParamValue('pathPcs','',@ischar);
		p.addParamValue('pathEmg','',@ischar);
		p.addParamValue('pathHdeeg','',@ischar);
		p.addParamValue('pathEvents','',@ischar);
		p.addParamValue('fNameFilters',{''},@iscell);
		p.addParamValue('outdir','',@ischar);
		p.addParamValue('pcsCuttingTime',[0 0],@(x) isnumeric(x) & numel(x)==2 );
		p.addParamValue('emgCuttingTime',[0 0],@(x) isnumeric(x) & numel(x)==2 );
		p.addParamValue('eegCuttingTime',[0 0],@(x) isnumeric(x) & numel(x)==2 );
		p.addParamValue('automaticDetection',true,@islogical);
		p.addParamValue('debug',false,@islogical);
		p.addParamValue('pcsRefChannel',2,@isnumeric);
		p.addParamValue('findTensPctg',80,@isnumeric);

		p.parse(varargin{:});
		recordingModalities = cell(1,3);

		if ~isempty(p.Results.pathHdeeg) 
				pathHdeeg = p.Results.pathHdeeg; 
				recordingModalities{1} = 'eeg';
		end
		if ~isempty(p.Results.pathEmg) 
				pathEmg = p.Results.pathEmg; 
				recordingModalities{2} = 'emg';
		end
		if ~isempty(p.Results.pathPcs) 
				pathPcs = p.Results.pathPcs; 
				recordingModalities{3} = 'pcs';	
		end

		if ~isempty(p.Results.pcsCuttingTime)
				pcsCuttingTime = p.Results.pcsCuttingTime;
		end

		if ~isempty(p.Results.emgCuttingTime)
				emgCuttingTime = p.Results.emgCuttingTime;
        end

        if ~isempty(p.Results.eegCuttingTime)
                eegCuttingTime = p.Results.eegCuttingTime;
        end
        
        globDebug = p.Results.debug;

		pathEvents = p.Results.pathEvents;
		fnameFilters = p.Results.fNameFilters;
		automaticDetection = p.Results.automaticDetection;
        findTensPctg = p.Results.findTensPctg;
        pcsRefChannel = p.Results.pcsRefChannel;
        
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

		recordingModalities(cellfun(@isempty,recordingModalities)) = [];
		recordingModalities = strjoin(recordingModalities,'_');

		switch recordingModalities
				case 'eeg_emg'
						%status = wlb_EEGEMGSynch(pathHdeeg,pathEmg);
				case 'eeg_emg_pcs'
						status = wlb_EEGEMGPCSSynch(pathHdeeg,pathEmg,pathPcs,outdir,pathEvents,...
												fnameFilters,'pcsCuttingTime',pcsCuttingTime,'emgCuttingTime',emgCuttingTime,...
												'automaticDetection',automaticDetection,'findTensPctg',findTensPctg,...
                                                'pcsRefChannel',pcsRefChannel);
				case 'emg_pcs'
						status = wlb_EMGPCSSynch(pathEmg,pathPcs,outdir,pathEvents,...
												fnameFilters,'pcsCuttingTime',pcsCuttingTime,'emgCuttingTime',emgCuttingTime,...
												'automaticDetection',automaticDetection,'findTensPctg',findTensPctg,...
                                                'pcsRefChannel',pcsRefChannel);
				case 'eeg_pcs'
						status = wlb_EEGPCSSynch(pathHdeeg,pathPcs,outdir,pathEvents,...
												fnameFilters,'pcsCuttingTime',pcsCuttingTime,'eegCuttingTime',eegCuttingTime,...
												'automaticDetection',automaticDetection,'findTensPctg',findTensPctg,...
                                                'pcsRefChannel',pcsRefChannel);
				otherwise
						error('Unsupported modality');
		end
end
