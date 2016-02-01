function outeegfilename = MTM_append_event(varargin)
%BN_TENSSYNCH synchronize multiple files in the same vhdr/eeg file
%	FLAG = BN_TENSSYNCH(VARARGIN) I'll edit it when it will be ready
%

% Edited 2014-09-22 by Gabriele Arnulfo <gabriele.arnulfo@gmail.com>
%
%%% DEFINE SUPPORTED FILE EXTENSIONS %%%%

p = inputParser;
p.addRequired('pathHdeeg',@ischar);
p.addRequired('pathEventFile',@ischar);
p.addRequired('outdir',@ischar);

p.addOptional('fnameFilters',{''},@iscell);
p.addOptional('ecg_channel','',@ischar);

p.parse(varargin{:});

pathHdeeg = p.Results.pathHdeeg;
pathEventFile = p.Results.pathEventFile;
fnameFilters = p.Results.fnameFilters;
ecg_channel  = p.Results.ecg_channel;
% get filenames within each path
eegFileNames = dir(fullfile(pathHdeeg,'*.eeg'));
eventFileNames = dir(fullfile(pathEventFile,'*.txt'));

if ~isempty(fnameFilters{1})
    eegFileNames = filterFnames(eegFileNames,fnameFilters);
    eventFileNames = filterFnames(eventFileNames,fnameFilters);
end
    
for fileIdx =1:numel(eegFileNames)
    
    eegFname = fullfile(pathHdeeg,eegFileNames(fileIdx).name);
    eventFname = fullfile(pathEventFile,eventFileNames(fileIdx).name);
    [eegHdr, eegData,vmrk_event] = wlb_readBrainvision( eegFname );  
    outeegfilename{fileIdx} = fullfile(p.Results.outdir,eegFileNames(fileIdx).name);
    write_brainvision_eeg(p.Results.outdir,eegHdr,eegData);
%     ecg_channel = find(strcmp(eegHdr.label,ecg_channel));
%     eegHdr.chantype(end-1:end)={'lfp','lfp'};
%     eegHdr.chantype(find(strcmp(eegHdr.label,'none'))) = {'Other'};
%     eegHdr.chantype(ecg_channel)={'ecg'};
%     write_brainvision_vhdr(p.Results.outdir,eegHdr);
%     event = importdata(eventFname);
%     move_onset_event = vmrk_event(2).samples;
% %     move_onset_event = vmrk_event(2).samples(find(strcmp(vmrk_event(2).type,'move_on')));
% 
%     visual_event = move_onset_event-round(event.data(:,4)'*eegHdr.Fs/1000);
%     move_offset_event = move_onset_event + round(event.data(:,5)'*eegHdr.Fs/1000);
%     return_event = move_onset_event + round(event.data(:,30)'*eegHdr.Fs);
%     vmrk_event(3:5) = vmrk_event(2);
%     vmrk_event(2).label = 'visual_tgt';
%     vmrk_event(3).label = 'move_on';
%     vmrk_event(4).label = 'move_off';
%     vmrk_event(5).label = 'return';
%     vmrk_event(2).type = repmat({'Response_vistgt'},length(visual_event),1);
%     vmrk_event(3).type = repmat({'Response_moveon'},length(move_onset_event),1);
%     vmrk_event(4).type = repmat({'Response_moveoff'},length(move_offset_event),1);
%     vmrk_event(5).type = repmat({'Response_return'},length(return_event),1);
%     vmrk_event(2).samples = (visual_event);
%     vmrk_event(3).samples = (move_onset_event);
%     vmrk_event(4).samples = (move_offset_event);
%     vmrk_event(5).samples = (return_event);
%     vmrk_event(2).epochs = ones(1,length(visual_event));
%     vmrk_event(3).epochs = ones(1,length(move_onset_event));
%     vmrk_event(4).epochs = ones(1,length(move_offset_event));
%     vmrk_event(5).epochs = ones(1,length(return_event));
%     vmrk_event(2).length = ones(1,length(visual_event));
%     vmrk_event(3).length = ones(1,length(move_onset_event));
%     vmrk_event(4).length = ones(1,length(move_offset_event));
%     vmrk_event(5).length = ones(1,length(return_event));
%     vmrk_event(2).chan_num = zeros(1,length(visual_event));
%     vmrk_event(3).chan_num = zeros(1,length(move_onset_event));
%     vmrk_event(4).chan_num = zeros(1,length(move_offset_event));
%     vmrk_event(5).chan_num = zeros(1,length(return_event));
% 
%     write_brainvision_vmrk(p.Results.outdir,eegHdr,vmrk_event);

end
end

function fnames = filterFnames(fnames,pattern)
%FILTERFNAMES Description
%	FNAME = FILTERFNAMES(FNAMES,PATTERN) Long description
%
		tmp = {fnames.name};
		mask= zeros(numel(tmp),numel(pattern));

		for el = 1:numel(tmp)
			mask(el,:) = ~cellfun(@isempty,regexp(tmp(el),pattern));
		end

		mask = logical(prod(mask,2));

		fnames = fnames(mask);
end

