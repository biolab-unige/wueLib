function eegFname = MTM_append_eeg(varargin)
%BN_TENSSYNCH synchronize multiple files in the same vhdr/eeg file
%	FLAG = BN_TENSSYNCH(VARARGIN) I'll edit it when it will be ready
%

% Edited 2014-09-22 by Gabriele Arnulfo <gabriele.arnulfo@gmail.com>
%
%%% DEFINE SUPPORTED FILE EXTENSIONS %%%%

p = inputParser;
p.addRequired('pathHdeeg',@ischar);
p.addRequired('outdir',@ischar);

p.addOptional('fnameFilters',{''},@iscell);
p.addOptional('ecg_channel','',@ischar);

p.parse(varargin{:});

pathHdeeg = p.Results.pathHdeeg;
fnameFilters = p.Results.fnameFilters;
ecg_channel  = p.Results.ecg_channel;
% get filenames within each path
eegFileNames = dir(fullfile(pathHdeeg,'*.eeg'));

if ~isempty(fnameFilters{1})
    eegFileNames = filterFnames(eegFileNames,fnameFilters);
end

    eegFname = fullfile(pathHdeeg,eegFileNames(1).name);
    [eegHdr, eegData,vmrk_event] = wlb_readBrainvision( eegFname );    
    ecg_channel = find(strcmp(eegHdr.label,ecg_channel));
    eegHdr.chantype(end-1:end)={'lfp','lfp'};
    eegHdr.chantype(ecg_channel)={'ecg'};
    write_brainvision_vhdr(pathHdeeg,eegHdr);
    
for fileIdx =2:numel(eegFileNames)
    
    eegFname2 = fullfile(pathHdeeg,eegFileNames(fileIdx).name);
        
    [eegHdr, eegData,vmrk_event] = concat_eeg_file(eegFname,eegFname2,ecg_channel);
   
    write_brainvision_eeg(p.Results.outdir, eegHdr, eegData);

    write_brainvision_vmrk(p.Results.outdir, eegHdr, vmrk_event);

    write_brainvision_vhdr(p.Results.outdir, eegHdr);

    eegFname = fullfile(p.Results.outdir,eegFileNames(1).name);

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

