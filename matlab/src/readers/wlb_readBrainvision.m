function [eegHdr, eegData, eegEve] = wlb_readBrainvision(filename)
%WLB_READBRAINVISION Description
%	[EEGHDR, EEGDATA, EEGEVE] = WLB_READBRAINVISION(FILENAME) Long description
%
	
	[p,f,~] = fileparts(filename);
	eegHdr = read_brainvision_vhdr(fullfile(p,strcat(f,'.vhdr')));
	eegData= read_brainvision_eeg(fullfile(p,strcat(f,'.eeg')),eegHdr,1,eegHdr.nSamples);
	eegEve = read_brainvision_vmrk(fullfile(p,strcat(f,'.vmrk')));

end
