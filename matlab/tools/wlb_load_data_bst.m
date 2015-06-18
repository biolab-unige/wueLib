function varargout = bn_load_data_bst(varargin)
% bn_load_data_bst

% Edited 2015-04-02 by Gabriele Arnulfo <gabriele.arnulfo@gmail.com>

path = uigetdir();

% read the xls with filenames that should be imported

% we create two conditions (ON/OFF)
% we create the subject
% we add here the trial# to respective condition/subject

files = dir(fullfile(path,	'*.eeg'));

for ii = 1:numel(files)
		filename = files(ii).name;
		subj = filename([2 3]);

		[sSubject, iSubject] = bst_get('Subject',subj);
	
		if(isempty(iSubject))
			[sSubject, iSubject] = db_add_subject(subj,[],0,0);
		end

		import_raw(fullfile(path,filename),'EEG-BRAINAMP',iSubject);
end

end
