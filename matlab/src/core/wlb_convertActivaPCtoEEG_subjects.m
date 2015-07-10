function status = wlb_convertActivaPCtoEEG_subejcts(fname)
%WLB_CONVERTACTIVAPCTOEEG_SUBEJCTS Description
%	STATUS = WLB_CONVERTACTIVAPCTOEEG_SUBEJCTS(FNAME) Long description
%
	logFid	= fopen('/media/gabri/My Passport/DBS/logFile.log','w');
	[parent, ~] = fileparts(fname);

	% we are in the data dir containing all the subjects
	[folder subjectNames] = buildFolderStruct(fname);	
	
	for sIdx = 1:numel(subjectNames)

			[visitFolders visitNames] = buildFolderStruct(folder{sIdx});

			for vIdx = 1:numel(visitNames)
			
				[ conditionFolders conditionNames] = buildFolderStruct(visitFolders{vIdx});

				for cIdx = 1:numel(conditionNames)

						wlb_convertActivaPCtoEEG(conditionFolders{cIdx},logFid);
				
				end; clear cIdx;

			end; clear vIdx;
	
	end; clear sIdx;

fclose(logFid);

end

function [folder names] = buildFolderStruct(folder)
%BUILDFOLDERSTRUCT Description
%	FOLDER = BUILDFOLDERSTRUCT(FOLDER) Long description
%
			parent = folder;
	
			folder = dir(folder);
			folder = clearFoldernames(folder);

			names = {folder.name};
			folder = cellfun(@fullfile,repmat({parent},size(names)),names,'Uni',false);

end


function folder = clearFoldernames(folder)
%CLEARFOLDERNAMES Description
%	folder = CLEARFOLDERNAMES(FOLDER) Long description
%

	mask1 = ~cellfun(@isempty,regexp({folder.name},'\w+'));
	mask2 = cellfun(@isempty,regexp({folder.name},'[.].+'));

	folder = folder(mask1 & mask2);

end
