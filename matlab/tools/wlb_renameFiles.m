% USAGE:
%
% 
[csvMatchFile,csvMatchFilePath] = uigetfile('/opt/data/*.csv','Select the csv matching file');

[subjectNames, experiments, drug, trial, fileNamesPcs,fileNamesEmg] = textread(fullfile(csvMatchFilePath,...
		csvMatchFile),'%s %s %s %d %s %s \n','delimiter',',','headerlines',1);

fid = fopen(fullfile(csvMatchFilePath,'conversion.log'),'w');
% for each pcs data
for file = 1:numel(fileNamesPcs)

		newFileNamePcs = strjoin([subjectNames(file),experiments(file),...
													drug(file),{strcat('trial',num2str(trial(file)))},'pcs'],'_');
		oldFileName_ = fullfile(lower(subjectNames{file}),fileNamesPcs{file});
		newFileNamePcs = fullfile(lower(subjectNames{file}),newFileNamePcs);

		extension = {'.txt','.xml'};
		for ext = 1:numel(extension)

				oldFileName = strcat(oldFileName_,extension{ext});
				newFileName = strcat(newFileNamePcs,extension{ext});
				[s, mesg] = movefile(fullfile(csvMatchFilePath,oldFileName),...
                    fullfile(csvMatchFilePath,newFileName));
				if(~s)
						fprintf(fid,mesg);
				end
		end
        
    newFileNameEmg= strjoin([subjectNames(file),experiments(file),...
													drug(file),{strcat('trial',num2str(trial(file)))},'emg'],'_');
		oldFileName_ = fullfile(lower(subjectNames{file}),fileNamesEmg{file});
		newFileNameEmg = fullfile(lower(subjectNames{file}),newFileNameEmg);

		oldFileName = strcat(oldFileName_,'.txt');
		newFileName = strcat(newFileNameEmg,'.txt');
		[s, mesg] = movefile(fullfile(csvMatchFilePath,oldFileName),...
								fullfile(csvMatchFilePath,newFileName));
		if(~s)
						fprintf(fid,mesg);
		end

end
fclose(fid);
