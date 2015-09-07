
[csv_match_file,csv_match_file_path] = uigetfile('/opt/data/*.csv','Select the csv matching file');

[subjectNames, experiments, drug, trial, fileNames_pcs,fileNames_emg] = textread(fullfile(csv_match_file_path,...
		csv_match_file),'%s %s %s %d %s %s\n','delimiter',',','headerlines',1);

fid = fopen(fullfile(csv_match_file_path,'conversion.log'),'w');

for file = 1:numel(fileNames_pcs)

		newFileName_pcs = strjoin([subjectNames(file),experiments(file),...
													drug(file),{strcat('trial',num2str(trial(file)))},'pcs'],'_');
		oldFileName_ = fullfile(lower(subjectNames{file}),fileNames_pcs{file});
		newFileName_pcs = fullfile(lower(subjectNames{file}),newFileName_pcs);

		extension = {'.txt','.xml'};
		for ext = 1:numel(extension)

				oldFileName = strcat(oldFileName_,extension{ext});
				newFileName = strcat(newFileName_pcs,extension{ext});
				[s, mesg] = movefile(fullfile(csv_match_file_path,oldFileName),...
                    fullfile(csv_match_file_path,newFileName));
				if(~s)
						fprintf(fid,mesg);
				end
		end
        
    newFileName_emg= strjoin([subjectNames(file),experiments(file),...
													drug(file),{strcat('trial',num2str(trial(file)))},'emg'],'_');
		oldFileName_ = fullfile(lower(subjectNames{file}),fileNames_emg{file});
		newFileName_emg = fullfile(lower(subjectNames{file}),newFileName_emg);

		oldFileName = strcat(oldFileName_,'.txt');
		newFileName = strcat(newFileName_emg,'.txt');
		[s, mesg] = movefile(fullfile(csv_match_file_path,oldFileName),...
								fullfile(csv_match_file_path,newFileName));
		if(~s)
						fprintf(fid,mesg);
		end

end
fclose(fid);
