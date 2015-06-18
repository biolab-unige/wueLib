
[subjectNames, experiments, drug, trial, fileNames] = textread('PCSgaitIniList.csv','%s %s %s %d %s\n',...
																										'delimiter',',','headerlines',1);

fid = fopen('./conversion.log','w');
for file = 1:numel(fileNames)

		newFileName_ = strjoin([subjectNames(file),experiments(file),...
													drug(file),{strcat('trial',num2str(trial(file)))}],'_');
		newFileName_ = regexprep(newFileName_,'\s','');
		oldFileName_ = fullfile(lower(subjectNames{file}),fileNames{file});
		newFileName_ = fullfile(lower(subjectNames{file}),newFileName_);

		extension = {'.txt','.xml'};

		for ext = 1:2

				oldFileName = strcat(oldFileName_,extension{ext});
				newFileName = strcat(newFileName_,extension{ext});
				[s, mesg] = movefile(oldFileName,newFileName);
				if(~s)
						fprintf(fid,mesg);
				end

		end
end
fclose(fid);
