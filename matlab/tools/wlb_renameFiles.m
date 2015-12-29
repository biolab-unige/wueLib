% USAGE:
%
%
[csvMatchFile,csvMatchFilePath] = uigetfile('/opt/data/*.csv','Select the csv matching file');

[subjectNames, experiments, drug, trial, fileNamesPcs,fileNamesEmg,fileNamesEeg] = textread(fullfile(csvMatchFilePath,...
    csvMatchFile),'%s %s %s %d %s %s %s','delimiter',',','headerlines',1);

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
    
    newFileNameEeg= strjoin([subjectNames(file),experiments(file),...
            drug(file),{strcat('trial',num2str(trial(file)))},'eeg'],'_');
    oldFileName_ = fullfile(lower(subjectNames{file}),fileNamesEeg{file});
    newFileNameEeg = fullfile(lower(subjectNames{file}),newFileNameEeg);

    extension = {'.eeg','.vhdr','.vmrk'};
    
    for ext = 1:numel(extension)
    
    oldFileName = strcat(oldFileName_,extension{ext});
    newFileName = strcat(newFileNameEeg,extension{ext});
    [s, mesg] = movefile(fullfile(csvMatchFilePath,oldFileName),...
        fullfile(csvMatchFilePath,newFileName));
    if(~s)
        fprintf(fid,mesg);
    end
    if(strcmp(extension{ext},'.vhdr'))
        [hdr] = read_brainvision_vhdr(fullfile(csvMatchFilePath,newFileName));
        newFileName = strcat(newFileNameEeg,'.eeg');
        hdr.DataFile = newFileName;
        if(not(strcmp(hdr.MarkerFile,'')))
            newFileName = strcat(newFileNameEeg,'.vmrk');
            hdr.MarkerFile = newFileName;
        end
        write_brainvision_vhdr(fullfile(csvMatchFilePath,lower(subjectNames{file})), hdr);
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
    else
        disp('found\n')
    end
    
end
fclose(fid);
