
files = dir(fullfile(pwd,'*.vhdr'));

for i = 1:numel(files)
    
    [dum,oldfilename] = fileparts(files(i).name);
    [c, matches] = strsplit(oldfilename,'_');
    subj = strrep(c{1},'W','wue');
    condition = (c{2});
    drug = [lower(c{4}) '-' lower(c{5})];
    trial = ['trial' num2str(str2num(c{7}))];
    type = 'pcsEeg';
    
    newfilename = strjoin({subj,condition,drug,trial,type},'_');
    
    movefile([oldfilename '.eeg'],[newfilename '.eeg']);
    movefile([oldfilename '.vmrk'],[newfilename '.vmrk']);
    movefile([oldfilename '.vhdr'],[newfilename '.vhdr']);
    
    [eegHdr, eegData, eegEve] = wlb_readBrainvision([newfilename '.vhdr']);
    write_brainvision_eeg(pwd,eegHdr,eegData);
    write_brainvision_vmrk(pwd,eegHdr,eegEve);
    write_brainvision_vhdr(pwd,eegHdr);
end

%     [dum,oldfilename] = fileparts(files(i).name);
%     [c, matches] = strsplit(oldfilename,'_');
%     subj = c{1};
%     condition = [(c{2}) '-' (c{3})];
%     if(strcmp(c{4},'off'))
%         drug = 'medoff-stimoff';
%     else
%         drug = 'medon-stimoff';
%     end
%     trial = c{5};
%     type = c{6};
%     
%     newfilename = strjoin({subj,condition,drug,trial,type},'_');
%     
%     movefile([oldfilename '.eeg'],[newfilename '.eeg']);
%     movefile([oldfilename '.vmrk'],[newfilename '.vmrk']);
%     movefile([oldfilename '.vhdr'],[newfilename '.vhdr']);
%     
%     [eegHdr, eegData, eegEve] = wlb_readBrainvision([newfilename '.vhdr']);
%     write_brainvision_eeg(pwd,eegHdr,eegData);
%     write_brainvision_vmrk(pwd,eegHdr,eegEve);
%     write_brainvision_vhdr(pwd,eegHdr);
% end