function popolate_brainstorm_db()

subjects = {'wue02','wue05','wue09','wue10','wue11'};
% subjects = {'wue02'};
cond = {'D1RANLEFT','D1RANRIGHT','D1ROTRIGHT','D2RANLEFT','D2RANRIGHT','D2ROTRIGHT'};
% cond = {'D1RANLEFT'};
ecg_channel = {'F4','F4','F4','F4','FFC3h'};

for iSub = 1:length(subjects)
        
    [sSubject, iSubject] = db_add_subject(subjects{iSub});
    for iCond = 1:length(cond)
        
        iStudies = db_add_condition(iSubject,cond{iCond});

        pathHdeeg = fullfile(pwd,subjects{iSub},'hdeeg_pcs_sync');
        pathEventFile = fullfile(pwd,subjects{iSub},'mtmdata');
        outdir = fullfile(pwd,subjects{iSub},'hdeeg_pcs_sync');
        eegFname = MTM_append_event_main('pathHdeeg',pathHdeeg,'pathEventFile',pathEventFile,'outdir',outdir,...
            'fNameFilters',{cond{iCond}},'ecg_channel',ecg_channel{iSub});     
        ImportOptions = db_template('ImportOptions');
        ImportOptions.CreateConditions = 0;
        ImportOptions.DisplayMessages = 0;
        ImportOptions.ChannelReplace = 0;
        OutputDataFile = import_raw(eegFname, 'EEG-BRAINAMP', iSubject,iStudies, ImportOptions);

    end
end
end

function OutputDataFile = import_raw(RawFiles, FileFormat, iSubject,iOutputStudy, ImportOptions)
% IMPORT_RAW: Create a link to a raw file in the Brainstorm database.
%
% USAGE:  OutputDataFile = import_raw(RawFiles, FileFormat, iSubject, ImportOptions)
%         OutputDataFile = import_raw(RawFiles, FileFormat, iSubject)
%         OutputDataFile = import_raw(RawFiles, FileFormat)
%         OutputDataFile = import_raw()
%
% INPUTS:
%     - RawFiles      : Full path to the file to import in database
%     - FileFormat    : String representing the file format (CTF, FIF, 4D, ...)
%     - iSubject      : Subject indice in which to import the raw file
%     - ImportOptions : Structure that describes how to import the recordings
%       => Fields used: ChannelAlign, ChannelReplace, DisplayMessages, EventsMode, EventsTrackMode

% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2015 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
% 
% Authors: Francois Tadel, 2009-2014

%% ===== PARSE INPUT =====
if (nargin < 4) || isempty(ImportOptions)
    ImportOptions = db_template('ImportOptions');
end
if (nargin < 3)
    iSubject = [];
end
if (nargin < 2)
    RawFiles = [];
    FileFormat = [];
end
% Force list of files to be a cell array
if ~isempty(RawFiles) && ischar(RawFiles)
    RawFiles = {RawFiles};
end
% Some verifications
if ~isempty(RawFiles) && isempty(FileFormat)
    error('If you pass the filenames in input, you must define also the FileFormat argument.');
end
% Get Protocol information
ProtocolInfo = bst_get('ProtocolInfo');
OutputDataFile = [];


%% ===== SELECT DATA FILE =====
% If file to load was not defined : open a dialog box to select it
if isempty(RawFiles) 
    % Get default import directory and formats
    LastUsedDirs = bst_get('LastUsedDirs');
    DefaultFormats = bst_get('DefaultFormats');
    % Get MRI file
    [RawFiles, FileFormat, FileFilter] = java_getfile( 'open', ...
        'Open raw EEG/MEG recordings...', ...  % Window title
        LastUsedDirs.ImportData, ...           % Last used directory
        'multiple', 'files_and_dirs', ...      % Selection mode
        bst_get('FileFilters', 'raw'), ...     % List of available file formats
        DefaultFormats.DataIn);
    % If no file was selected: exit
    if isempty(RawFiles)
        return
    end
    % Save default import directory
    LastUsedDirs.ImportData = bst_fileparts(RawFiles{1});
    bst_set('LastUsedDirs', LastUsedDirs);
    % Save default import format
    DefaultFormats.DataIn = FileFormat;
    bst_set('DefaultFormats',  DefaultFormats);
    % Process the selected directories :
    %    1) If they are .ds/ directory with .meg4 and .res4 files : keep them as "files to open"
    %    2) Else : add all the data files they contains (subdirectories included)
    RawFiles = io_expand_filenames(FileFilter, RawFiles);
    if isempty(RawFiles)
        bst_error(['No data ' FileFormat ' file in the selected directories.'], 'Open raw EEG/MEG recordings...', 0);
        return
    end

    % ===== SUB-CATEGORIES IN FILE FORMAT =====
    if strcmpi(FileFormat, 'EEG-NEUROSCAN')
        [tmp, tmp, fileExt] = bst_fileparts(RawFiles{1});
        % Switch between different Neuroscan formats
        switch (lower(fileExt))
            case '.cnt',  FileFormat = 'EEG-NEUROSCAN-CNT';
            case '.eeg',  FileFormat = 'EEG-NEUROSCAN-EEG';
            case '.avg',  FileFormat = 'EEG-NEUROSCAN-AVG';
            case '.dat',  FileFormat = 'EEG-NEUROSCAN-DAT';
        end
    end
end


%% ===== IMPORT =====
% Loop on the files to import
for iFile = 1:length(RawFiles)
    % ===== OPENING FILE =====
    bst_progress('start', 'Open raw EEG/MEG recordings', 'Reading file header...');
    % Open file
    [sFile, ChannelMat, errMsg] = in_fopen(RawFiles{iFile}, FileFormat, ImportOptions);
    if isempty(sFile)
        bst_progress('stop');
        return;
    end
    % Yokogawa non-registered warning
    if ~isempty(errMsg) && ImportOptions.DisplayMessages
        java_dialog('warning', errMsg, 'Open raw EEG/MEG recordings');
    end

    % ===== OUTPUT STUDY =====
    % Get short filename
    [tmp, fBase] = bst_fileparts(RawFiles{iFile});
    % Build output condition name
    if isfield(sFile, 'condition') && ~isempty(sFile.condition)
        ConditionName = ['@raw' sFile.condition];
    else
        ConditionName = ['@raw' fBase];
    end
    % Output subject
    if isempty(iSubject)
        % Get default subject
        SubjectName = 'NewSubject';
        [sSubject, iSubject] = bst_get('Subject', SubjectName, 1);
        % If subject does not exist yet: create it
        if isempty(sSubject)
            [sSubject, iSubject] = db_add_subject(SubjectName);
        end
        % If subject cannot be created
        if isempty(sSubject)
            bst_error(['Could not create subject "' SubjectName '"'], 'Open raw EEG/MEG recordings');
            return;
        end
    else
        % Get specified subject
        sSubject = bst_get('Subject', iSubject, 1);
    end
    % Do not allow automatic registration with head points when using the default anatomy
    if (sSubject.UseDefaultAnat)
        ImportOptions.ChannelAlign = 0;
    end

    % If condition already exists
    [sExistStudy, iExistStudy] = bst_get('StudyWithCondition', bst_fullfile(sSubject.Name, file_standardize(ConditionName, 1)));
    if ~isempty(sExistStudy) && ~isempty(sExistStudy.Data)
        % Need to check if the raw file is the same or they are two files with the same name in different folders
        % Get the raw data files
        iRaw = find(strcmpi({sExistStudy.Data.DataType}, 'raw'));
        if ~isempty(iRaw)
            % Load data description
            DataFile = sExistStudy.Data(iRaw).FileName;
            DataMat = in_bst_data(DataFile);
            % If same filenames: cannot link it again in the database
            LinkFile = DataMat.F.filename;
            minLength = min(length(LinkFile), length(RawFiles{iFile}));
            if file_compare(LinkFile(1:minLength), RawFiles{iFile}(1:minLength))
                %bst_error('This file is already available in the explorer.', 'Open raw EEG/MEG recordings', 0);
                panel_protocols('SelectNode', [], 'rawdata', iExistStudy, iRaw );
                OutputDataFile = DataFile;
                continue;
            % Else: Create a condition with a different name
            else
                % Add a numeric tag at the end of the condition name
                curPath = bst_fullfile(ProtocolInfo.STUDIES, bst_fileparts(sExistStudy.FileName));
                curPath = file_unique(curPath);
                [tmp__, ConditionName] = bst_fileparts(curPath, 1);
            end
        end
    end
    % Create output condition
%     iOutputStudy = db_add_condition(sSubject.Name, ConditionName);
    if isempty(iOutputStudy)
        error('Folder could not be created : "%s/%s".', bst_fileparts(sSubject.FileName), ConditionName);
    end
    % Get output study
    sOutputStudy = bst_get('Study', iOutputStudy);
    % Get the study in which the channel file has to be saved
    [sChannel, iChannelStudy] = bst_get('ChannelForStudy', iOutputStudy);
    
    % ===== CREATE DEFAULT CHANNEL FILE =====
    if isempty(ChannelMat)
        ChannelMat = db_template('channelmat');
        ChannelMat.Comment = [sFile.device ' channels'];
        ChannelMat.Channel = repmat(db_template('channeldesc'), [1, length(sFile.channelflag)]);
        % For each channel
        for i = 1:length(ChannelMat.Channel)
            if (length(ChannelMat.Channel) > 99)
                ChannelMat.Channel(i).Name = sprintf('E%03d', i);
            else
                ChannelMat.Channel(i).Name = sprintf('E%02d', i);
            end
            ChannelMat.Channel(i).Type    = 'EEG';
            ChannelMat.Channel(i).Loc     = [0; 0; 0];
            ChannelMat.Channel(i).Orient  = [];
            ChannelMat.Channel(i).Weight  = 1;
            ChannelMat.Channel(i).Comment = [];
        end
        % Add channel file to database
        db_set_channel(iChannelStudy, ChannelMat, 0, 0);

    % ===== SAVE LOADED CHANNEL FILE =====
    else
        % Add history field to channel structure
        ChannelMat = bst_history('add', ChannelMat, 'import', ['Link to file: ' RawFiles{iFile} ' (Format: ' FileFormat ')']);
        % Remove fiducials only from polhemus and ascii files
        isRemoveFid = ismember(FileFormat, {'MEGDRAW', 'POLHEMUS', 'ASCII_XYZ', 'ASCII_NXYZ', 'ASCII_XYZN', 'ASCII_NXY', 'ASCII_XY', 'ASCII_NTP', 'ASCII_TP'});
        % Detect auxiliary EEG channels
        ChannelMat = channel_detect_type(ChannelMat, 0, isRemoveFid);
        % Do not align data coming from Brainstorm exported files (already aligned)
        if strcmpi(FileFormat, 'BST-BIN')
            ImportOptions.ChannelAlign = 0;
        end
        % Add channel file to database
        [ChannelFile, ChannelMat, ImportOptions.ChannelReplace, ImportOptions.ChannelAlign, Modality] = db_set_channel(iChannelStudy, ChannelMat, ImportOptions.ChannelReplace, ImportOptions.ChannelAlign);
        % Display the registration if this was skipped in db_set_channel
        if (ImportOptions.ChannelAlign == 0) && (ImportOptions.DisplayMessages) && ~isempty(Modality)
            bst_memory('UnloadAll', 'Forced');
            channel_align_manual(ChannelFile, Modality, 0);
        end
    end
    
    % ===== SAVE DATA FILE =====
    % Build output filename
    OutputDataFile = bst_fullfile(ProtocolInfo.STUDIES, bst_fileparts(sOutputStudy.FileName), ['data_0raw_' fBase '.mat']);
    % Build output structure
    DataMat = db_template('DataMat');
    DataMat.F           = sFile;
    DataMat.Comment     = 'Link to raw file';
    DataMat.ChannelFlag = sFile.channelflag;
    DataMat.Time        = sFile.prop.times;
    DataMat.DataType    = 'raw';
    DataMat.Device      = sFile.device;
    % Compumedics: add start time to the file comment
    if strcmpi(sFile.format, 'EEG-COMPUMEDICS-PFS') && isfield(sFile.header, 'rda_startstr') && ~isempty(sFile.header.rda_startstr)
        DataMat.Comment = [DataMat.Comment ' [' sFile.header.rda_startstr ']'];
    end
    % Add history field
    DataMat = bst_history('add', DataMat, 'import', ['Link to raw file: ' RawFiles{iFile}]);
    % Save file on hard drive
    bst_save(OutputDataFile, DataMat, 'v6');
    % Add file to database
    sOutputStudy = db_add_data(iOutputStudy, OutputDataFile, DataMat);

    % ===== UPDATE DATABASE =====
    % Update links
    db_links('Study', iOutputStudy);
    % Refresh both data node and channel node
    iUpdateStudies = unique([iOutputStudy, iChannelStudy]);
    panel_protocols('UpdateNode', 'Study', iUpdateStudies);
end

% If something was imported
if ~isempty(iOutputStudy)
    % Select the data study node
    panel_protocols('SelectStudyNode', iOutputStudy);
    % Save database
    db_save();
end

bst_progress('stop');
end

