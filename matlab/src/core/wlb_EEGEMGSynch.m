function status = wlb_EEGEMGSynch(varargin)
%BN_TENSSYNCH synchronize multiple files in the same vhdr/eeg file
%	FLAG = BN_TENSSYNCH(VARARGIN) I'll edit it when it will be ready
%

% Edited 2014-09-22 by Gabriele Arnulfo <gabriele.arnulfo@gmail.com>
%
%%% DEFINE SUPPORTED FILE EXTENSIONS %%%%
global globDebug;

try_manualTENS = false;

p = inputParser;
p.addRequired('pathEmg',@ischar);
p.addRequired('pathHdeeg',@ischar);
p.addRequired('outdir',@ischar);

p.addOptional('pathEvents','',@ischar)
p.addOptional('fnameFilters',{''},@iscell);

p.addOptional('eegCuttingTime',[0 0],@(x) isnumeric(x) & numel(x)==2);
p.addOptional('emgCuttingTime',[0 0],@(x) isnumeric(x) & numel(x)==2);
p.addOptional('automaticDetection',true,@islogical);
p.addOptional('find_tens_pctg',80,@isnumeric);

p.parse(varargin{:});

pathEmg = p.Results.pathEmg;
pathHdeeg = p.Results.pathHdeeg;
pathEvents = p.Results.pathEvents;
fnameFilters = p.Results.fnameFilters;
pctg=p.Results.find_tens_pctg;
% get filenames within each path
emgFileNames = dir(fullfile(pathEmg,'*emg.*'));
eegFileNames = dir(fullfile(pathHdeeg,'*eeg.eeg'));

eveFileNames = dir(fullfile(pathEvents,'*event*.csv'));

assert(wlb_checkDataConsistency(eegFileNames,emgFileNames),...
    'Number of eeg files and emg files do not match');

if ~isempty(fnameFilters{1})
    emgFileNames = wlb_filterFnames(emgFileNames,fnameFilters);
    eegFileNames = wlb_filterFnames(eegFileNames,fnameFilters);
end

for fileIdx =1:numel(eegFileNames)
    close all
    synched_file = false;
    discarded_file = false;
    
    while(not(synched_file)&&not(discarded_file))
        % files has been renamed to match
        % subject_studyName_drug_trial##_modality.**
        drugCondition = regexp(eegFileNames(fileIdx).name,'_','split');
        
        % pick drug string
        drugCondition = drugCondition{3};
        
        emgFname = fullfile(pathEmg,emgFileNames(fileIdx).name);
        eegFname = fullfile(pathHdeeg,eegFileNames(fileIdx).name);
        
        trialIdx = cell2mat(regexp(eegFname,'trial\d+','match'));
        if isempty( trialIdx )
            % this means that we are considering sets of trials in a
            % single eeg data sesssion
            setIdx = cell2mat(regexp(eegFname,'set\d+','match'));
            setIdx = str2double(cell2mat(regexp(setIdx,'\d+','match')));
            eveFname = eveFileNames(~cellfun(@isempty,regexp({eveFileNames.name},drugCondition)));
            eveFname = eveFname(~cellfun(@isempty,regexp({eveFname.name},strcat('set',num2str(setIdx)))));
            eveFname = fullfile(pathEvents,eveFname.name);
            
        else
            % we here separate trials in single files
            trialIdx = str2double(cell2mat(regexp(trialIdx,'\d+','match')));
            eveFname = eveFileNames(~cellfun(@isempty,regexp({eveFileNames.name},drugCondition)));
            eveFname = fullfile(pathEvents,eveFname.name);
        end
        
        fprintf('Synch:\t%s\n\t%s\n\t%s\n',eegFname,emgFname,eveFname);
        
        eventsInfo = [];
        if ~isempty(p.Results.pathEvents)
            eventsInfo = wlb_readExternalEventFile( eveFname, trialIdx );
        end
        
        %				try
        % read pcs header file
        [emgHdr, emgData, emgEvent] = wlb_readEMG_wue( emgFname );
        [eegHdr, eegData,eegEvent] = wlb_readBrainvision( eegFname );
        
        [emgData,emgHdr,emgEvent] = wlb_cutInitialSamplesData(emgData,emgHdr,emgEvent,p.Results.emgCuttingTime,emgHdr.freq);
        [eegData,eegHdr,eegEvent] = wlb_cutInitialSamplesData(eegData,eegHdr,eegEvent,p.Results.eegCuttingTime,eegHdr.Fs);
        
        % define TENS channels
        emgChIdx = find(strcmp(emgHdr.chantype,'pulse'));
        
        if isempty(emgChIdx)
            error('EMG pulse channel not found');
        end
        
        % extract TENS channels
        eegCh = mean(eegData(:,:));
        emgCh = emgData(emgChIdx,:);
                
        % sometimes EMG or PCS data have been cutted at
        % the begginig or at the end, we then pad with
        % random samples around the mean
        padValsEeg = zeros(1,4*eegHdr.Fs);
        padValsEmg = zeros(1,4*emgHdr.freq);
        
        eegCh = [padValsEeg, eegCh, padValsEeg];
        emgCh = [padValsEmg, emgCh, padValsEmg];
        
        % find time windows containing TENS
        eegLocs = wlb_findTENSArtefact(eegCh,eegHdr.Fs,pctg,'EEG');
        emgLocs = wlb_findTENSArtefact(emgCh,eegHdr.Fs,pctg,'EMG');
        
        method = min([length(eegLocs)/2,length(emgLocs)/2]);
        
        % actually compute t0 for all channels
        dataCell = [{eegCh},{emgCh}];
        
        if p.Results.automaticDetection && not(try_manualTENS)
            t0 = cellfun(@wlb_findTInit,dataCell,{eegLocs,emgLocs},...
                {eegHdr.Fs,emgHdr.freq},...
                {method,method},'uni',false);
        else
            t0 = cellfun(@wlb_manualTENS,dataCell,'uni',false);
            method = min(length(t0{1}),length(t0{2}));
            try_manualTENS = false;
        end
                
        if( method == 2 )
            t0 = reshape([t0{:}],2,2)';
        else
            t0 = [t0{:}]';
            t0(:,2) = [length(eegCh) length(emgCh)];
        end
        
        % remove the padded samples
        offsets = [eegHdr.Fs*4; emgHdr.freq*4];
        t0 = t0 - repmat(offsets,[1 2]);
                                       
        emgChannels = size(emgData,1);
                
        IMUevent_indx = find(cellfun(@(x)not(isempty(x)),strfind({eegEvent.label},'R128')));
        IMUevent = eegEvent(IMUevent_indx);
        
        if((not(isempty(emgEvent)))&&(not(isempty(IMUevent))))
            emgIMUtimes = emgEvent.times;
            eegIMUtimes = [IMUevent.samples]/eegHdr.Fs;
            
            eegIMUlocs = IMUevent.samples;
            emgIMUlocs = emgEvent.samples;

            emg_prepadding = max([0,eegIMUlocs(1)-emgIMUlocs(1)]);
            eeg_prepadding = max([0,emgIMUlocs(1)-eegIMUlocs(1)]);

            emgDataOut = [zeros(size(emgData,1),emg_prepadding) emgData];
            eegDataOut = [zeros(size(eegData,1),eeg_prepadding) eegData];
            
            t0(1,:) = t0(1,:) + eeg_prepadding;
            t0(2,:) = t0(2,:) + emg_prepadding;          

            emg_postpadding = max([0,size(eegDataOut,2)-size(emgDataOut,2)]);
            eeg_postpadding = max([0,size(emgDataOut,2)-size(eegDataOut,2)]);
            
            emgDataOut = [emgDataOut zeros(size(emgDataOut,1),emg_postpadding)];
            eegDataOut = [eegDataOut zeros(size(eegDataOut,1),eeg_postpadding)];
            
            dataOut 	= [eegDataOut; emgDataOut];
            
            dataOut = dataOut(:,1+emg_prepadding:end - emg_postpadding);
            t0 = t0 - emg_prepadding;

        else
            
            emg_prepadding = max([0,t0(1,1)-t0(2,1)]);
            eeg_prepadding = max([0,t0(2,1)-t0(1,1)]);

            emgDataOut = [zeros(size(emgData,1),emg_prepadding) emgData];
            eegDataOut = [zeros(size(eegData,1),eeg_prepadding) eegData];
                         
            t0(1,:) = t0(1,:) + eeg_prepadding;
            t0(2,:) = t0(2,:) + emg_prepadding;

            emg_postpadding = max([0,size(eegDataOut,2)-size(emgDataOut,2)]);
            eeg_postpadding = max([0,size(emgDataOut,2)-size(eegDataOut,2)]);
            
            emgDataOut = [emgDataOut zeros(size(emgDataOut,1),emg_postpadding)];
            eegDataOut = [eegDataOut zeros(size(eegDataOut,1),eeg_postpadding)];
            
            dataOut 	= [eegDataOut; emgDataOut];
            dataOut = dataOut(:,1+emg_prepadding:end - emg_postpadding);
            t0 = t0 - emg_prepadding;
            
        end
          
        emgHdr.chanunit(1:emgChannels) = {'mV'};
        
        outHdr.chanunit = [eegHdr.chanunit, emgHdr.chanunit ];
        
        if globDebug
            figure, clf
            subplot(411)
            plot([-100:100]/fs*1000,mean(dataOut(1:eegHdr.NumberOfChannels,(-100:100) + t0(1,1)),1)/max(mean(dataOut(1:eegHdr.NumberOfChannels,(-100:100) + t0(1,1)),1))),
            hold on, plot([-100:100]/fs*1000,dataOut(eegHdr.NumberOfChannels+emgChIdx,(-100:100) + t0(1,1))/max(dataOut(eegHdr.NumberOfChannels+emgChIdx,(-100:100) + t0(1,1))),'k');
            if(method == 2)
            subplot(412)
            plot([-100:100]/fs*1000,mean(dataOut(1:eegHdr.NumberOfChannels,(-100:100) + t0(1,2)),1)/max(mean(dataOut(1:eegHdr.NumberOfChannels,(-100:100) + t0(1,2)),1))),
            hold on,plot([-100:100]/fs*1000,dataOut(eegHdr.NumberOfChannels+emgChIdx,(-100:100) + t0(1,2))/max(dataOut(eegHdr.NumberOfChannels+emgChIdx,(-100:100) + t0(1,2))),'k');
            end
            
            h1 = subplot(413);
            plot([0:size(dataOut,2)-1]/fs,mean(dataOut(1:eegHdr.NumberOfChannels,:),1)/max(mean(dataOut(1:eegHdr.NumberOfChannels,:),1)));
            hold on, plot((t0(1,1:method)-1)/fs,mean(dataOut(1:eegHdr.NumberOfChannels,t0(1,1:method)),1)/max(mean(dataOut(1:eegHdr.NumberOfChannels,:),1)),'*g');
            h2 = subplot(414);
            plot([0:size(dataOut,2)-1]/fs,dataOut(eegHdr.NumberOfChannels+emgChIdx,:)/max(dataOut(eegHdr.NumberOfChannels+emgChIdx,:)));
            hold on, plot((t0(2,1:method)-1)/fs,dataOut(eegHdr.NumberOfChannels+emgChIdx,t0(2,1:method))/max(dataOut(eegHdr.NumberOfChannels+emgChIdx,:)),'*g');
            drawnow
            title('Synched data');
            linkaxes([h1 h2],'x');
        end
        
        % update header info
        outHdr.label 		= [eegHdr.label;  emgHdr.labels];
        outHdr.nChans		= eegHdr.NumberOfChannels + emgChannels;
        outHdr.NumberOfChannels = outHdr.nChans;
        outHdr.chantype	= [eegHdr.chantype';emgHdr.chantype'];
        outHdr.Fs				= eegHdr.Fs;
                
        emgPosStruct = struct('type','emg',...
            'labels',emgHdr.labels,...
            'sph_theta_besa',-134,...
            'sph_phi_besa',-45);
        
        outHdr.layout.pos = [eegHdr.layout.pos,  emgPosStruct'];
        
        % this is the only supported data format
        outHdr.DataFormat      = 'BINARY';
        outHdr.DataOrientation = 'MULTIPLEXED';
        outHdr.BinaryFormat    = 'IEEE_FLOAT_32';
        
        % no additional calibration needed, since float32
        outHdr.resolution      = ones(size(outHdr.label));
        
        str = input('Do you want to write data (Y), try manual TENS (M) or discard file (N)? Y/M/N','s');
        if(lower(str)=='y')% write data
            [~, fnameEeg, ~] = fileparts(eegFname);
            filename = strcat(fnameEeg,'EegEmg');
            
            outHdr.DataFile = strcat(filename,'.eeg');
            outHdr.MarkerFile = '';
            
            write_brainvision_eeg(p.Results.outdir, outHdr, dataOut);
            
            if(not(isempty(eventsInfo)))
                eventsInfo = rmfield(eventsInfo,'times');
            end
            eventsInfo = [eegEvent eventsInfo];
            
            if(~isempty(eventsInfo))
                outHdr.MarkerFile = strcat(filename,'.vmrk');
                write_brainvision_vmrk(p.Results.outdir, outHdr, eventsInfo);
            end
            
            outHdr.good_chan = ones(outHdr.NumberOfChannels,1);
            write_brainvision_vhdr(p.Results.outdir, outHdr);
            synched_file = true;
        elseif(lower(str)=='m')% write data
            try_manualTENS = true;
        else
            discarded_file = true;
            disp('File discarded');
        end
    end
end % for files

status = 0;


end % function





