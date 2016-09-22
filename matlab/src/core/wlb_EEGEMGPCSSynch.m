function status = wlb_EEGEMGPCSSynch(varargin)
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
p.addRequired('pathPcs',@ischar);
p.addRequired('pathHdeeg',@ischar);
p.addRequired('outdir',@ischar);

p.addOptional('pathEvents','',@ischar)
p.addOptional('fnameFilters',{''},@iscell);

p.addOptional('pcsRefChannel',1,@isnumeric);
p.addOptional('emgRefChannel','pulse',@ischar);
p.addOptional('pcsCuttingTime',[0 0],@(x) isnumeric(x) & numel(x)==2 );
p.addOptional('eegCuttingTime',[0 0],@(x) isnumeric(x) & numel(x)==2);
p.addOptional('emgCuttingTime',[0 0],@(x) isnumeric(x) & numel(x)==2);
p.addOptional('automaticDetection',true,@islogical);
p.addOptional('find_tens_pctg',80,@isnumeric);

p.parse(varargin{:});

pathPcs = p.Results.pathPcs;
pathEmg = p.Results.pathEmg;
pathHdeeg = p.Results.pathHdeeg;
pathEvents = p.Results.pathEvents;
fnameFilters = p.Results.fnameFilters;
pctg=p.Results.find_tens_pctg;
% get filenames within each path
pcsFileNames = dir(fullfile(pathPcs,'*pcs.xml'));
emgFileNames = dir(fullfile(pathEmg,'*emg.*'));
eegFileNames = dir(fullfile(pathHdeeg,'*eeg.eeg'));

eveFileNames = dir(fullfile(pathEvents,'*event*.csv'));

assert(wlb_checkDataConsistency(pcsFileNames,emgFileNames),...
    'Number of emg files and pcs files do not match');

assert(wlb_checkDataConsistency(pcsFileNames,eegFileNames),...
    'Number of hdeeg files and pcs files do not match');


if ~isempty(fnameFilters{1})
    pcsFileNames = wlb_filterFnames(pcsFileNames,fnameFilters);
    emgFileNames = wlb_filterFnames(emgFileNames,fnameFilters);
    eegFileNames = wlb_filterFnames(eegFileNames,fnameFilters);
end

for fileIdx =1:numel(pcsFileNames)
    close all
    synched_file = false;
    discarded_file = false;
    
    while(not(synched_file)&&not(discarded_file))
        % files has been renamed to match
        % subject_studyName_drug_trial##_modality.**
        drugCondition = regexp(pcsFileNames(fileIdx).name,'_','split');
        
        % pick drug string
        drugCondition = drugCondition{3};
        
        pcsFname = fullfile(pathPcs,pcsFileNames(fileIdx).name);
        emgFname = fullfile(pathEmg,emgFileNames(fileIdx).name);
        eegFname = fullfile(pathHdeeg,eegFileNames(fileIdx).name);
        
        trialIdx = cell2mat(regexp(pcsFname,'trial\d+','match'));
        if isempty( trialIdx )
            % this means that we are considering sets of trials in a
            % single eeg data sesssion
            setIdx = cell2mat(regexp(pcsFname,'set\d+','match'));
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
        
        fprintf('Synch:\t%s\n\t%s\n\t%s\n\t%s\n',pcsFname,emgFname,eegFname,eveFname);
                
        %				try
        % read pcs header file
        [pcsHdr, pcsData] = wlb_readActivaPC( pcsFname );
        [emgHdr, emgData, emgEvent] = wlb_readEMG_wue( emgFname );
        [eegHdr, eegData,eegEvent] = wlb_readBrainvision( eegFname );
        
        eventsInfo = [];
        if ~isempty(p.Results.pathEvents)
            eventsInfo = wlb_readExternalEventFile(emgHdr.freq eveFname, trialIdx );
        end

        
        pcsData = wlb_cutInitialSamplesData(pcsData,p.Results.pcsCuttingTime,pcsHdr.SenseChannelConfig.TDSampleRate);
        emgData = wlb_cutInitialSamplesData(emgData,p.Results.emgCuttingTime,emgHdr.freq);
        eegData = wlb_cutInitialSamplesData(eegData,p.Results.eegCuttingTime,eegHdr.Fs);
        
        % define TENS channels
        pcsChIdx = p.Results.pcsRefChannel;
        emgChIdx = find(strcmp(emgHdr.chantype,p.Results.emgRefChannel));
        
        if isempty(pcsChIdx )
            error('Invalid PC+S reference channel');
        end
        if isempty(emgChIdx)
            error('EMG pulse channel not found');
        end
        
        % extract TENS channels
        pcsCh = pcsData(pcsChIdx,:);
        eegCh = mean(eegData(:,:));
        emgCh = emgData(emgChIdx,:);
        
        % pick central window
        pcsWnd = round(numel(pcsCh)/2)+[-1/2*pcsHdr.SenseChannelConfig.TDSampleRate:...
            1/2*pcsHdr.SenseChannelConfig.TDSampleRate];
        eegWnd = round(numel(eegCh)/2)+[-1/2*eegHdr.Fs:1/2*eegHdr.Fs];
        emgWnd = round(numel(emgCh)/2)+[-1/2*emgHdr.freq:1/2*emgHdr.freq];
        
        pcsStd = std(pcsCh(pcsWnd));
        eegStd = std(eegCh(eegWnd));
        emgStd = std(emgCh(emgWnd));
        
%         % sometimes EMG or PCS data have been cutted at
%         % the begginig or at the end, we then pad with
%         % random samples around the mean
%         padValsPcs = pcsStd.*zeros(1,4*pcsHdr.SenseChannelConfig.TDSampleRate);
%         padValsEeg = eegStd.*zeros(1,4*eegHdr.Fs);
%         padValsEmg = emgStd.*zeros(1,4*emgHdr.freq);
%         
%         pcsCh = [padValsPcs, pcsCh, padValsPcs];
%         eegCh = [padValsEeg, eegCh, padValsEeg];
%         emgCh = [padValsEmg, emgCh, padValsEmg];
        % find time windows containing TENS
        eegLocs = wlb_findTENSArtefact(eegCh,eegHdr.Fs,pctg,'EEG');
        pcsLocs = wlb_findTENSArtefact(pcsCh,pcsHdr.SenseChannelConfig.TDSampleRate,pctg,'PC+S');
        emgLocs = wlb_findTENSArtefact(emgCh,eegHdr.Fs,pctg,'EMG');
        
        method = ([length(eegLocs)/2,length(pcsLocs)/2,length(emgLocs)/2]);
        
        % actually compute t0 for all channels
        dataCell = [{eegCh},{pcsCh},{emgCh}];
        
        if p.Results.automaticDetection && not(try_manualTENS)
            tinit = cellfun(@wlb_findTInit,dataCell,{eegLocs,pcsLocs,emgLocs},...
                {eegHdr.Fs,pcsHdr.SenseChannelConfig.TDSampleRate,emgHdr.freq},...
                {method(1), method(2),method(3)},'uni',false);
        else
            tinit = cellfun(@wlb_manualTENS,dataCell,'uni',false);
            method = [length(tinit{1}),length(tinit{2}),length(tinit{3})];
            try_manualTENS = false;
        end
               
        t0 = NaN(3,2);
        
        if( method(1) == 2 && method(2) == 2)
            t0(1,:) = tinit{1};
            t0(2,:) = tinit{2};
            % estimate the correct fs for PCS
            pcsFs = (t0(2,2)-t0(2,1))* eegHdr.Fs /(t0(1,2)-t0(1,1));
        elseif( method(1) == 1 && method(2) == 1)
            t0(1,:) = tinit{1};
            t0(2,:) = tinit{2};
            t0(1:2,2) = NaN;
            if pcsHdr.SenseChannelConfig.TDSampleRate > 422
                pcsFs = 793.65;
            else
                pcsFs = 422;
            end
        else
            try_manualTENS = true;
            disp('Trying Manual Tens identification');
            continue;
        end
        
        t0(3,1:method(3)) = tinit{3};
        
%         % remove the padded samples
%         offsets = [eegHdr.Fs*4; pcsHdr.SenseChannelConfig.TDSampleRate*4;emgHdr.freq*4];
%         t0 = t0 - repmat(offsets,[1 2]);
               
        % downsample pcs data to integer sampling frequency and eeg accordingly
        fs 				= 400;
        pcsData         = wlb_ResampleCascade(pcsData,fs,pcsFs);
        eegData 		= wlb_ResampleCascade(eegData,fs,eegHdr.Fs);
        emgData 		= wlb_ResampleCascade(emgData,fs,emgHdr.freq);
        
        % each t0 row represent eeg,pcs data before
        % we have to recompute the exact point in time after
        % resampling
        t0(1,:) = (round(t0(1,:)/eegHdr.Fs*fs));
        t0(2,:) = (round(t0(2,:)/pcsFs*fs));
        t0(3,:) = (round(t0(3,:)/emgHdr.freq*fs));
        
        % t0 has start and end samples for each channel
        % cut data to have the same number of samples
        % We compute the smallest offset possible in the initial
        % samples prior each TENS artefact
        eegDataOut = [zeros(size(eegData,1),max([0,t0(2,1)-t0(1,1)])) eegData];
        pcsDataOut = [zeros(size(pcsData,1),max([0,t0(1,1)-t0(2,1)])) pcsData];
        
        eegDataOut = [eegDataOut zeros(size(eegDataOut,1),max([0,size(pcsDataOut,2)-size(eegDataOut,2)]))];
        pcsDataOut = [pcsDataOut zeros(size(pcsDataOut,1),max([0,size(eegDataOut,2)-size(pcsDataOut,2)]))];
  
        dataOut 	= [eegDataOut;pcsDataOut];
        
        temp = t0(1,:) + max([0,t0(2,1)-t0(1,1)]);
        t0(2,:) = t0(2,:) + max([0,t0(1,1)-t0(2,1)]);
        t0(1,:) = temp;
        
        pcsChannels	= size(pcsData,1);
        emgChannels = size(emgData,1);
        
        % if we recorded only one side with PCs
        % we just add an empty row to prevent channel size
        % incosisntencies across modalities and to prevent errors
        % in cases we mixed recording techinques such as one trial with 1 STN
        % and another with 2 STN
        if(pcsChannels == 1)
            dataOut  = [dataOut; zeros(1,size(dataOut,2))];
            pcsChannels = 2;
        end
        
        IMUevent_indx = find(cellfun(@(x)not(isempty(x)),strfind({eegEvent.label},'R128')));
        IMUevent = eegEvent(IMUevent_indx);
        
        if((not(isempty(emgEvent)))&&(not(isempty(IMUevent))))
            emgIMUtimes = emgEvent.times;
            eegIMUtimes = [IMUevent.samples]/eegHdr.Fs;
            
            eegIMUlocs = round(eegIMUtimes*fs) + max([0,t0(2,1)-t0(1,1)]);
            emgIMUlocs = round(emgIMUtimes*fs);

            emg_prepadding = max([0,eegIMUlocs(1)-emgIMUlocs(1)]);
            dataOut_prepadding = max([0,emgIMUlocs(1)-eegIMUlocs(1)]);

            emgDataOut = [zeros(size(emgData,1),emg_prepadding) emgData];
            dataOut = [zeros(size(dataOut,1),dataOut_prepadding) dataOut];
%             emgDataOut = emgData(:,[emgIMUlocs(1):emgIMUlocs(1)+min(diff(eegIMUlocs),diff(emgIMUlocs))]);
%             dataOut = dataOut(:,[eegIMUlocs(1):eegIMUlocs(1)+min(diff(eegIMUlocs),diff(emgIMUlocs))]);

            t0([1 2],:) = t0([1 2],:) + dataOut_prepadding;
            t0(3,:) = t0(3,:) + emg_prepadding;          
%             t0([1 2],:) = t0([1 2],:) - eegIMUlocs(1)+1;
%             t0(3,:) = t0(3,:) - emgIMUlocs(1)+1;

            emg_postpadding = max([0,size(dataOut,2)-size(emgDataOut,2)]);
            dataOut_postpadding = max([0,size(emgDataOut,2)-size(dataOut,2)]);
            
            emgDataOut = [emgDataOut zeros(size(emgDataOut,1),emg_postpadding)];
            dataOut = [dataOut zeros(size(dataOut,1),dataOut_postpadding)];
            
            dataOut 	= [dataOut; emgDataOut];
            
            dataOut = dataOut(:,1+emg_prepadding:end - 0*emg_postpadding);
            t0 = t0 - emg_prepadding;

        else
            if(method(3) == method(1))
                emg_prepadding = max([0,t0(1,1)-t0(3,1)]);
                dataOut_prepadding = max([0,t0(3,1)-t0(1,1)]);
                
                emgDataOut = [zeros(size(emgData,1),emg_prepadding) emgData];
                dataOut = [zeros(size(dataOut,1),dataOut_prepadding) dataOut];
                
                temp = t0(3,:) + emg_prepadding;
                t0([1 2],:) = t0([1 2],:) + dataOut_prepadding;
                t0(3,:) = temp;
                
                emg_postpadding = max([0,size(dataOut,2)-size(emgDataOut,2)]);
                dataOut_postpadding = max([0,size(emgDataOut,2)-size(dataOut,2)]);
                
                emgDataOut = [emgDataOut zeros(size(emgDataOut,1),emg_postpadding)];
                dataOut = [dataOut zeros(size(dataOut,1),dataOut_postpadding)];
                
                dataOut 	= [dataOut; emgDataOut];
                dataOut = dataOut(:,1+emg_prepadding:end - emg_postpadding);
                t0 = t0 - emg_prepadding;
            else
                if(method(1) == 2)
                    
                    TENS_idx = input('Which EEG/PCS TENS do you want to use? 1/2');

                    emg_prepadding = max([0,t0(1,TENS_idx)-t0(3,1)]);
                    dataOut_prepadding = max([0,t0(3,1)-t0(1,TENS_idx)]);
                    
                    emgDataOut = [zeros(size(emgData,1),emg_prepadding) emgData];
                    dataOut = [zeros(size(dataOut,1),dataOut_prepadding) dataOut];
                    
                    temp = t0(3,:) + emg_prepadding;
                    t0([1 2],:) = t0([1 2],:) + dataOut_prepadding;
                    t0(3,:) = temp;
                    
                    emg_postpadding = max([0,size(dataOut,2)-size(emgDataOut,2)]);
                    dataOut_postpadding = max([0,size(emgDataOut,2)-size(dataOut,2)]);
                    
                    emgDataOut = [emgDataOut zeros(size(emgDataOut,1),emg_postpadding)];
                    dataOut = [dataOut zeros(size(dataOut,1),dataOut_postpadding)];
                    
                    dataOut 	= [dataOut; emgDataOut];
                    dataOut = dataOut(:,1+emg_prepadding:end - emg_postpadding);
                    t0 = t0 - emg_prepadding;
                else
                    try_manualTENS = true;
                    disp('Trying Manual TENS identification');
                    continue;
                end
            end
        end
                  
        % also compute the sample indices for each events with new sampling freq
        for evIdx = 1:numel(eventsInfo)
            if (eventsInfo(evIdx).times >=0)
                eventsInfo(evIdx).samples = round(eventsInfo(evIdx).times * fs);
            else
                eventsInfo(evIdx).samples = round(eventsInfo(evIdx).times * fs) + size(dataOut,2);
            end
        end
        
        pcsHdr.labels(end+1) = {'none'};
        pcsHdr.chanunit(1:2) = {'mV'};
        emgHdr.chanunit(1:emgChannels) = {'mV'};
        
        outHdr.chanunit = [eegHdr.chanunit,pcsHdr.chanunit, emgHdr.chanunit ];
        
        if globDebug
            figure, clf
            subplot(511)
            plot([-100:100]/fs*1000,mean(dataOut(1:eegHdr.NumberOfChannels,(-100:100) + t0(1,1)),1)/max(mean(dataOut(1:eegHdr.NumberOfChannels,(-100:100) + t0(1,1)),1))),
            hold on, plot([-100:100]/fs*1000,dataOut(eegHdr.NumberOfChannels+pcsChIdx,(-100:100) + t0(1,1))/max(dataOut(eegHdr.NumberOfChannels+pcsChIdx,(-100:100) + t0(1,1))),'r');
            plot([-100:100]/fs*1000,dataOut(eegHdr.NumberOfChannels+pcsChannels+emgChIdx,(-100:100) + t0(1,1))/max(dataOut(eegHdr.NumberOfChannels+pcsChannels+emgChIdx,(-100:100) + t0(1,1))),'k');
            if(method(1) == 2)
                subplot(512)
                hold on,plot([-100:100]/fs*1000,mean(dataOut(1:eegHdr.NumberOfChannels,(-100:100) + t0(1,2)),1)/max(mean(dataOut(1:eegHdr.NumberOfChannels,(-100:100) + t0(1,2)),1))),
                plot([-100:100]/fs*1000,dataOut(eegHdr.NumberOfChannels+pcsChIdx,(-100:100) + t0(1,2))/max(dataOut(eegHdr.NumberOfChannels+pcsChIdx,(-100:100) + t0(1,2))),'r');
                plot([-100:100]/fs*1000,dataOut(eegHdr.NumberOfChannels+pcsChannels+emgChIdx,(-100:100) + t0(1,2))/max(dataOut(eegHdr.NumberOfChannels+pcsChannels+emgChIdx,(-100:100) + t0(1,2))),'k');
            end
            
            h1 = subplot(513);
            plot([0:size(dataOut,2)-1]/fs,mean(dataOut(1:eegHdr.NumberOfChannels,:),1)/max(mean(dataOut(1:eegHdr.NumberOfChannels,:),1)));
            hold on, plot((t0(1,1:method(1))-1)/fs,mean(dataOut(1:eegHdr.NumberOfChannels,t0(1,1:method(1))),1)/max(mean(dataOut(1:eegHdr.NumberOfChannels,:),1)),'*g');
            h2 = subplot(514);
            plot([0:size(dataOut,2)-1]/fs,dataOut(eegHdr.NumberOfChannels+pcsChIdx,:)/max(dataOut(eegHdr.NumberOfChannels+pcsChIdx,:)));
            hold on, plot((t0(2,1:method(2))-1)/fs,dataOut(eegHdr.NumberOfChannels+pcsChIdx,t0(2,1:method(2)))/max(dataOut(eegHdr.NumberOfChannels+pcsChIdx,:)),'*g');
            h3 = subplot(515);
            plot([0:size(dataOut,2)-1]/fs,dataOut(eegHdr.NumberOfChannels+pcsChannels+emgChIdx,:)/max(dataOut(eegHdr.NumberOfChannels+pcsChannels+emgChIdx,:)));
            hold on, plot((t0(3,1:method(3))-1)/fs,dataOut(eegHdr.NumberOfChannels+pcsChannels+emgChIdx,t0(3,1:method(3)))/max(dataOut(eegHdr.NumberOfChannels+pcsChannels+emgChIdx,:)),'*g');
            drawnow
            title('Synched data');
            linkaxes([h1 h2 h3],'x');
        end
        
        % update header info
        outHdr.label 		= [eegHdr.label;pcsHdr.labels';  emgHdr.labels];
        outHdr.nChans		= eegHdr.NumberOfChannels + pcsChannels + emgChannels;
        outHdr.NumberOfChannels = outHdr.nChans;
        outHdr.chantype	= [eegHdr.chantype';pcsHdr.chantype;emgHdr.chantype'];
        outHdr.Fs				= fs;
        
        stnPosStruct = struct('type','stn',...
            'labels',pcsHdr.labels,...
            'sph_theta_besa',-134,...
            'sph_phi_besa',-45);
        
        emgPosStruct = struct('type','emg',...
            'labels',emgHdr.labels,...
            'sph_theta_besa',-134,...
            'sph_phi_besa',-45);
        
        outHdr.layout.pos = [eegHdr.layout.pos, stnPosStruct,  emgPosStruct'];
        
        % this is the only supported data format
        outHdr.DataFormat      = 'BINARY';
        outHdr.DataOrientation = 'MULTIPLEXED';
        outHdr.BinaryFormat    = 'IEEE_FLOAT_32';
        
        % no additional calibration needed, since float32
        outHdr.resolution      = ones(size(outHdr.label));
        
        str = input('Do you want to write data (Y), try manual TENS (M) or discard file (N)? Y/M/N','s');
        if(lower(str)=='y')% write data
            [~, fnamePcs, ~] = fileparts(pcsFname);
            filename = strcat(fnamePcs,'EegPcsEmg');
            
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





