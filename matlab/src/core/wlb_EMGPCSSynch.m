function status = wlb_EMGPCSSynch(varargin)
%WLBTENSSYNCH synchronize multiple files in the same vhdr/eeg file
%	FLAG = WLBTENSSYNCH(VARARGIN) I'll edit it when it will be ready

% Edited 2015-06-15 by Gabriele Arnulfo <gabriele.arnulfo@gmail.com>
global globDebug;

p = inputParser;
p.addRequired('pathEmg',@ischar);
p.addRequired('pathPcs',@ischar);
p.addRequired('outdir',@ischar);

p.addOptional('pathEvents','',@ischar)
p.addOptional('fnameFilters',{''},@iscell);

p.addOptional('pcsRefChannel',1,@isnumeric);
p.addOptional('pcsCuttingTime',[0 0],@(x) isnumeric(x) & numel(x)==2 );
p.addOptional('emgCuttingTime',[0 0],@(x) isnumeric(x) & numel(x)==2);
p.addOptional('automaticDetection',true,@islogical);
p.addOptional('findTensPctg',80,@isnumeric);

p.parse(varargin{:});

pathPcs = p.Results.pathPcs;
pathEmg = p.Results.pathEmg;
pathEvents = p.Results.pathEvents;
fnameFilters = p.Results.fnameFilters;
pctg=p.Results.findTensPctg;
% get filenames within each path
pcsFileNames = dir(fullfile(pathPcs,'*pcs.xml'));
emgFileNames = dir(fullfile(pathEmg,'*emg.txt'));

eveFileNames = dir(fullfile(pathEvents,'*event*.csv'));

% check whether we have the same number of files
assert(wlb_checkDataConsistency(pcsFileNames,emgFileNames));

if ~isempty(fnameFilters{1})
    pcsFileNames = wlb_filterFnames(pcsFileNames,fnameFilters);
    emgFileNames = wlb_filterFnames(emgFileNames,fnameFilters);
end

if isempty(pcsFileNames)
    error('Empty Files check directory');
end


for fileIdx = 1 : numel(pcsFileNames)
    
    pcsFname = fullfile(pathPcs,pcsFileNames(fileIdx).name);
    emgFname = fullfile(pathEmg,emgFileNames(fileIdx).name);
    
    % files has been renamed to match
    % subjectStudyNameDrugTrial##Modality.**
    drugCondition = regexp(pcsFname,'_','split');
    
    % pick drug string
    drugCondition = drugCondition{  ~cellfun(@isempty,(regexp(drugCondition,'(off|on)'))) };
    
    trialIdx = cell2mat(regexp(pcsFname,'trial\d+','match'));
    trialIdx = str2double(cell2mat(regexp(trialIdx,'\d+','match')));
    eveFname = eveFileNames(~cellfun(@isempty,regexp({eveFileNames.name},drugCondition)));
    eveFname = fullfile(pathEvents,eveFname.name);
    
    fprintf('Synch:\t%s\n\t%s\n\t%s\n',pcsFname,emgFname,eveFname);
    
    eventsInfo = [];
    if ~isempty(p.Results.pathEvents)
        eventsInfo = wlb_readExternalEventFile( eveFname, trialIdx );
    end
    
    try
        % read pcs header file
        [pcsHdr, pcsData] = wlb_readActivaPC( pcsFname );
        [emgHdr, emgData] = wlb_readEMG_wue( emgFname );
        
        % cut data if needed
        pcsData = pcsData(:,(p.Results.pcsCuttingTime(1)*...
            pcsHdr.SenseChannelConfig.TDSampleRate)+1:end-(p.Results.pcsCuttingTime(2)*...
            pcsHdr.SenseChannelConfig.TDSampleRate));
        emgData = emgData(:,(p.Results.emgCuttingTime(1)*emgHdr.freq)+1:end-...
            (p.Results.emgCuttingTime(2)*emgHdr.freq));
        
        % pick the first channel
        pcsChIdx = p.Results.pcsRefChannel;
        emgChIdx = find(wlb_ismemberWildcards(lower(emgHdr.labels),...
            {'artefa.t_pulse','tens_pulse'})==1 );
        
        if isempty(pcsChIdx ) || isempty(emgChIdx)
            error('Invalid reference channel');
        end
        
        pcsCh = pcsData(pcsChIdx,:);
        emgCh = emgData(emgChIdx,:);
        
        if globDebug
            figure(1000),clf
            hold on
            plot(pcsCh.*1e4)
            plot(emgCh,'r')
            title('original data');
            legend([{'pcs'},{'emg'}])
        end
        
        
        % search for the TENS artefact
        pcsLocs = wlb_findTENSArtefact(pcsData(pcsChIdx,:),pcsHdr.SenseChannelConfig.TDSampleRate,pctg);
        emgLocs = wlb_findTENSArtefact(emgData(emgChIdx,:),emgHdr.freq,pctg);
        
        method = min([length(pcsLocs)/2,length(emgLocs)/2]);
        pcsLocs = pcsLocs(1:method*2);
        emgLocs = emgLocs(1:method*2);
        
        % actually compute t0 for all channels
        dataCell = [{pcsCh},{emgCh}];
        if p.Results.automaticDetection
            t0 = cellfun(@wlb_findTInit,dataCell,{pcsLocs,emgLocs},...
                {pcsHdr.SenseChannelConfig.TDSampleRate,emgHdr.freq},...
                {method, method},'uni',false);
        else
            t0 = cellfun(@wlb_manualTENS,dataCell,{method,method},'uni',false);
        end
        
        if( method == 2 )
            
            t0 = reshape([t0{:}],2,2)';
            % estimate the correct fs for PCS from EMG
            pcsFs = ((t0(1,2)-t0(1,1))* emgHdr.freq) /(t0(2,2)-t0(2,1));
        else
            t0 = [t0{:}]';
            % use default sampling frequency
            t0(:,2) = [length(pcsCh) length(emgCh)];
            if pcsHdr.SenseChannelConfig.TDSampleRate > 422
                pcsFs = 793.65;
            else
                pcsFs = 422;
            end
        end
        
        % downsample pcs data to integer sampling frequency and eeg accordingly
        fs 							= 400;
        [pcsData,~,~] 	= wlb_resampleCascade(pcsData,fs,pcsFs);
        emgData 				= resample(emgData',fs,emgHdr.freq)';
        
        % each t0 row represent eeg,pcs,emg data 
        % we have to recompute the exact point in time after
        % resampling
        t0(2,:) = (round(t0(2,:)/emgHdr.freq*pcsFs));
        t0 			= round(t0/pcsFs*fs);
        
        % also compute the sample indices fo each events with new sampling freq
        for evIdx = 1:numel(eventsInfo)
            
%            if (eventsInfo(evIdx).times >=0)
                eventsInfo(evIdx).samples = round(eventsInfo(evIdx).times * fs)+ min(t0(:,1));
                eventsInfo(evIdx).times	 = eventsInfo(evIdx).times + min(t0(:,1))/fs;
%            else
%                eventsInfo(evIdx).samples = round(eventsInfo(evIdx).times * fs)+ min(t0(:,1));
%                eventsInfo(evIdx).times	 = eventsInfo(evIdx).times + min(t0(:,1))/fs;
%            end
            
        end
        
        
        % t0 has start and end samples foreach channel
        % cut data to have the same number of samples
        pcsDataOut 	= pcsData(:,max(1,t0(1,1)-t0(2,1)):end);
        emgDataOut 	= emgData(:,max(1,t0(2,1)-t0(1,1)):end);
        
        finalSize  	= min([size(pcsDataOut,2),size(emgDataOut,2)]);
        
        pcsDataOut 	= pcsDataOut(:,1:finalSize);
        emgDataOut 	= emgDataOut(:,1:finalSize);
        
        dataOut 	 	= [pcsDataOut;emgDataOut];
        
        pcsChannels	= size(pcsData,1);
        emgChannels = size(emgData,1);
        
        if(pcsChannels == 1)
            dataOut  = [dataOut; zeros(1,size(dataOut,2))];
            pcsChannels = 2;
        end
        
        pcsHdr.labels(end+1) 	= {'none'};
        pcsHdr.chanUnits(pcsChannels) = {'mV'};
        outHdr.chanunit = [pcsHdr.chanUnits, emgHdr.units];
        
        if globDebug
            wndPlot = -100:100;
            figure(666), clf
            subplot(211)
            hold on, plot(pcsData(pcsChIdx,wndPlot + t0(1)).*1e4,'r');
            plot(emgData(emgChIdx,wndPlot + t0(2)),'k');
            legend([{'pcs'},{'emg'}])
            
            if(method == 2)
                subplot(212)
                hold on, plot(pcsData(pcsChIdx,wndPlot + t0(3)).*1e4,'r');
                plot(emgData(emgChIdx,wndPlot+ t0(4)),'k');
            end
            drawnow
            
            figure(999),clf
            hold on
            plot(pcsDataOut(pcsChIdx,:).*1e4)
            plot(emgDataOut(emgChIdx,:),'r')
            title('data synched');
            legend([{'pcs'},{'emg'}])
        end
        
        pcsHdr.labels = pcsHdr.labels([1,2]);
        pcsHdr.chanUnits = pcsHdr.chanUnits([1,2]);
        pcsChannels = 2;
        
        % update header info
        outHdr.label 	= [pcsHdr.labels';emgHdr.labels'];
        outHdr.nChans = pcsChannels + emgChannels;
        outHdr.NumberOfChannels = outHdr.nChans;
        outHdr.Fs			= fs;
        outHdr.chanunit(outHdr.nChans) = {'mV'};
        
        stnPosStruct = struct('type','stn',...
            'labels',pcsHdr.labels,...
            'sph_theta_besa',-134,...
            'sph_phi_besa',-45);
        
        
        emgPosStruct = struct('type','emg',...
            'labels',emgHdr.labels,...
            'sph_theta_besa',-134,...
            'sph_phi_besa',-45);
        
        outHdr.layout.pos = [stnPosStruct, emgPosStruct];
        outHdr.chantype(outHdr.nChans) = {'other'};
        
        % this is the only supported data format
        outHdr.DataFormat      = 'BINARY';
        outHdr.DataOrientation = 'MULTIPLEXED';
        outHdr.BinaryFormat    = 'IEEE_FLOAT_32';
        
        % no additional calibration needed, since float32
        outHdr.resolution      = ones(size(outHdr.label));
        
    catch ME
        fprintf('%s\n',ME.message);
        continue;
    end
    % write data
    [~, fnamePcs, ~] = fileparts(pcsFname);
    filename = strcat(fnamePcs,'Emg');
    
    outHdr.DataFile = strcat(filename,'.eeg');
    outHdr.MarkerFile = '';
    
    write_brainvision_eeg(p.Results.outdir, outHdr, dataOut);
    if(~isempty(eventsInfo))
        outHdr.MarkerFile = strcat(filename,'.vmrk');
        write_brainvision_vmrk(p.Results.outdir, outHdr, eventsInfo);
    end
    
    write_brainvision_vhdr(p.Results.outdir, outHdr);
    
    
end % for files

status = 0;

end % function

