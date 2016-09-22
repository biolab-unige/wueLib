function status = wlb_EMGPCSSynch(varargin)
%BN_TENSSYNCH synchronize multiple files in the same vhdr/eeg file
%	FLAG = BN_TENSSYNCH(VARARGIN) I'll edit it when it will be ready
%

% Edited 2014-09-22 by Gabriele Arnulfo <gabriele.arnulfo@gmail.com>
%
%%% DEFINE SUPPORTED FILE EXTENSIONS %%%%
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
    p.addOptional('find_tens_pctg',80,@isnumeric);

		p.parse(varargin{:});

		pathPcs = p.Results.pathPcs;
		pathEmg = p.Results.pathEmg;
		pathEvents = p.Results.pathEvents;
		fnameFilters = p.Results.fnameFilters;
        pctg=p.Results.find_tens_pctg;
		% get filenames within each path
		pcsFileNames = dir(fullfile(pathPcs,'*pcs.xml'));
		emgFileNames = dir(fullfile(pathEmg,'*emg.*'));
		eveFileNames = dir(fullfile(pathEvents,'*event*.csv'));

		assert(checkDataConsistency(pcsFileNames,emgFileNames),...
				'Number of emg files and pcs files do not match');

		if ~isempty(fnameFilters{1})
				pcsFileNames = filterFnames(pcsFileNames,fnameFilters);
				emgFileNames = filterFnames(emgFileNames,fnameFilters);
		end

		for fileIdx =1:numel(pcsFileNames)

                % files has been renamed to match 
				% subject_studyName_drug_trial##_modality.**
				drugCondition = regexp(pcsFileNames(fileIdx).name,'_','split');

				% pick drug string
				drugCondition = drugCondition{3};

                pcsFname = fullfile(pathPcs,pcsFileNames(fileIdx).name);
				emgFname = fullfile(pathEmg,emgFileNames(fileIdx).name);

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

				fprintf('Synch:\t%s\n\t%s\n\t%s\n\t%s\n',pcsFname,emgFname,eveFname);

				eventsInfo = [];
				if ~isempty(p.Results.pathEvents)
					eventsInfo = wlb_readExternalEventFile( eveFname, trialIdx );
                end
               
                %				try
                % read pcs header file
                [pcsHdr, pcsData] = wlb_readActivaPC( pcsFname );
                [emgHdr, emgData] = wlb_readEMG_wue( emgFname );

                pcsData = cutInitialSamplesData(pcsData,[],[],p.Results.pcsCuttingTime,pcsHdr.SenseChannelConfig.TDSampleRate);
                emgData = cutInitialSamplesData(emgData,[],[],p.Results.emgCuttingTime,emgHdr.freq);
                
                % define TENS channels
                pcsChIdx = p.Results.pcsRefChannel;
                emgChIdx = find(ismemberWildcards(lower(emgHdr.labels),...
                    {'artefa.t_pulse','tens_pulse'})==1 );

                if isempty(pcsChIdx ) || isempty(emgChIdx)
                    error('Invalid reference channel');
                end

                pcsCh = pcsData(pcsChIdx,:);
                emgCh = emgData(emgChIdx,:);
                
                % pick central window
                pcsWnd = round(numel(pcsCh)/2)+[-1/2*pcsHdr.SenseChannelConfig.TDSampleRate:...
                    1/2*pcsHdr.SenseChannelConfig.TDSampleRate];
                emgWnd = round(numel(emgCh)/2)+[-1/2*emgHdr.freq:1/2*emgHdr.freq];
                
                pcsStd = std(pcsCh(pcsWnd));
                emgStd = std(emgCh(emgWnd));
                
                % sometimes EMG or PCS data have been cutted at
                % the begginig or at the end, we then pad with
                % random samples around the mean
                padValsPcs = pcsStd.*randn(1,4*pcsHdr.SenseChannelConfig.TDSampleRate);
                padValsEmg = emgStd.*randn(1,4*emgHdr.freq);
                
                pcsCh = [padValsPcs, pcsCh, padValsPcs];
                emgCh = [padValsEmg, emgCh, padValsEmg];
                
                % find time windows containing TENS
                pcsLocs = findTENSArtefact(pcsCh,pcsHdr.SenseChannelConfig.TDSampleRate,pctg);
                emgLocs = findTENSArtefact(emgCh,emgHdr.freq,pctg);

                method = min([length(emgLocs)/2,length(pcsLocs)/2]);
                
                % actually compute t0 for all channels
                dataCell = [{emgCh},{pcsCh}];
                
                if p.Results.automaticDetection
                    t0 = cellfun(@findTInit,dataCell,{emgLocs,pcsLocs},...
                        {emgHdr.freq,pcsHdr.SenseChannelConfig.TDSampleRate},...
                        {method, method},'uni',false);
                else
                    t0 = cellfun(@manualTENS,dataCell,{method,method},'uni',false);
                end
                
                
                if( method == 2 )
                    t0 = reshape([t0{:}],2,2)';
                    % estimate the correct fs for PCS
                    pcsFs = (t0(2,2)-t0(2,1))* emgHdr.freq /(t0(1,2)-t0(1,1));
                else
                    t0 = [t0{:}]';
                    t0(:,2) = [length(emgCh) length(pcsCh)];
                    if pcsHdr.SenseChannelConfig.TDSampleRate > 422
                        pcsFs = 793.65;
                    else
                        pcsFs = 422;
                    end
                end
                
                % remove the padded samples
                offsets = [emgHdr.freq*4; pcsHdr.SenseChannelConfig.TDSampleRate*4];
                t0 = t0 - repmat(offsets,[1 2]);
                
                
                % downsample pcs data to integer sampling frequency and eeg accordingly
                fs 				= 400;
                [pcsData,~,~]   = ResampleCascade(pcsData,fs,pcsFs);
                emgData 		= resample(emgData',fs,emgHdr.freq)';
                
                % each t0 row represent eeg,pcs data before
                % we have to recompute the exact point in time after
                % resampling
                t0(1,:) = (round(t0(1,:)/emgHdr.freq*pcsFs));
                t0 		= round(t0/pcsFs*fs);
                
                % t0 has start and end samples foreach channel
                % cut data to have the same number of samples
                emgDataOut 	= emgData(:,max(1,t0(1,1)-t0(2,1)):end);
                pcsDataOut 	= pcsData(:,max(1,t0(2,1)-t0(1,1)):end);
                
                finalSize   	= min([size(pcsDataOut,2),size(emgDataOut,2)]);
                
                pcsDataOut 	= pcsDataOut(:,1:finalSize);
                emgDataOut 	= emgDataOut(:,1:finalSize);
                
                dataOut 	 		= [pcsDataOut;emgDataOut];
                
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

                                
                % also compute the sample indices for each events with new sampling freq
                for evIdx = 1:numel(eventsInfo)
                    if (eventsInfo(evIdx).times >=0)
                        eventsInfo(evIdx).samples = round(eventsInfo(evIdx).times * fs)+ min(t0(:,1));
                    else
                        eventsInfo(evIdx).samples = round(eventsInfo(evIdx).times * fs)+ max(t0(:,1));
                    end
                end
                                
                pcsHdr.labels(end+1) = {'none'};
                pcsHdr.chanunit(1:2) = {'mV'};
                emgHdr.chanunit(1:emgChannels) = {'mV'};
                
                outHdr.chanunit = [pcsHdr.chanunit, emgHdr.chanunit ];
                
                if globDebug
                    figure, clf
                    subplot(411)
                    plot([-100:100]/400*1000,emgData(emgChIdx,(-100:100) + t0(1))),
                    hold on, plot([-100:100]/400*1000,100*pcsData(pcsChIdx,(-100:100) + t0(2)).*100,'r');
                    
                    if(method == 2)
                        subplot(412)
                        plot([-100:100]/400*1000,emgData(emgChIdx,(-100:100) + t0(3))),
                        hold on, plot([-100:100]/400*1000,100*pcsData(pcsChIdx,(-100:100) + t0(4)).*100,'r');
                    end
                    subplot(413),
                    plot([0:size(dataOut,2)-1]/400,dataOut([pcsChIdx],:)');
                    subplot(414),
                    plot([0:size(dataOut,2)-1]/400,dataOut(pcsChannels +emgChIdx,:)');
                    drawnow
                    title('Synched data');
                end
              
                % update header info
                outHdr.label 		= [pcsHdr.labels';  emgHdr.labels'];
                outHdr.nChans		= pcsChannels + emgChannels;
                outHdr.NumberOfChannels = outHdr.nChans;
                outHdr.chantype	= [pcsHdr.chantype;emgHdr.chantype'];
                outHdr.Fs				= fs;
                
                stnPosStruct = struct('type','stn',...
                    'labels',pcsHdr.labels,...
                    'sph_theta_besa',-134,...
                    'sph_phi_besa',-45);
                
                emgPosStruct = struct('type','emg',...
                    'labels',emgHdr.labels,...
                    'sph_theta_besa',-134,...
                    'sph_phi_besa',-45);
                
                outHdr.layout.pos = [stnPosStruct,  emgPosStruct];
                
                % this is the only supported data format
                outHdr.DataFormat      = 'BINARY';
                outHdr.DataOrientation = 'MULTIPLEXED';
                outHdr.BinaryFormat    = 'IEEE_FLOAT_32';
                
                % no additional calibration needed, since float32
                outHdr.resolution      = ones(size(outHdr.label));
                
                %				catch ME
                %						fprintf('%s\n',ME.message);
                %						continue;
                %				end
                
                % write data
                [~, fnamePcs, ~] = fileparts(pcsFname);
                filename = strcat(fnamePcs,'EmgPcs');
                
                outHdr.DataFile = strcat(filename,'.eeg');
                outHdr.MarkerFile = '';
                
                write_brainvision_eeg(p.Results.outdir, outHdr, dataOut);
                
                eventsInfo = rmfield(eventsInfo,'times');
                eventsInfo = [eventsInfo];
                
                if(~isempty(eventsInfo))
                    outHdr.MarkerFile = strcat(filename,'.vmrk');
                    write_brainvision_vmrk(p.Results.outdir, outHdr, eventsInfo);
                end
                
                outHdr.good_chan = ones(outHdr.NumberOfChannels,1);
                write_brainvision_vhdr(p.Results.outdir, outHdr);
                
        end % for files
        
        status = 0;


end % function

function tau = findTInit(D,locs,~,chunks)
%FINDTINIT Description
%	TAU = FINDTINIT(D,LOCS,CHUNKS) Long description
%

% number of levels used for WICA and wavelet family
numLvl    = 8;
vfilter   = 'db1';
tau 			= nan(2,1);


for chunk = 0:chunks-1
    
    artefactIdx	= round(locs(1+(2*chunk))):round(locs(2+(2*chunk))) ;
    
    artefactIdx(artefactIdx > numel(D)) = [];
    artefactIdx(artefactIdx < 1) = [];
    data 				= D(artefactIdx);
    
    nsamples    = length(data);
    offset      = ceil(nsamples/(2^numLvl)) * (2^numLvl);
    
    data(nsamples+1:offset) = zeros(1,offset-nsamples);
    
    
    [thr] = 0.5*ddencmp('den','wv',data);
    
    [swa, swd] = swt(data,numLvl, vfilter);
    
    swd(abs(swd) < thr) = 0;
    swd(3:numLvl,:) = zeros(size(swd(3:numLvl,:)));
    swa = zeros(size(swa));
    
    out = iswt(swa,swd,vfilter)';
    out = out(1:nsamples);
    thresh = 2*std(abs(out));
    
    out( abs(out) <  std(abs(out)) ) = 0;
    
    [~, maxLocs] = findpeaks(abs((out)),'MINPEAKHEIGHT',thresh);
    dataLocs = max(maxLocs);
    
    tau(chunk+1) = dataLocs + artefactIdx(1)-1;
    
end
tau(isnan(tau)) = [];
figure,
plot(D) , hold on, plot(tau,D(tau),'rx');

end

function [x,Pfac,Qfac] = ResampleCascade(x,NewRate,OldRate,Method)
    % Default method: 'resample'
    if (nargin < 4)
        Method = 'resample';
    end
    % Common factors
    [P,Q] = rat(NewRate/OldRate);
    % We want to upsample by P and downsample by Q to achieve the new rate
    % But big numbers cause problems.
    Pfac = factor(P);
    Qfac = factor(Q);
    % Longest number of factors
    iFacs = max(length(Pfac),length(Qfac));
    % Pad the shorter one to have unity factors
    Pfac((length(Pfac)+1):iFacs) = 1;
    Qfac((length(Qfac)+1):iFacs) = 1;

    % So now we have two factorization lists of the same length, and
    % prod(Pfac) / prod(Qfac) = P/Q.
    Pfac = sort(Pfac,'descend'); % upsample largest first
    Qfac = sort(Qfac,'ascend'); % downsample smallest rates first
    Rates = Pfac./Qfac;  % rates per step
    CRate = cumprod(Rates); % cumulative resampling rates

    % We can't go below min(1,P/Q) without losing information. Because of low-pass filtering, 
		% don't be too precise
    Problem = CRate < (0.9 * P/Q);
    if any(Problem)
        fprintf(1, 'RESAMPLE> Warning: Desired rate is %.f\n', P/Q);
    end
    if any(Pfac > 10)
        disp(['RESAMPLE> Warning: Upsampling by more than 10 in the cascades, P = ' sprintf('%d ', Pfac)]);
    end
    if any(Qfac > 10)
        disp(['RESAMPLE> Warning: Downsampling by more than 10 in the cascades, Q = ' sprintf('%d ', Qfac)]);
    end

    % ===== RESAMPLING =====
    switch Method
        % Decimate/interp inputs cannot be vectorized
        case 'decimate'
            % Initialize output parameters
            lenResmp = ceil(size(x,2) * prod(Pfac) / prod(Qfac));
            nRow = size(x,1);
            xResmp = zeros(nRow, lenResmp);
            % Loop on factors and rows
            for iRow = 1:size(x,1)
                xTmp = x(iRow,:);
                for i = 1:iFacs
                    xTmp = decimate(interp(xTmp, Pfac(i)), Qfac(i));
                end
                xResmp(iRow,:) = xTmp;
            end
            x = xResmp;
        % Resample takes vectorized inputs
        case 'resample'
            for i = 1:iFacs
                x = resample(x', Pfac(i), Qfac(i))';
            end
    end
end

function fnames = filterFnames(fnames,pattern)
%FILTERFNAMES Description
%	FNAME = FILTERFNAMES(FNAMES,PATTERN) Long description
%
		tmp = {fnames.name};
		mask= zeros(numel(tmp),numel(pattern));

		for el = 1:numel(tmp)
			mask(el,:) = ~cellfun(@isempty,regexp(tmp(el),pattern));
		end

		mask = logical(prod(mask,2));

		fnames = fnames(mask);
end

function bool = checkDataConsistency(fnameMod1, fnameMod2)
%CHECKDATACONSISTENCY Description
%	BOOL = CHECKDATACONSISTENCY(FNAMEMOD1, FNAMEMOD2) Long description
%

		fnameMod1 = {fnameMod1.name};
		fnameMod2 = {fnameMod2.name};

		[~,mod1,~] = cellfun(@fileparts,fnameMod1,'uni',false);
		[~,mod2,~] = cellfun(@fileparts,fnameMod2,'uni',false);

		nMod1Files = numel(mod1);
		nMod2Files = numel(mod2);

		nMatchingFiles = sum(ismember(mod1,mod2));

		if nMod1Files == nMod2Files || nMatchingFiles == nMod1Files
				bool = true;
		else
				bool = false;
		end
end


function newloc = findTENSArtefact(data, fs,pctg)
%FINDTENSARTEFACT data [1xN] time samples
%	LOCS = FINDTENSARTEFACT(DATA) Long description
	  global globDebug	
		if globDebug
				figure, plot(data.*1e-2), hold on;
                title('findTENSArtefact');
		end
		data = abs(data);

		dataBp = wlb_bandpass_fft(data, fs, 90, 110,1,1,[]);
		dataBp(:,:,2) = wlb_bandpass_fft(data, fs, 190, 200, 1,1,[]);
        
		if fs > 700			
            
			dataBp(:,:,3) = wlb_bandpass_fft(data, fs, 290, 310, 1,1,[]);
			
		end
		dataBp = wlb_bandpass_fft(mean(abs(dataBp(:,:,:)),3), fs, 0.001, 1, 1,1,[]);

		thr = prctile(dataBp,pctg);
		while(thr<0)
				pctg = pctg +1;
				thr = prctile(dataBp,pctg);        
		end
        
		[pks,pksLocs] = findpeaks(dataBp,'MINPEAKHEIGHT',thr,'MINPEAKDISTANCE',10*fs);
		if(length(pksLocs)>2)
				% in general we should have only two tens 
				% but we might end up with some weird artefacts


				[val,I] = sort(pks,'descend');
				pks 		= val(1:2);
				pksLocs	= pksLocs(I);
				pksLocs	= pksLocs(1:2);

		end

		[locs] = find(abs(diff((dataBp>thr))));
		dataBpDer = -savitzkyGolayFilt(dataBp,4,1,11);
		locsDer = dataBpDer(locs);
		
		if(locsDer(1)<0),locs(1)=[];end
		if(locsDer(end)>0),locs(end)=[];end

		if globDebug
				plot(dataBp,'r');
				plot(locs,dataBp(locs),'ko');
		end

	 
		locs = reshape(locs,2,numel(locs)/2);
		
		duration = diff(locs);
		
		locs = locs(:,duration>0.8*fs);
        
		newloc = [];
		
%  		if(length(locs)/numel(pks)~=2) %ci sono piu di 2 intersezioni in uno o entrambi i picchi
				for i=1:numel(pks)
						loc = locs( abs(locs-pksLocs(i))<6*fs ); %cerco le intersezioni piu vicine al picco
						startLoc = find(loc-pksLocs(i)<0,1,'first');
						endLoc = find(loc-pksLocs(i)>0,1,'last');
						loc = [loc(startLoc) loc(endLoc)];
						if ~isempty( loc )
							newloc = [newloc loc];
						end
 				end
%  		else
%  				newloc = locs;
%  		end
%           
		newloc = sort(newloc,'ascend');

		if globDebug
				plot(newloc,dataBp(newloc),'kx','MarkerSize',5);
		end


end

function tau = manualTENS(D,method)
%MANUALTENS Description
%	TAU = MANUALTENS(D,locs) Long description
%
	warnMessage = sprintf('You should pick only %d point(s)',method);
	warndlg(warnMessage);
	f1 = figure;
% 	imf = ceemdan(D,0.0002,20,100,1);
	plot(D)%bsxfun(@plus,imf,max(imf(:))*(1:(size(imf,1)))')')
	addCrossair(f1);

	waitfor(f1)

	tau = evalin('base','cursorValue');
	

end


function mask = ismemberWildcards(stringsIn, patterns)
%ISMEMBERWILDCARDS Description
%	MASK = ISMEMBERWILDCARDS(STRINGSIN, PATTERNS) Long description
%
	
	nStrings = numel(stringsIn);
	nPatterns= numel(patterns);

	mask = cellfun(@(x) regexp(x,patterns),stringsIn,'Uni',false);
	mask = [mask{:}];
	mask = reshape(mask,[nPatterns,nStrings]);
	mask(cellfun(@isempty,mask)) = {0};

	mask = cell2mat(mask);

	mask = sum(mask) >= 1;
end

function [data,hdr,event] = cutInitialSamplesData(data,hdr,event,offset,fs)
%CUTINITIALSAMPLESDATA Description
%	DATA = CUTINITIALSAMPLESDATA(DATA,OFFSET,FS) Long description
%

data = data(:,round(offset(1)*fs)+1:end-round(offset(2)*fs));
if(not(isempty(hdr)))
    hdr.nSamples = size(data,2);
end
if(not(isempty(event)))
    for i=1:numel(event)
        event(i).samples = event(i).samples - round(offset(1)*fs);
    end
end
end
