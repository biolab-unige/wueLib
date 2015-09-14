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

		p.addOptional('pcsRefChannel',2,@isnumeric);
		p.addOptional('pcsCuttingTime',0,@isnumeric);
		p.addOptional('emgCuttingTime',0,@isnumeric);
		p.addOptional('automaticDetection',true,@islogical);

		p.parse(varargin{:});

		pathPcs = p.Results.pathPcs;
		pathEmg = p.Results.pathEmg;
		pathEvents = p.Results.pathEvents;
		fnameFilters = p.Results.fnameFilters;

		% get filenames within each path
		pcsFileNames = dir(fullfile(pathPcs,'*pcs.xml'));
		emgFileNames = dir(fullfile(pathEmg,'*emg.txt'));

		eveFileNames = dir(fullfile(pathEvents,'*event*.csv'));

		% check whether we have the same number of files 
		assert(checkDataConsistency(pcsFileNames,emgFileNames));
		
		if ~isempty(fnameFilters{1})
				pcsFileNames = filterFnames(pcsFileNames,fnameFilters);
				emgFileNames = filterFnames(emgFileNames,fnameFilters);
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
				drugCondition = drugCondition{3};

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
						pcsData = pcsData(:,(p.Results.pcsCuttingTime*...
								pcsHdr.SenseChannelConfig.TDSampleRate)+1:end);
						emgData = emgData(:,(p.Results.emgCuttingTime*emgHdr.freq)+1:end);

						% pick the first channel
						pcsChIdx = p.Results.pcsRefChannel;
						emgChIdx = find(ismemberWildcards(lower(emgHdr.labels),...
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
						pcsLocs = findTENSArtefact(pcsData(pcsChIdx,:),pcsHdr.SenseChannelConfig.TDSampleRate);
						emgLocs = findTENSArtefact(emgData(emgChIdx,:),emgHdr.freq);

						method = min([length(pcsLocs)/2,length(emgLocs)/2]);
						pcsLocs = pcsLocs(1:method*2);
						emgLocs = emgLocs(1:method*2);
						
						% actually compute t0 for all channels
						dataCell = [{pcsCh},{emgCh}];
						if p.Results.automaticDetection	
								t0 = cellfun(@findTInit,dataCell,{pcsLocs,emgLocs},...
											{pcsHdr.SenseChannelConfig.TDSampleRate,emgHdr.freq},...
											{method, method},'uni',false);
						else
								t0 = cellfun(@manualTENS,dataCell,{method,method},'uni',false);
						end
							 
						if( method == 2 )				
												
   							t0 = reshape([t0{:}],2,2)';
								% estimate the correct fs for PCS from EMG
								pcsFs = ((t0(1,2)-t0(1,1))* emgHdr.freq) /(t0(2,2)-t0(2,1));
						else
								t0 = [t0{:}]';
								% use default
								t0(:,2) = [length(pcsCh) length(emgCh)];
								if pcsHdr.SenseChannelConfig.TDSampleRate > 422
									pcsFs = 793.65;
								else
									pcsFs = 422;
								end
						end
							
						% downsample pcs data to integer sampling frequency and eeg accordingly
						fs 							= 400;
						[pcsData,~,~] 	= ResampleCascade(pcsData,fs,pcsFs);
						emgData 				= resample(emgData',fs,emgHdr.freq)';
						
						% each t0 row represent eeg,pcs,emg data before
						% we have to recompute the exact point in time after
						% resampling
						t0(2,:) = (round(t0(2,:)/emgHdr.freq*pcsFs));
						t0 			= round(t0/pcsFs*fs);
															
						% also compute the sample indices fo each events with new sampling freq
						for evIdx = 1:numel(eventsInfo)
							
							if (eventsInfo(evIdx).times >=0)
								eventsInfo(evIdx).samples = round(eventsInfo(evIdx).times * fs)+ min(t0(:,1));						
							  eventsInfo(evIdx).times	 = eventsInfo(evIdx).times + min(t0(:,1))/fs;
							else
								eventsInfo(evIdx).samples = round(eventsInfo(evIdx).times * fs)+ max(t0(:,1));
								eventsInfo(evIdx).times	 = eventsInfo(evIdx).times + max(t0(:,1))/fs;
							end

						end


						% t0 has start and end samples foreach channel
						% cut data to have the same number of samples
						pcsDataOut 	= pcsData(:,max(1,t0(1,1)-t0(2,1)):end);
						emgDataOut 	= emgData(:,max(1,t0(2,1)-t0(1,1)):end);
						
						finalSize   	= min([size(pcsDataOut,2),size(emgDataOut,2)]);
																												
						pcsDataOut 	= pcsDataOut(:,1:finalSize);
						emgDataOut 	= emgDataOut(:,1:finalSize);
																												
						dataOut 	 		= [pcsDataOut;emgDataOut];
						
						pcsChannels	= size(pcsData,1);
						emgChannels  = size(emgData,1);
						
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
						outHdr.nChans 	= pcsChannels + emgChannels;
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

		if nMod1Files == nMod2Files || nMatchingFiles == Mod1Files
				bool = true;
		else
				bool = false;
		end
end


function newloc = findTENSArtefact(data, fs)
%FINDTENSARTEFACT data [1xN] time samples
%	LOCS = FINDTENSARTEFACT(DATA) Long description
	  global globDebug	
		if globDebug
				figure, plot(data.*1e-2), hold on;
		end
		data = abs(data);

		dataBp = wlb_bandpass_fft(data, fs, 90, 110,1,1,[]);
		dataBp(:,:,2) = wlb_bandpass_fft(data, fs, 190, 200, 1,1,[]);
        
		if fs > 700			
            
			dataBp(:,:,3) = wlb_bandpass_fft(data, fs, 290, 310, 1,1,[]);
			
		end
		dataBp = wlb_bandpass_fft(mean(abs(dataBp(:,:,:)),3), fs, 0.001, 1, 1,1,[]);

    pctg = 80;
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
		
 		if(length(locs)/numel(pks)~=2) %ci sono piu di 2 intersezioni in uno o entrambi i picchi
				for i=1:numel(pks)
						loc = locs( abs(locs-pksLocs(i))<5*fs ); %cerco le intersezioni piu vicine al picco
						startLoc = find(loc-pksLocs(i)<0,1,'first');
						endLoc = find(loc-pksLocs(i)>0,1,'last');
						loc = [loc(startLoc) loc(endLoc)];
						if ~isempty( loc )
							newloc = [newloc loc];
						end
				end
 		else
 				newloc = locs;
 		end
          
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
	imf = ceemdan(D,0.0002,20,100,1);
	plot(bsxfun(@plus,imf,max(imf(:))*(1:(size(imf,1)))')')
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
