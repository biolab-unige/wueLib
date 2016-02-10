function status = wlb_EEGEMGPCSSynch(varargin)
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
		p.addRequired('pathHdeeg',@ischar);
		p.addRequired('outdir',@ischar);

		p.addOptional('pathEvents','',@ischar)
		p.addOptional('fnameFilters',{''},@iscell);

    p.addOptional('pcsRefChannel',1,@isnumeric);
    p.addOptional('pcsCuttingTime',[0 0],@(x) isnumeric(x) & numel(x)==2 );
    p.addOptional('eegCuttingTime',[0 0],@(x) isnumeric(x) & numel(x)==2);
    p.addOptional('emgCuttingTime',[0 0],@(x) isnumeric(x) & numel(x)==2);
    p.addOptional('automaticDetection',true,@islogical);
    p.addOptional('findTensPctg',80,@isnumeric);

		p.parse(varargin{:});

		pathPcs = p.Results.pathPcs;
		pathEmg = p.Results.pathEmg;
		pathHdeeg = p.Results.pathHdeeg;
		pathEvents = p.Results.pathEvents;
		fnameFilters = p.Results.fnameFilters;
    pctg=p.Results.findTensPctg;
		% get filenames within each path
		pcsFileNames = dir(fullfile(pathPcs,'*pcs.xml'));
		emgFileNames = dir(fullfile(pathEmg,'*emg.txt'));
		eegFileNames = dir(fullfile(pathHdeeg,'*eeg.eeg'));

		eveFileNames = dir(fullfile(pathEvents,'*event*.csv'));

		assert(checkDataConsistency(pcsFileNames,emgFileNames),...
				'Number of emg files and pcs files do not match');

		assert(checkDataConsistency(pcsFileNames,eegFileNames),...
				'Number of hdeeg files and pcs files do not match');


		if ~isempty(fnameFilters{1})
				pcsFileNames = wlb_filterFnames(pcsFileNames,fnameFilters);
				emgFileNames = wlb_filterFnames(emgFileNames,fnameFilters);
				eegFileNames = wlb_filterFnames(eegFileNames,fnameFilters);
		end

		for fileIdx =1:numel(pcsFileNames)

				pcsFname = fullfile(pathPcs,pcsFileNames(fileIdx).name);
				emgFname = fullfile(pathEmg,emgFileNames(fileIdx).name);
				eegFname = fullfile(pathHdeeg,eegFileNames(fileIdx).name);

				% files has been renamed to match 
				% subject_studyName_drug_trial##_modality.**
				drugCondition = regexp(pcsFname,'_','split');

				% pick drug string
    		drugCondition = drugCondition{  ~cellfun(@isempty,(regexp(drugCondition,'(off|on)'))) };

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

				eventsInfo = [];
				if ~isempty(p.Results.pathEvents)
					eventsInfo = wlb_readExternalEventFile( eveFname, trialIdx );
				end

				try
						% read pcs header file
						[pcsHdr, pcsData] = wlb_readActivaPC( pcsFname );
						[emgHdr, emgData] = wlb_readEMG_wue( emgFname );
						[eegHdr, eegData] = wlb_readBrainvision( eegFname );
					 
						% cut data if needed
						pcsData = wlb_cutInitialSamplesData(pcsData,p.Results.pcsCuttingTime,eegHdr.SenseChannelConfig.TDSampleRate);
						emgData = wlb_cutInitialSamplesData(emgData,p.Results.emgCuttingTime,eegHdr.freq);
						eegData = wlb_cutInitialSamplesData(eegData,p.Results.eegCuttingTime,eegHdr.SamplingInterval);
                        
                        

						% define TENS channels
						pcsChIdx = p.Results.pcsRefChannel;
						emgChIdx = find(wlb_ismemberWildcards(lower(emgHdr.labels),{'art.fa.t_pulse','tens_pulse'})==1 );
						eegChIdx = 1;

						if isempty(pcsChIdx ) || isempty(emgChIdx) || isempty(eegChIdx)
								error('Invalid reference channel');
						end	
					
						% extract TENS channels
						pcsCh = pcsData(pcsChIdx,:);
						emgCh = emgData(emgChIdx,:);
						eegCh = eegData(eegChIdx,:);
						
						% pick central window
						pcsWnd = round(numel(pcsCh)/2)-[1/2*pcsHdr.SenseChannelConfig.TDSampleRate:...
								1/2*pcsHdr.SenseChannelConfig.TDSampleRate];
						emgWnd = round(numel(emgCh)/2)-[1/2*emgHdr.freq:1/2*emgHdr.freq];
						eegWnd = round(numel(eegCh)/2)-[1/2*eegHdr.SamplingInterval:1/2*eegHdr.SamplingInterval];
						
						pcsStd = std(pcsCh(pcsWnd));
						emgStd = std(emgCh(emgWnd));
						eegStd = std(eegCh(eegWnd));
						
						% sometimes EMG or PCS data have been cutted at
						% the begginig or at the end, we then pad with
						% random samples around the mean
						padValsPcs = pcsStd.*randn(1,4*pcsHdr.SenseChannelConfig.TDSampleRate);
						padValsEmg = emgStd.*randn(1,4*emgHdr.freq);
						padValsEeg = eegStd.*randn(1,4*eegHdr.SamplingInterval);
						
						pcsCh = [padValsPcs, pcsCh, padValsPcs];
						emgCh = [padValsEmg, emgCh, padValsEmg];
						eegCh = [padValsEeg, eegCh, padValsEeg];
			
						% find time windows containing TENS
						pcsLocs = wlb_findTENSArtefact(pcsCh,pcsHdr.SenseChannelConfig.TDSampleRate,pctg);
						eegLocs = wlb_findTENSArtefact(eegCh,eegHdr.SamplingInterval,pctg);
						emgLocs = wlb_findTENSArtefact(emgCh,emgHdr.freq,pctg);
						
						method = min([length(eegLocs)/2,length(pcsLocs)/2,length(emgLocs)/2]);
						
						% actually compute t0 for all channels
						dataCell = [{eegCh},{pcsCh},{emgCh}];
						
						if p.Results.automaticDetection	
								t0 = cellfun(@wlb_findTInit,dataCell,{eegLocs,pcsLocs,emgLocs},...
										{eegHdr.Fs,pcsHdr.SenseChannelConfig.TDSampleRate,emgHdr.freq},...
										{method, method,method},'UniformOutput',false);
						else
								t0 = cellfun(@manualTENS,dataCell,{method,method,method},'UniformOutput',false);
						end

										
						if( method == 2 )
								t0 = reshape([t0{:}],2,3)';
								% estimate the correct fs for PCS
								pcsFs = (t0(2,2)-t0(2,1))* eegHdr.Fs /(t0(1,2)-t0(1,1));
						else
								t0 = [t0{:}]';
								t0(:,2) = [length(eegCh) length(pcsCh) length(emgCh)];
								if pcsHdr.SenseChannelConfig.TDSampleRate > 422
										pcsFs = 793.65;
								else
										pcsFs = 422;
								end
						end
						
						% remove the padded samples
						offsets = [eegHdr.SamplingInterval*4; pcsHdr.SenseChannelConfig.TDSampleRate*4;...
												emgHdr.freq*4 ];
						t0 = t0 - repmat(offsets,[1 2]);
					
						
						% downsample pcs data to integer sampling frequency and eeg accordingly
						fs 				= 400;
						[pcsData,~,~]   = wlb_resampleCascade(pcsData,fs,pcsFs);
						eegData 		= resample(eegData',fs,eegHdr.Fs)';
						emgData 		= resample(emgData',fs,emgHdr.freq)';
						
						% each t0 row represent eeg,pcs,emg data before
						% we have to recompute the exact point in time after
						% resampling
						t0(1,:) = (round(t0(1,:)/eegHdr.Fs*pcsFs));
						t0(3,:) = (round(t0(3,:)/emgHdr.freq*pcsFs));
						t0 		= round(t0/pcsFs*fs);
									 
						% also compute the sample indices for each events with new sampling freq
						for evIdx = 1:numel(eventsInfo)
							
							if (eventsInfo(evIdx).times >=0)
								eventsInfo(evIdx).samples = round(eventsInfo(evIdx).times * fs)+ min(t0(:,1));						
                                eventsInfo(evIdx).times	 = eventsInfo(evIdx).times + min(t0(:,1))/fs;

							else
								eventsInfo(evIdx).samples = round(eventsInfo(evIdx).times * fs)+ max(t0(:,1));
								eventsInfo(evIdx).times	 = eventsInfo(evIdx).times + max(t0(:,1))/fs;
							end

						end


						% t0 has start and end samples for each channel
						% cut data to have the same number of samples
						% We compute the smallest offset possible in the initial
						% samples prior each TENS artefact
						eegDataOut 	= eegData(:,max([1,t0(1,1)-t0(2,1),t0(1,1)-t0(3,1)]):end);
						pcsDataOut 	= pcsData(:,max([1,t0(2,1)-t0(1,1),t0(2,1)-t0(3,1)]):end);
						emgDataOut 	= emgData(:,max([1,t0(3,1)-t0(1,1),t0(3,1)-t0(2,1)]):end);
						
						% the we determine the minimun data length
						finalSize  	= min([size(eegDataOut,2),...
                      size(pcsDataOut,2),....
											size(emgDataOut,2)]);
						% cut accordingly each modality and pack everything
						% in a single data matrix
						dataOut 	= [eegDataOut(:,1:finalSize); ...
										emgDataOut(:,1:finalSize);
                    pcsDataOut(:,1:finalSize)];
						
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
						
						pcsHdr.labels(end+1) = {'none'};
						pcsHdr.chanunit(1:2) = {'mV'};
						emgHdr.chanunit(1:emgChannels) = {'mV'};

						outHdr.chanunit = [eegHdr.chanunit, emgHdr.chanunit, pcsHdr.chanunit];

						if globDebug	
								figure, clf
								subplot(311)
								plot(eegData(eegChIdx,(-100:100) + t0(1))),
								hold on, plot(pcsData(pcsChIdx,(-100:100) + t0(2)).*10,'r');
								plot(emgData(emgChIdx,(-100:100) + t0(3)),'k');

								if(method == 2)
										subplot(312)
										plot(eegData(eegChIdx,(-100:100) + t0(4))),
										hold on, plot(pcsData(pcsChIdx,(100:100) + t0(5)).*10,'r');
										plot(emgData(emgChIdx,(-100:100) + t0(6)),'k');
								end
								subplot(313),
								plot(dataOut([eegChIdx end-1:end],:)');
								drawnow
                title('Synched data');
						end
						
%						if(eventTimeoffset)
%								mtmdataFile=fullfile(eventOffsetDir, [strjoin({fnameEegParts{2},...
%										fnameEegParts{end-1}},'_'),'.txt']);
%
%								[out, delimiter, headerlines] = importdata(mtmdataFile);
%								offset = out.data(:,4);
%								iMrk = find(strcmpi({eegEvent.label}, eventMarkerLabel));
%								if(length(offset)== size(eegEvent(iMrk).samples,2))
%										for offsetIdx = 1:length(offset)
%										eegEvent(iMrk).samples(:,offsetIdx)=...
%												round(eegEvent(iMrk).samples(:,offsetIdx)/...
%												eegHdr.Fs*fs)+ repmat(round(offset(offsetIdx)*1e-3*fs),...
%												size(eegEvent(iMrk).samples,1),1)-dataWnd(1,1)...
%												-eegCuttingtime*fs;
%										end
%								else 
%										error('Offset length differs from event length');
%								end
%						end
						
						% update header info
						outHdr.label 		= [eegHdr.label;  emgHdr.labels';pcsHdr.labels'];
						outHdr.nChans		= eegHdr.NumberOfChannels + pcsChannels + emgChannels;
						outHdr.NumberOfChannels = outHdr.nChans;
						outHdr.chantype	= eegHdr.chantype;
						outHdr.Fs				= fs;
						
						stnPosStruct = struct('type','stn',...
										'labels',pcsHdr.labels,...
										'sph_theta_besa',-134,...
										'sph_phi_besa',-45);

						emgPosStruct = struct('type','emg',...
										'labels',emgHdr.labels,...
										'sph_theta_besa',-134,...
										'sph_phi_besa',-45);
						
						outHdr.layout.pos = [eegHdr.layout.pos,  emgPosStruct, stnPosStruct];
						outHdr.chantype(end:outHdr.nChans) = {'other'};
						
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
				filename = strcat(fnamePcs,'EmgEeg');
				
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
