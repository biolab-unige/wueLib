function flag = bn_EEGEMGPCSSynch(path_pcs,path_hdeeg,path_emg)
%BN_TENSSYNCH synchronize multiple files in the same vhdr/eeg file
%	FLAG = BN_TENSSYNCH(VARARGIN) I'll edit it when it will be ready
%

% Edited 2014-09-22 by Gabriele Arnulfo <gabriele.arnulfo@gmail.com>
%
%%% DEFINE SUPPORTED FILE EXTENSIONS %%%%
ext_hdeeg 				= '*.eeg';
ext_pcs 					= '*.xml';
ext_emg 					= '*.txt';
event_timeoffset 	= 0;

fname_pcs 				= dir(fullfile(path_pcs,ext_pcs));
fname_hdeeg 			= dir(fullfile(path_hdeeg,ext_hdeeg));
fname_emg		 			= dir(fullfile(path_emg,ext_emg));

assert(checkDataConsistency(fname_hdeeg,fname_pcs),...
		'Number of hdeeg files and pcs files do not match');

assert(checkDataConsistency(fname_hdeeg,fname_emg),...
		'Number of hdeeg files and emg files do not match');

fname_hdeeg = filterFnames(fname_hdeeg,'reachandgrasp');
fname_hdeeg = filterFnames(fname_hdeeg,'trial2');

for i =1:length(fname_hdeeg)
    disp([num2str(i) ' - ' fname_hdeeg(i).name])

    % read eeg header file
    fname_eeg_parts = strsplit(fname_hdeeg(i).name,{'_','.'});
    fname_pcs = [strjoin([fname_eeg_parts([1 2]) {'OFF'}...
        fname_eeg_parts([4:end-1])],'_'), '.xml'];

    fname_emg = [strjoin([fname_eeg_parts([1 2]) {'OFF'}...
        fname_eeg_parts([4:end-1])],'_'), '.txt'];
        
    eeg_hdr  = read_brainvision_vhdr(fullfile(path_hdeeg,...
				strcat(strtok(fname_hdeeg(i).name,'.'),'.vhdr')));
    eeg_data = read_brainvision_eeg(fullfile(path_hdeeg,...
				strcat(strtok(fname_hdeeg(i).name,'.'),'.eeg')),...
        eeg_hdr,0,eeg_hdr.nSamples);

		[p f e] = fileparts(fname_hdeeg(i).name);

%    eeg_event = read_brainvision_vmrk(fullfile(path_hdeeg,eeg_hdr.MarkerFile));
		eeg_event = bn_readExternalEventFile(fullfile(path_hdeeg,strcat(f,'.mrk')));
    
   
    [pcs_hdr, pcs_data] = bn_readActivaPC( fullfile(path_pcs,fname_pcs));
		[emg_hdr, emg_data] = bn_readEMG_wue( fullfile(path_emg,fname_emg));
    
    % pick the first channel
    pcs_data = pcs_data'.*1e3;
    pcs_data = pcs_data(:,pcs_cuttingtime*793:end);
    eeg_data = eeg_data(:,eeg_cuttingtime*eeg_hdr.Fs:end);

    pcs_ch_idx = 1;
    eeg_ch_idx = 10;
		emg_ch_idx = 3;

		pcs_locs = findTENSArtefact(pcs_data(pcs_ch_idx,:));
		eeg_locs = findTENSArtefact(eeg_data(eeg_ch_idx,:));
		emg_locs = findTENSArtefact(emg_data(emg_ch_idx,:));
    
%    pcs_ch = pcs_data(pcs_ch_idx,:);
%    pcs_ch_bp = process_bandpass('Compute', pcs_ch, 793, 90, 110, [],1);
%    pcs_ch_bp(:,:,2) = process_bandpass('Compute', pcs_ch, 793, 190, 210, [],1);
%    pcs_ch_bp(:,:,3) = process_bandpass('Compute', pcs_ch, 793, 290, 310, [],1);
%    pcs_ch_bp = process_bandpass('Compute', ...
%				mean(abs(pcs_ch_bp(:,:,:)),3), 793, .001, 1, [],1);
%
%    thr = 3*std(pcs_ch_bp);
%    [pks,locs] = findpeaks(pcs_ch_bp,'MINPEAKHEIGHT',thr,'MINPEAKDISTANCE',10*793);
%    if(length(locs)>2)
%        [val,idx] = sort(pks,'descend');
%        thr = val(3)+2*eps;
%    end
%    
%    [pcs_locs] = find(abs(diff((pcs_ch_bp>thr))));
%    if(length(pcs_locs)==3)
%        start_end = find(diff(pcs_locs)>5*793);
%        if (start_end == 1)
%            pcs_locs = [1 pcs_locs];
%        else
%            pcs_locs = pcs_locs(1:2);
%        end
%    end
%    
%    
%    eeg_ch = eeg_data(eeg_ch_idx,:);
%    eeg_ch_bp = process_bandpass('Compute', eeg_ch, ...
%				eeg_hdr.Fs, 90, 110, [],1);
%
%    eeg_ch_bp(:,:,2) = process_bandpass('Compute', eeg_ch,...
%				eeg_hdr.Fs , 190, 210, [],1);
%
%    eeg_ch_bp(:,:,3) = process_bandpass('Compute', eeg_ch, ....
%				eeg_hdr.Fs, 290, 310, [],1);
%
%    eeg_ch_bp = process_bandpass('Compute', ...
%				mean(abs(eeg_ch_bp(:,:,:)),3), eeg_hdr.Fs, 0.001, 1, [],1);
%
%    thr = 3*std(eeg_ch_bp);
%    [pks,locs] = findpeaks(eeg_ch_bp,'MINPEAKHEIGHT',thr,...
%				'MINPEAKDISTANCE',10*eeg_hdr.Fs);
%
%    if(length(locs)>2)
%        [val,idx] = sort(pks,'descend');
%        thr = val(3)+2*eps;
%    end
%    
%    [eeg_locs] = find(abs(diff((eeg_ch_bp>thr))));
%        if(length(eeg_locs)==3)
%        start_end = find(diff(eeg_locs)>5*eeg_hdr.Fs);
%        if (start_end == 1)
%            eeg_locs = [1 eeg_locs];
%        else
%            eeg_locs = eeg_locs(1:2);
%        end
%    end
%
%    emg_ch = emg_data(emg_ch_idx,:);
%    emg_ch_bp = process_bandpass('Compute', emg_ch, ...
%				emg_hdr.freq, 90, 110, [],1);
%
%    emg_ch_bp(:,:,2) = process_bandpass('Compute', emg_ch,...
%				emg_hdr.freq , 190, 210, [],1);
%
%    emg_ch_bp(:,:,3) = process_bandpass('Compute', emg_ch, ....
%				emg_hdr.freq, 290, 310, [],1);
%
%    emg_ch_bp = process_bandpass('Compute', ...
%				mean(abs(emg_ch_bp(:,:,:)),3), emg_hdr.freq, 0.001, 1, [],1);
%
%    thr = 3*std(emg_ch_bp);
%    [pks,locs] = findpeaks(emg_ch_bp,'MINPEAKHEIGHT',thr,...
%				'MINPEAKDISTANCE',10*emg_hdr.freq);
%
%    if(length(locs)>2)
%        [val,idx] = sort(pks,'descend');
%        thr = val(3)+2*eps;
%    end
%    
%    [emg_locs] = find(abs(diff((emg_ch_bp>thr))));
%        if(length(emg_locs)==3)
%        start_end = find(diff(emg_locs)>5*emg_hdr.freq);
%        if (start_end == 1)
%            emg_locs = [1 emg_locs];
%        else
%            emg_locs = emg_locs(1:2);
%        end
%    end

    method = min([length(eeg_locs)/2,length(pcs_locs)/2,length(emg_locs)/2])
    
    % actually compute t0 for all channels
    data_cell = [{eeg_ch},{pcs_ch},{emg_ch}];
    
    t0 = cellfun(@find_t_init,data_cell,{eeg_locs,pcs_locs,emg_locs},...
					{eeg_hdr.Fs,pcs_hdr.SenseChannelConfig.TDSampleRate,emg_hdr.freq},...
					{method, method,method},'uni',false);

    t0 = reshape([t0{:}],2,3)';
       
    if( method == 2 )
        % estimate the correct fs for PCS
        pcs_fs = (t0(2,2)-t0(2,1))* eeg_hdr.Fs /(t0(1,2)-t0(1,1))
    else
        t0(:,2) = [length(eeg_ch) length(pcs_ch) length(emg_ch)];
        pcs_fs = 793.65
    end
    
    
    % downsample pcs data to integer sampling frequency and eeg accordingly
    fs 			= 500;
    [pcs_data,~,~] = ResampleCascade(pcs_data,fs,pcs_fs);
    eeg_data 	= resample(eeg_data',fs,eeg_hdr.Fs)';
    emg_data 	= resample(emg_data',fs,emg_hdr.freq)';
    
		% each t0 row represent eeg,pcs,emg data before
		% we have to recompute the exact point in time after
		% resampling
    t0(1,:) = (round(t0(1,:)/eeg_hdr.Fs*pcs_fs));
    t0(3,:) = (round(t0(3,:)/emg_hdr.freq*pcs_fs));
    t0 			= round(t0/pcs_fs*fs);
           
		% also compute the sample indices for each events with new sampling freq
		for evIdx = 1:numel(eeg_event)
				eeg_event(evIdx).samples = round(eeg_event(evIdx).times * fs)...
						+ min(t0(:,1));

				eeg_event(evIdx).times	 = eeg_event(evIdx).times + min(t0(:,1))/fs;
		end


    % t0 has start and end samples for each channel
    % cut data to have the same number of samples
    onset  			= min(t0(:,1));
    off 		 		= [1 -1;1 -1; 1 -1] * onset ;
    data_wnd		= t0 - off;
        
    eeg_data_out 	= eeg_data(:,1+data_wnd(1,1):...
				min([data_wnd(1,2),length(eeg_data)]));

    pcs_data_out 	= pcs_data(:,1+data_wnd(2,1):...
				min([data_wnd(2,2),length(pcs_data)]));

    emg_data_out 	= emg_data(:,1+data_wnd(3,1):...
				min([data_wnd(3,2),length(emg_data)]));
    
    final_size   	= min([size(eeg_data_out,2),...
												size(pcs_data_out,2),....
												size(emg_data_out,2)]);
    data_out 	 		= [eeg_data_out(:,1:final_size); ...
										 pcs_data_out(:,1:final_size);...
										emg_data_out(:,1:final_size)];
    
    pcs_channels	= size(pcs_data,1);
		emg_channels  = size(emg_data,1);
    
    if(pcs_channels == 1)
        data_out  = [data_out; zeros(1,size(data_out,2))];
				pcs_channels = 2;
    end
    
    pcs_hdr.labels(end+1) 	= {'none'};
    eeg_hdr.chanunit([127 128]) = {'muV'};
		eeg_hdr.chanunit = [eeg_hdr.chanunit, emg_hdr.units];
    
    figure(1), clf
    subplot(311)
    plot(eeg_data(eeg_ch_idx,[-100:100] + t0(1))),
		hold on, plot(pcs_data(pcs_ch_idx,[-100:100] + t0(2)).*10,'r');
		plot(emg_data(emg_ch_idx,[-100:100] + t0(3)),'k');

    if(method == 2)
        subplot(312)
        plot(eeg_data(eeg_ch_idx,[-100:100] + t0(4))),
				hold on, plot(pcs_data(pcs_ch_idx,[-100:100] + t0(5)).*10,'r');
				plot(emg_data(emg_ch_idx,[-100:100] + t0(6)),'k');
    end
    subplot(313),
    plot(data_out([eeg_ch_idx end-1:end],:)');
    drawnow
    
    if(event_timeoffset)
        mtmdata_file=fullfile(event_offset_dir, [strjoin({fname_eeg_parts{2},...
						fname_eeg_parts{end-1}},'_'),'.txt']);

        [out, delimiter, headerlines] = importdata(mtmdata_file);
        offset = out.data(:,4);
        iMrk = find(strcmpi({eeg_event.label}, event_marker_label));
        if(length(offset)== size(eeg_event(iMrk).samples,2))
            for offset_idx = 1:length(offset)
            eeg_event(iMrk).samples(:,offset_idx)=...
								round(eeg_event(iMrk).samples(:,offset_idx)/...
                eeg_hdr.Fs*fs)+ repmat(round(offset(offset_idx)*1e-3*fs),...
                size(eeg_event(iMrk).samples,1),1)-data_wnd(1,1)...
								-eeg_cuttingtime*fs;
            end
        else 
            error('Offset length differs from event length');
        end
    end
    
    % update header info
    eeg_hdr.label 	= [eeg_hdr.label; pcs_hdr.labels';emg_hdr.labels'];
    eeg_hdr.nChans 	= eeg_hdr.NumberOfChannels + pcs_channels + emg_channels;
    eeg_hdr.NumberOfChannels = eeg_hdr.nChans;
    eeg_hdr.Fs			= fs;
    
    for ii = 1:2

        stn_pos_struct(ii)  = struct('type','stn',...
            'labels',pcs_hdr.labels(ii),...
            'sph_theta_besa',-134,...
            'sph_phi_besa',-45);
    end

		for ii = 1:emg_channels

        emg_pos_struct(ii)  = struct('type','emg',...
            'labels',emg_hdr.labels(ii),...
            'sph_theta_besa',-134,...
            'sph_phi_besa',-45);

		end
    
    eeg_hdr.layout.pos = [eeg_hdr.layout.pos, stn_pos_struct, emg_pos_struct];
		eeg_hdr.chantype(end:eeg_hdr.nChans) = {'other'};
    
    % this is the only supported data format
    eeg_hdr.DataFormat      = 'BINARY';
    eeg_hdr.DataOrientation = 'MULTIPLEXED';
    eeg_hdr.BinaryFormat    = 'IEEE_FLOAT_32';

		% no additional calibration needed, since float32
    eeg_hdr.resolution      = ones(size(eeg_hdr.label));      
    % write data
    filename = strjoin({fname_eeg_parts{1:2},[fname_eeg_parts{3},'-PCS'],...
				fname_eeg_parts{4:end-1}},'_');

    eeg_hdr.DataFile = strcat(filename,'.eeg');
    eeg_hdr.MarkerFile = strcat(filename,'.vmrk');
    
    write_brainvision_eeg(outdir, eeg_hdr, data_out);
    write_brainvision_vmrk(outdir, eeg_hdr, eeg_event);
    write_brainvision_vhdr(outdir, eeg_hdr);
end

end

function tau = find_t_init(D,locs,Fs,chunks)

% number of levels used for WICA and wavelet family
num_lvl   = 8;
vfilter   = 'db1';

% below we separate the portions containing the TENS artefacts

tau = [0 0];

for chunk = 0:chunks-1
    
    artefact_duration = locs(2+2*chunk)-locs(1+2*chunk);
    if((locs(2+2*chunk)+artefact_duration)<=length(D))
        data0 				= D(locs(1+2*chunk):(locs(2+2*chunk)+artefact_duration));
    else
        data0 				= D(locs(1+2*chunk):end);
    end
    data 				= data0;

    nsamples    = length(data);
    offset      = ceil(nsamples/(2^num_lvl)) * (2^num_lvl);
    
    data(nsamples+1:offset) = zeros(1,offset-nsamples);
    
    [thr, ~, ~] = ddencmp('den','wv',data);
    
    [swa, swd] = swt(data,num_lvl, vfilter);

    swd(abs(swd) < thr) = 0;
    swd(3:num_lvl,:) = zeros(size(swd(3:num_lvl,:)));
    swa = zeros(size(swa));
    
    out = iswt(swa,swd,vfilter)';
    out = out(1:nsamples);
    thresh = 2*std(abs(out((artefact_duration+1):end)));
    
    out( abs(out) <  std(out) ) = 0;
    
    [~, max_locs] = findpeaks((out),'MINPEAKHEIGHT',thresh);
    [~, min_locs] = findpeaks(-(out),'MINPEAKHEIGHT',thresh);  

    data_locs = max(max(max_locs),max(min_locs));
    
%     figure,
%     subplot(2,1,1), plot(t,out,t(data_locs),out(data_locs),'rx');
%     subplot(2,1,2), plot(t,data0,t(data_locs),data0(data_locs),'rx');
    
    tau(chunk+1) = data_locs + locs(1+2*chunk)-1;
end

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

    % We can't go below min(1,P/Q) without losing information. Because of low-pass filtering, don't be too precise
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
            len_resmp = ceil(size(x,2) * prod(Pfac) / prod(Qfac));
            nRow = size(x,1);
            x_resmp = zeros(nRow, len_resmp);
            % Loop on factors and rows
            for iRow = 1:size(x,1)
                x_tmp = x(iRow,:);
                for i = 1:iFacs
                    x_tmp = decimate(interp(x_tmp, Pfac(i)), Qfac(i));
                end
                x_resmp(iRow,:) = x_tmp;
            end
            x = x_resmp;
        % Resample takes vectorized inputs
        case 'resample'
            for i = 1:iFacs
                x = resample(x', Pfac(i), Qfac(i))';
            end
    end
end

function fname = filterFnames(fname,pattern)

		tmp = [{fname.name}];
		mask = ~cellfun(@isempty,regexp(tmp,pattern));
		fname = fname(mask);
end

function bool = checkDataConsistency(fname_hdeeg, fname_pcs)

		fname_hdeeg = [{fname_hdeeg.name}];
		fname_pcs = [{fname_pcs.name}];

		[p,hdeeg,e] = cellfun(@fileparts,fname_hdeeg,'uni',false);
		[p,pcs,e] 	= cellfun(@fileparts,fname_pcs,'uni',false);

		nEEGFiles = numel(hdeeg);
		nPCSFiles = numel(pcs);

		nMatchingFiles = sum(ismember(hdeeg,pcs));

		if nEEGFiles == nPCSFiles || nMatchingFiles == nEEGFiles
				bool = true;
		else
				bool = false;
		end
end


function locs = findTENSArtefact(data)
%FINDTENSARTEFACT data [1xN] time samples
%	LOCS = FINDTENSARTEFACT(DATA) Long description
%
%

		data_bp = process_bandpass('Compute', data, 793, 90, 110, [],1);
		data_bp(:,:,2) = process_bandpass('Compute', data, 793, 190, 210, [],1);
		data_bp(:,:,3) = process_bandpass('Compute', data, 793, 290, 310, [],1);
		data_bp = process_bandpass('Compute', ...
				mean(abs(data_bp(:,:,:)),3), 793, .001, 1, [],1);

		thr = 3*std(data_bp);
		[pks,locs] = findpeaks(data_bp,'MINPEAKHEIGHT',thr,'MINPEAKDISTANCE',10*793);
		if(length(locs)>2)
				[val,idx] = sort(pks,'descend');
				thr = val(3)+2*eps;
		end

		[locs] = find(abs(diff((data_bp>thr))));
		if(length(locs)==3)
				start_end = find(diff(locs)>5*793);
				if (start_end == 1)
						locs = [1 locs];
				else
						locs = locs(1:2);
				end
		end



end


