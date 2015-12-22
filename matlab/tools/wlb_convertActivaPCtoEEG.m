function wlb_convertActivaPCtoEEG(folder,logFid)
%WLB_CONVERTACTIVAPCTOEEG Append to a single BrainAmp file all the PC+s data
%	WLB_CONVERTACTIVAPCTOEEG(FOLDER) 
%


		fprintf(logFid,'Reading in %s\n',folder);
		fprintf('Reading in %s\n',folder);
		fnames = dir(fullfile(folder,'*xml'));
		
		dataOut = cell(numel(fnames),1);
		labels  = cell(numel(fnames),1);

		for file = 1:numel(fnames)
				try
					[hdr, data]	= wlb_readActivaPC(fullfile(folder,fnames(file).name));
				catch ME
						fprintf(logFid,ME.message);
						continue
				end

				durSeconds(file)  	= sum(str2double(regexp(char(...
						hdr.RecordingDuration),':','split')).*[3600,60,1]);
				fs(file)						= hdr.SenseChannelConfig.TDSampleRate;
				labels(file)				= {hdr.labels};
				dataOut(file)				= {data'};
                figure(1),clf,plot(data')
		end

		realFs = unique(fs);
		if numel(realFs) > 1
				error(' Inconsistent sampling frequency between selected files\n');
		end

		clear file hdr data;

		data = [dataOut{:}]';

		out_hdr.label = [labels{:}];
		out_hdr.Fs 	  = realFs;


		out_hdr.nChans 	= numel(out_hdr.label);
		out_hdr.NumberOfChannels = out_hdr.nChans;
		out_hdr.Fs			= realFs;
		out_hdr.chanunit(out_hdr.nChans) = {'mV'};

		stn_pos_struct = struct('type','stn',...
						'labels',out_hdr.label,...
						'sph_theta_besa',-134,...
						'sph_phi_besa',-45);

		
		out_hdr.layout.pos = stn_pos_struct;
		out_hdr.chantype(out_hdr.nChans) = {'other'};
		
		% this is the only supported data format
		out_hdr.DataFormat      = 'BINARY';
		out_hdr.DataOrientation = 'MULTIPLEXED';
		out_hdr.BinaryFormat    = 'IEEE_FLOAT_32';

		% no additional calibration needed, since float32
		out_hdr.resolution      = ones(size(out_hdr.label));      

		% write data
		stripPaths = regexp(folder,'/orig/','split');

		filename = regexprep(stripPaths{2},'/','_');


		
		out_hdr.DataFile = strcat(filename,'.eeg');
		out_hdr.MarkerFile = '';
		fprintf(logFid,'Writing to %s\n',filename);
						
		write_brainvision_eeg(folder, out_hdr, data);
		write_brainvision_vhdr(folder, out_hdr);


end
