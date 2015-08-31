function APM2Brainamp(Filename)

apmdata = APMReadData(Filename);

hdr.DataFormat      = 'BINARY';
hdr.DataOrientation = 'MULTIPLEXED';
hdr.BinaryFormat    = 'IEEE_FLOAT_32';

[p,f,x] = fileparts(Filename);
% write data
hdr.DataFile = strcat(f,'.eeg');
hdr.MarkerFile = [];

split = strsplit(f,'_');
channel = split{7};
channel = channel(2:end);

hdr.nChans 	= 0.5*(length(channel));
for i = 0:hdr.nChans-1
    hdr.label{i+1} 	= channel((1:2)+i*2);
end

hdr.Fs	    = apmdata.channels(1).sampling_frequency;
hdr.resolution = ones(hdr.nChans,1);
hdr.chantype = repmat({'eeg'}, size(hdr.label));
hdr.chanunit = repmat({'uV'},  size(hdr.label));

maxlength = min([apmdata.channels(1:hdr.nChans).continuous_samples]);
volt_calib = [apmdata.channels.voltage_calibration]';
for i = 1:hdr.nChans
    data_out(i,:) = apmdata.channels(i).continuous(1:maxlength);
end

data_out =  bsxfun(@times,data_out,volt_calib(1:hdr.nChans));

write_brainvision_eeg(p, hdr, data_out);
write_brainvision_vhdr(p, hdr);

