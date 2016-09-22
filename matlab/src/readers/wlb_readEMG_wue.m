function [hdr, data,eventStruct]= wlb_readEMG_wue(filename)

eventStruct = [];

[p,f,e] = fileparts(filename);

if(strcmp(e,'.tdf'))
    [startTime,frequency,emgMap,labels,data] = tdfReadDataEmg (filename);
    hdr.nSamples = size(data,2);   
    hdr.freq = frequency;
    for i=1:numel(labels)
       C = strsplit(labels{i},' ') ;
       if((isempty(C{end})))
           C{end} = '2';
       end
       labels{i} = cellfun(@(x)x(1),C);
    end
    
    no_emg_idx = find(cellfun(@(x)not(isempty(x)),strfind(labels,'ES')));
    [pxx,f] = pmtm(data(no_emg_idx,:)',3,[100:2:hdr.freq/2],hdr.freq);
    [dum,idx] = sort(sum(pxx));
    
    pulse_idx = no_emg_idx(idx(2));
    labels(pulse_idx) = strcat(labels(pulse_idx),'_pulse');
    hdr.chantype(pulse_idx) = {'pulse'};
    ecg_idx = no_emg_idx(idx(1));
    labels(ecg_idx) = strcat(labels(ecg_idx),'_ecg');
    hdr.chantype(ecg_idx) = {'ecg'};
    emg_idx = setdiff([1:numel(labels)],[pulse_idx ecg_idx]);
    labels(emg_idx) = strcat(labels(emg_idx),'_emg');
    hdr.labels = labels;
    hdr.chantype(emg_idx) = {'emg'};
    
    [startTime,frequency,gpMap,labels,gpData] = tdfReadDataGenPurpose (filename);

    IMU_channel = find(cellfun(@(x)not(isempty(x)),strfind(num2cell(labels,2),'IMU')));
    if(not(isempty(IMU_channel)))
        IMU_data = gpData(IMU_channel,:);
        [pks,rise_locs] = findpeaks(diff(IMU_data),'Threshold',3);
        [pks,fall_locs] = findpeaks(-diff(IMU_data),'Threshold',3);
        locs = rise_locs + 0*round(fall_locs - rise_locs)/2;
        if(not(isempty(locs)))
            locs = locs +1;
            nEvents = length(locs);
            eventStruct.type	= repmat({'custom'},[1 nEvents]);
            eventStruct.label 	= 'IMU';
            eventStruct.samples = round(locs/frequency*hdr.freq);
            eventStruct.times = (eventStruct.samples/hdr.freq);
            eventStruct.epochs 	= ones(1,nEvents);
            eventStruct.length 	= ones(1,nEvents);
            eventStruct.chan_num = zeros(1,nEvents);
        end
    end
elseif(strcmp(e,'.txt'))
    [hdr, data]= txtReadDataEmg(filename);
    for i=1:numel(hdr.labels)
       C = strsplit(hdr.labels{i},'_') ;
       labels{i} = cellfun(@(x)x(1),C(1:end-1));
       hdr.chantype{i} = C{end};
    end
    hdr.nSamples = size(data,2);
else
    error('file emg not supported');
end