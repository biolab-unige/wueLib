function write_brainvision_vmrk(path,hdr,event)
%WRITE_BRAINVISION_VMRK Description
%	 = WRITE_BRAINVISION_VMRK(PATH,HDR,EVENT) Long description
%

markerfile = fullfile(path, hdr.MarkerFile);
datafile   = hdr.DataFile;

% open the header file and write the ascii header information
fid = fopen(markerfile, 'w');

fprintf(fid, 'Brain Vision Data Exchange Marker File, Version 1.0\r\n');
fprintf(fid, '\r\n');
fprintf(fid, '[Common Infos]\r\n');
fprintf(fid, 'Codepage=UTF-8\r\n');
fprintf(fid, 'DataFile=%s\r\n', datafile);
fprintf(fid, '\r\n');
fprintf(fid, '[Marker Infos]\r\n');
fprintf(fid, '; Each entry: Mk<Marker number>=<Type>,<Description>,<Position in data points>,\r\n');
fprintf(fid, '; <Size in data points>, <Channel number (0 = marker is related to all channels)>\r\n');
fprintf(fid, '; Fields are delimited by commas, some fields might be omitted (empty).\r\n');
fprintf(fid, '; Commas in type or description text are coded as "\\1".\r\n');

idx = 1;
for iEvt=1:length(event)
    for iEpc=1:length(event(iEvt).epochs)        
        fprintf(fid, 'Mk%d=%s,%s,%lu,%lu,%lu\r\n',...
						idx,event(iEvt).type{iEpc}, event(iEvt).label, ...
						event(iEvt).samples(iEpc),event(iEvt).length(iEpc),...
						event(iEvt).chan_num(iEpc));

        idx = idx +1;
    end
end

fclose(fid);
