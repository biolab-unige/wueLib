function write_brainvision_vhdr(path, hdr)

datafile   = hdr.DataFile;
markerfile = hdr.MarkerFile;
[p, f, x] = fileparts(datafile);
headerfile = fullfile(path, [f '.vhdr']);

% open the header file and write the ascii header information
fid = fopen(headerfile, 'w');
fprintf(fid, 'Brain Vision Data Exchange Header File Version 1.0\r\n');
fprintf(fid, '; Data created by FieldTrip\r\n');
fprintf(fid, '\r\n');
fprintf(fid, '[Common Infos]\r\n');

fprintf(fid, 'DataFile=%s\r\n', datafile);
if ~isempty(markerfile)
  fprintf(fid, 'MarkerFile=%s\r\n',      markerfile);
end
fprintf(fid, 'DataFormat=%s\r\n',        hdr.DataFormat);
fprintf(fid, 'DataOrientation=%s\r\n',   hdr.DataOrientation);
fprintf(fid, 'NumberOfChannels=%d\r\n',  hdr.NumberOfChannels);
% Sampling interval in microseconds
fprintf(fid, 'SamplingInterval=%d\r\n',  round(1e6/hdr.Fs));
fprintf(fid, '\r\n');
fprintf(fid, '[Binary Infos]\r\n');
fprintf(fid, 'BinaryFormat=%s\r\n',      hdr.BinaryFormat);
fprintf(fid, '\r\n');
fprintf(fid, '[Channel Infos]\r\n');
% Each entry: Ch<Channel number>=<Name>,<Reference channel name>,<Resolution in microvolts>,<Future extensions>...
% Fields are delimited by commas, some fields might be omitted (empty).
% Commas in channel names should be coded as "\1", but are not supported here
for i=1:hdr.NumberOfChannels
  fprintf(fid, 'Ch%d=%s,,%g,%s,%s\r\n', i, hdr.label{i}, hdr.resolution(i),hdr.chanunit{i},hdr.chantype{i});
end

fprintf(fid, '\r\n');
fprintf(fid, '[Coordinates]\r\n');
fprintf(fid, '; Electrode Position File: /dev/null \r\n');

for i=1:hdr.NumberOfChannels
  fprintf(fid, 'Ch%d=%d,%d,%d\r\n', i, 1, hdr.layout.pos(i).sph_theta_besa, hdr.layout.pos(i).sph_phi_besa);
end

fclose(fid);