function [hdr, data]= wlb_readEMG_wue(filename)
% wlb_readEMG_wue

% Edited 2014-09-19 by Gabriele Arnulfo <gabriele.arnulfo@gmail.com>

[p,f,~] = fileparts(filename);

fid = fopen(fullfile(p,[f,'.txt']),'r');

hdr.person 	= fgetl(fid);
hdr.sex			= fgetl(fid);
hdr.born		= fgetl(fid);
hdr.code    = fgetl(fid);
fgetl(fid);
hdr.record  = fgetl(fid);
% jump  Application Generic
fgetl(fid);
hdr.creation= fgetl(fid); %  03/07/2014 11:53:07
% jump Exercises   {
fgetl(fid);
% jump #   Name    Start,sec   Length,sec  Start time
fgetl(fid);
hdr.exercises = fgetl(fid); %1   Combined    0.000   23.800  11:52:34.38
%jump closing bracket
fgetl(fid);
fgetl(fid);
hdr.freq = sscanf(fgetl(fid),'%*s%d');
hdr.samples = sscanf(fgetl(fid),'%*s%d');
fgetl(fid);
stringS = fgetl(fid);
labels  = regexp(regexprep(stringS,'\s','_'),',m[V|s]_','split');
units	  = regexp(stringS,'m[V|s]','match');


hdr.n_chan = numel(labels);
hdr.units = units;
hdr.labels = labels;

% jump two lines to remove Calibr field which seems to be 
% useless
fgetl(fid);
fgetl(fid);

data = nan(hdr.samples,hdr.n_chan+1);

for lines=1:hdr.samples
		data(lines,:) = sscanf(fgetl(fid),'%f')';
end

data = data(:,2:end)';
fclose(fid);
end
