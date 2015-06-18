function [hdr data]= bn_readEMG_wue(filename)
% bn_readEMG_wue

% Edited 2014-09-19 by Gabriele Arnulfo <gabriele.arnulfo@gmail.com>

fid = fopen(filename,'r');

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
labels = textscan(fgetl(fid),'%s','delimiter','\t');
labels = labels{1};

labels = regexp(labels,',','split');

hdr.n_chan = numel(labels);

for l = 1:hdr.n_chan
		hdr.labels(l) = regexprep(labels{l}(1),' ','_');
		if( numel(labels{l}) == 2)
			hdr.units(l)  = labels{l}(2);
		else
			hdr.units(l)  = {'unk'};
		end
end
fgetl(fid);
fgetl(fid);

data = nan(hdr.samples,hdr.n_chan);

for lines=1:hdr.samples
		data(lines,:) = sscanf(fgetl(fid),'%f')';
end

data = data';
fclose(fid);
end
