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
stringS = textscan(fgetl(fid),'%s','delimiter','\t');
units	  = regexp(stringS{1},'m[V|s]','match');
units(cellfun(@isempty,units)) = {'unk'};
units = [units{:}];
labels = stringS{1};

miscPattern = [{'time'},{'artef.*'},{'digital'},{'tens'}];
mask = cellfun(@(x)(regexp(x,miscPattern)),lower(labels),'Unif',false);
mask = [mask{:}];
mask(cellfun(@isempty,mask)) = {0};

% labels x patterns
mask = reshape(mask,[numel(miscPattern),numel(labels)])';
mask = cell2mat(mask);
mask = logical(sum(mask,2));

labels(mask) = strcat(labels(mask),'_pulse');
labels(~mask) = strcat(labels(~mask),'_emg');

mask = strfind(lower(labels),'ekg');
mask(cellfun(@isempty,mask)) = {0};
mask = logical([mask{:}]);
labels(mask) = {'EKG'};
labels = regexprep(labels,',m[V|s]','');
labels = regexprep(labels,'\s+','_');


% jump two lines to remove Calibr field which seems to be 
% useless
fgetl(fid);
fgetl(fid);

data = nan(hdr.samples,numel(labels));

for lines=1:hdr.samples
		data(lines,:) = sscanf(fgetl(fid),'%f')';
end

data = data';
fclose(fid);

channelZeros = find(mean(data,2)==0);
data( channelZeros,:) = [];

% for some weird reasons we mihgt end up with an empyt
% labels due to random white spaces in the line
labels(channelZeros) = [];
units(channelZeros) = [];

hdr.n_chan = numel(labels);
hdr.units = units;
hdr.labels = labels';

end
