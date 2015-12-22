function [event] = read_brainvision_vmrk(filename)

% READ_BRAINVISION_VMRK reads the markers and latencies
% it returns the stimulus/response code and latency in ms.
%
% Use as
%   [stim, resp, segment, timezero] = read_brainvision_vmrk(filename)
%
% This function needs to read the header from a separate file and
% assumes that it is located at the same location.
%
% See also READ_BRAINVISION_VHDR, READ_BRAINVISION_EEG

% original M. Schulte 31.07.2003
% modifications R. Oostenveld 14.08.2003
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: read_brainvision_vmrk.m 7123 2012-12-06 21:21:38Z roboos $


[p, f, x] = fileparts(filename);
filename = fullfile(p,[f '.vmrk']);


fid=fopen(filename,'rt');
if fid==-1,
    error('cannot open marker file')
end

while ~feof(fid)
    line=fgetl(fid);
    % pause
    if ~isempty(sscanf(line,'Mk%d'))
        A = tokenize(line,',');
        index = sscanf(A{1},'Mk%d=');
        B = tokenize(A{1},'=');
        epochs = 1;
        event(index).epochs = epochs;
        event(index).type{epochs} = B{2};
        event(index).label = A{2};
        event(index).samples(epochs) = str2double(A{3});
        event(index).length(epochs) = str2double(A{4});
        event(index).chan_num(epochs) = str2double(A{5});
    end
end

fclose(fid);
