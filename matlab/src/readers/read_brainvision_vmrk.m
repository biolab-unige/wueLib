function events = read_brainvision_vmrk(EventFile)
% IN_EVENTS_BRAINAMP: Open a BrainVision BrainAmp .vmrk file.
%
% OUTPUT:
%    - events(i): array of structures with following fields (one structure per event type)
%        |- label   : Identifier of event #i
%        |- samples : Array of unique time indices for event #i in the corresponding raw file
%        |- times   : Array of unique time latencies (in seconds) for event #i in the corresponding raw file
%                     => Not defined for files read from -eve.fif files

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
%
% Copyright (c)2000-2015 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, 2012

% Open and read file
fid = fopen(EventFile,'r');
% Markers list
Markers = {};
uniqueEvt = [];
isMarkerSection = 0;
% Read file line by line
if fid~=-1
    while 1
        % Read one line
        newLine = fgetl(fid);
        if ~ischar(newLine)
            break;
        end
        % Lines to skip
        if isempty(newLine)
            continue;
        elseif ~isempty(strfind(newLine, '[Marker Infos]'))
            isMarkerSection = 1;
        elseif ~isMarkerSection || ismember(newLine(1), {'[', ';', char(10), char(13)}) || ~any(newLine == '=')
            continue;
        end
        % Split around the '=' and ','
				argLine = strtrim(regexp(newLine,'[=|,]','split'));
%        argLine = strtrim(str_split(newLine, '=,', 0));
        if (length(argLine) < 6) || (length(argLine{1}) < 2)
            continue;
        end
        % Markers start with 'Mk'
        if ~strcmpi(argLine{1}(1:2), 'Mk')
            continue;
        end
        % Marker label
        if ~isempty(argLine{3})
            mlabel = argLine{3};
        else
            mlabel = 'Mk';
        end
        % Add markers entry: {name, type, start, length, channel num}
        Markers(end+1,:) = {mlabel, argLine{2}, str2num(argLine{4}), str2num(argLine{5}), str2num(argLine{6})};
    end
    fclose(fid);
    
    uniqueEvt = unique(Markers(:,1)');
end
% Close file

% List of events
% Initialize returned structure
events = repmat(struct('label',{},'type',{},'epochs',{},'samples',{},...
    'length',{},'chan_num',{}), [1, length(uniqueEvt)]);

% Create events list
for iEvt = 1:length(uniqueEvt)
    % Find all the occurrences of event #iEvt
    iMrk = find(strcmpi(Markers(:,1)', uniqueEvt{iEvt}));
    % Add event structure
    events(iEvt).label   = uniqueEvt{iEvt};
    events(iEvt).type   = Markers(iMrk,2);
    events(iEvt).epochs  = ones(1, length(iMrk));
    events(iEvt).samples = [Markers{iMrk,3}];
    events(iEvt).length = [Markers{iMrk,4}];
    events(iEvt).chan_num = [Markers{iMrk,5}];
end

end
