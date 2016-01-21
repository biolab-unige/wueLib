function write_brainvision_eeg(path,hdr, dat)

% WRITE_BRAINVISION_EEG exports continuous EEG data to a BrainVision *.eeg
% and corresponding *.vhdr file. The samples in the exported file are
% multiplexed and stored in ieee-le float32 format.
%
% Use as
%   write_brainvision_eeg(filename, hdr, dat)
%
% See also READ_BRAINVISION_EEG, READ_BRAINVISION_VHDR, READ_BRAINVISION_VMRK

% Copyright (C) 2007, Robert Oostenveld
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
% $Id$

if length(size(dat))>2
  ntrl  = size(dat,1);
  nchan = size(dat,2);
  nsmp  = size(dat,3);
else
  nchan = size(dat,1);
  nsmp  = size(dat,2);
end

if hdr.NumberOfChannels~=nchan
  error('number of channels in in header does not match with the data');
end

% determine the filenames
datafile   = fullfile(path, hdr.DataFile);

% open the data file and write the binary data
fid = fopen(datafile, 'wb', 'ieee-le');
if length(size(dat))>2
  warning('writing segmented data as if it were continuous');
  for i=1:ntrl
    fwrite(fid, squeeze(dat(i,:,:)), 'float32');
  end
else
  fwrite(fid, dat, 'float32');
end

fclose(fid);


