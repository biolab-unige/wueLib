function [out] = wlb_cutInitialSamplesData(data,offset,fs)
%CUTINITIALSAMPLESDATA Description
%	DATA = CUTINITIALSAMPLESDATA(DATA,OFFSET,FS) Long description
%
out= zeros(size(data));
out(:,round(offset(1)*fs)+1:end-round(offset(2)*fs)) = data(:,round(offset(1)*fs)+1:end-round(offset(2)*fs));
