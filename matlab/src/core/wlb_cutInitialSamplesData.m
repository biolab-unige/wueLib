function data = cutInitialSamplesData(data,offset,fs)
%CUTINITIALSAMPLESDATA Description
%	DATA = CUTINITIALSAMPLESDATA(DATA,OFFSET,FS) Long description
%
	data = data(:,(offset(1)*fs)+1:end-(offset(2)*fs));

end
