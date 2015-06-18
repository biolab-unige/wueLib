function [hdr, data] = bn_readActivaPC(filename)
%BN_READACTIVAPC read PC+s xml/txt files 
%	[HDR, DATA] = BN_READACTIVAPC(FILENAME)
%	data [ch,time] represents pcs data in mV 
%

% Edited 2014-09-22 by Gabriele Arnulfo <gabriele.arnulfo@gmail.com>
import javax.xml.xpath.*
factory = XPathFactory.newInstance;
xpath = factory.newXPath;

switch(class(filename))
	case 'cell'
		filename = filename{:};
	otherwise
		filename = filename(:)';
end

[p,f,e] = fileparts(filename); 

xDoc = xmlread(fullfile(p,strcat(f,'.xml')));
docNode = xDoc.getDocumentElement;

Attributes = [{'PatientID'} {'RecordingDuration'} ...
			{'INSTimeStamp'} {'SPTimeStamp'} ...
			{'SenseChannelConfig/TDSampleRate'} ... 
			{'SenseChannelConfig/Channel1/ChannelType'} ...
			{'SenseChannelConfig/Channel1/PlusInput'} ...
			{'SenseChannelConfig/Channel1/MinusInput'} ... 
			{'SenseChannelConfig/Channel2/ChannelType'} ...
			{'SenseChannelConfig/Channel2/PlusInput'} ...
			{'SenseChannelConfig/Channel2/MinusInput'} ...
			{'SenseChannelConfig/Channel3/ChannelType'} ...
			{'SenseChannelConfig/Channel3/PlusInput'} ...
			{'SenseChannelConfig/Channel3/MinusInput'} ...
			{'SenseChannelConfig/Channel4/ChannelType'} ...
			{'SenseChannelConfig/Channel4/PlusInput'} ....
			{'SenseChannelConfig/Channel4/MinusInput'}];

for ele = 1:numel(Attributes)
	expression = xpath.compile(Attributes{ele});
	nestedNames = regexp(Attributes{ele},'/','split');
	str2eval = '';
	for n = 1:numel(nestedNames)
		str2eval = strcat(str2eval,'.(nestedNames{',num2str(n),'})');
	end
	
	node = expression.evaluate(docNode,XPathConstants.NODE);
	eval(strcat('hdr',str2eval,' = 	node.item(0).getData;'));
end

% we now should extract the labels from current hdr
chan_num = 1;
for ii = 1:4
	chan = strcat('Channel',num2str(ii));
	if( strcmp(hdr.SenseChannelConfig.(chan).ChannelType,'TD'))
		ch_index = str2num(cell2mat(regexp(hdr.SenseChannelConfig.(chan).PlusInput.toCharArray','\d+','match')));
		if( ch_index > 6)
			hemisphere(chan_num) = {'l'};
		else
				hemisphere(chan_num) = {'r'};
		end
		labels(chan_num) = {strcat(hdr.SenseChannelConfig.(chan).PlusInput.toCharArray','+',...
				hdr.SenseChannelConfig.(chan).MinusInput.toCharArray','-')};
		chan_num = chan_num + 1;
	end
end

hdr.labels = labels;
hdr.recordingChannels = chan_num;
hdr.SenseChannelConfig.TDSampleRate = cellfun(@str2num,(regexp(hdr.SenseChannelConfig.TDSampleRate.toCharArray','\d+','match')));
hdr.hemi_flag = hemisphere;

% read
data = importdata(fullfile(p,filesep,strcat(f,'.txt')));
data = data(:,sum(data,1) > 0); 

data = data';

end
