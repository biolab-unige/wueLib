function newloc = wlb_findTENSArtefact(data, fs,pctg,modality)
%FINDTENSARTEFACT data [1xN] time samples
%	LOCS = FINDTENSARTEFACT(DATA) Long description
global globDebug
if globDebug
    figure, subplot(2,2,[1 2]),plot([0:length(data)-1]/fs,data/max(data)), hold on;
    title([modality ' findTENSArtefact']);
end
% data = abs(data);

[pxx,f] = pmtm(data,3,[60:2:fs/2],fs);
dev_std = std(pxx);
[pks,locs] = findpeaks(pxx);
[val,idx] = sort(pks,'descend');

TENSfreq = f(locs(idx));

for i=1:3
    dataBp(:,:,i) = wlb_bandpass_fft(data, fs,TENSfreq(i)-10, min(TENSfreq(i)+10,fs/2-20), 1,1,[]);
end
dataBp = wlb_bandpass_fft(mean(abs(dataBp(:,:,:)),3), fs, 0.001, 1, 1,1,[]);

thr = prctile(dataBp,pctg);
while(thr<0)
    pctg = pctg +1;
    thr = prctile(dataBp,pctg);
end

[pks,pksLocs] = findpeaks(dataBp,'MINPEAKHEIGHT',thr,'MINPEAKDISTANCE',10*fs);
if(length(pksLocs)>2)
    % in general we should have only two tens
    % but we might end up with some weird artefacts
       
    [val,I] = sort(pks,'descend');
    pks 		= val(1:2);
    pksLocs	= pksLocs(I);
    pksLocs	= pksLocs(1:2);
    
end

[locs] = find(abs(diff((dataBp>thr))));
dataBpDer = -savitzkyGolayFilt(dataBp,4,1,11);
locsDer = dataBpDer(locs);

if(locsDer(1)<0),locs(1)=[];end
if(locsDer(end)>0),locs(end)=[];end

if globDebug
    subplot(2,2,[1 2]),plot([0:length(dataBp)-1]/fs,dataBp/max(dataBp),'r');
    plot((locs-1)/fs,dataBp(locs)/max(dataBp),'ko');
end

locs = reshape(locs,2,numel(locs)/2);

duration = diff(locs);

locs = locs(:,duration>0.8*fs);

newloc = [];

%  		if(length(locs)/numel(pks)~=2) %ci sono piu di 2 intersezioni in uno o entrambi i picchi
for i=1:numel(pks)
    loc = locs( abs(locs-pksLocs(i))<15*fs ); %cerco le intersezioni piu vicine al picco
    startLoc = find(loc-pksLocs(i)<0,1,'first');
    endLoc = find(loc-pksLocs(i)>0,1,'last');
    loc = [loc(startLoc) loc(endLoc)];
    if ~isempty( loc )
        newloc = [newloc loc];
    end
end
%  		else
%  				newloc = locs;
%  		end
%
newloc = sort(newloc,'ascend');

if globDebug
    subplot(2,2,[1 2]),plot((newloc-1)/fs,dataBp(newloc)/max(dataBp),'kx','MarkerSize',5);
    subplot(2,2,3),plot(((newloc(1):newloc(2))-1)/fs,data(newloc(1):newloc(2))/max(data));
    xlim(([newloc(1) newloc(2)]-1)/fs)
    if(length(newloc)>2)
        subplot(2,2,4),plot(((newloc(3):newloc(4))-1)/fs,data(newloc(3):newloc(4))/max(data));
        xlim(([newloc(3) newloc(4)]-1)/fs)

    end
end


end
