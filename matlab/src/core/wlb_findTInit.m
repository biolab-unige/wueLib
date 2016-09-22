function tau = wlb_findTInit(D,locs,freq,chunks)
%FINDTINIT Description
%	TAU = FINDTINIT(D,LOCS,CHUNKS) Long description
%

% number of levels used for WICA and wavelet family
numLvl    = 10;
vfilter   = 'db1';
tau       = nan(2,1);

for chunk = 0:chunks-1
    
    artefactIdx	= round(locs(1+(2*chunk))):round(locs(2+(2*chunk))) ;
    
    artefactIdx(artefactIdx > numel(D)) = [];
    artefactIdx(artefactIdx < 1) = [];
    data 				= D(artefactIdx);
    
    nsamples    = length(data);
    offset      = ceil(nsamples/(2^numLvl)) * (2^numLvl);
    
    data(nsamples+1:offset) = zeros(1,offset-nsamples);
    [thr] = 0.5*ddencmp('den','wv',data);
                
    [swa, swd] = swt(data,numLvl, vfilter);
    
    swd(abs(swd) < thr) = 0;
    swd([3:numLvl],:) = zeros(size(swd([3:numLvl],:)));
    swa = zeros(size(swa));
    
    out = iswt(swa,swd,vfilter)';
    out = out(1:nsamples);

%     imf=emd(data);   
%     out = sum(imf(1:2,:),1);

    thresh = std((out));
       
    out( abs(out) <  thresh ) = 0;
    
    [~, maxLocs] = findpeaks(((out)),'MINPEAKHEIGHT',thresh);
    datamaxLocs = max(maxLocs);
    [~, minLocs] = findpeaks((-(out)),'MINPEAKHEIGHT',thresh);
    dataminLocs = max(minLocs);
    
    tau(chunk+1) = max(datamaxLocs,dataminLocs) + artefactIdx(1)-1;
    
end
tau(isnan(tau)) = [];
figure,
subplot(211)
plot(([-100:100])/freq*1000,D(tau(1) + [-100:100])) , hold on, plot(0,D(tau(1)),'rx');
xlim(([-100 100])/freq*1000)
if(length(tau)==2)
subplot(212)
plot(([-100:100])/freq*1000,D(tau(2) + [-100:100])) , hold on, plot(0,D(tau(2)),'rx');
xlim(([-100 100])/freq*1000)
end
title('t0 points');

end