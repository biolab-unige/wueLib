function tau = findTInit(D,locs,~,chunks)
%FINDTINIT Description
%	TAU = FINDTINIT(D,LOCS,CHUNKS) Long description
%

		% number of levels used for WICA and wavelet family
		numLvl    = 8;
		vfilter   = 'db1';
		tau 			= nan(2,1);


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
				swd(3:numLvl,:) = zeros(size(swd(3:numLvl,:)));
				swa = zeros(size(swa));
				
				out = iswt(swa,swd,vfilter)';
				out = out(1:nsamples);
				thresh = 2*std(abs(out));
				
				out( abs(out) <  std(abs(out)) ) = 0;
				
				[~, maxLocs] = findpeaks(abs((out)),'MINPEAKHEIGHT',thresh);
				dataLocs = max(maxLocs);				
			
				tau(chunk+1) = dataLocs + artefactIdx(1)-1;

		end
		tau(isnan(tau)) = [];
		figure,
		plot(D) , hold on, plot(tau,D(tau),'rx');
        title('t0 points');

end

