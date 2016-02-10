function newloc = findTENSArtefact(data, fs,pctg)
%FINDTENSARTEFACT data [1xN] time samples
%	LOCS = FINDTENSARTEFACT(DATA) Long description
	  global globDebug	
		if globDebug
				figure, plot(data.*1e-2), hold on;
                title('findTENSArtefact');
		end
		data = abs(data);

		dataBp = wlb_bandpass_fft(data, fs, 90, 110,1,1,[]);
		dataBp(:,:,2) = wlb_bandpass_fft(data, fs, 190, 200, 1,1,[]);
        
		if fs > 700			
            
			dataBp(:,:,3) = wlb_bandpass_fft(data, fs, 290, 310, 1,1,[]);
			
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
				plot(dataBp,'r');
				plot(locs,dataBp(locs),'ko');
		end

	 
		locs = reshape(locs,2,numel(locs)/2);
		
		duration = diff(locs);
		
		locs = locs(:,duration>0.8*fs);
        
		newloc = [];
		
%  		if(length(locs)/numel(pks)~=2) %ci sono piu di 2 intersezioni in uno o entrambi i picchi
				for i=1:numel(pks)
						loc = locs( abs(locs-pksLocs(i))<5*fs ); %cerco le intersezioni piu vicine al picco
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
				plot(newloc,dataBp(newloc),'kx','MarkerSize',5);
		end


end

