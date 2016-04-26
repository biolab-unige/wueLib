function tau = manualTENS(D,method)
%MANUALTENS Description
%	TAU = MANUALTENS(D,locs) Long description
%
	warnMessage = sprintf('You should pick only %d point(s)',method);
	warndlg(warnMessage);
	f1 = figure;
	imf = ceemdan(D,0.0002,20,100,1);
	plot(bsxfun(@plus,imf,max(imf(:))*(1:(size(imf,1)))')')
	addCrossair(f1);

	waitfor(f1)

	tau = evalin('base','cursorValue');
	size(tau)
	

end
