function t0 = wlb_readTensTiming(filename,trial)
%WLB_READTENSTIMING Description
%	T0 = WLB_READTENSTIMING(FILENAME) Long description
%

[pcsFirstTens,emgFirstTens,pcsLastTens,emgLastTens] = textread(filename,'%d %d %d %d\n','delimiter',',');

t0 = [pcsFirstTens(trial), pcsLastTens(trial);emgFirstTens(trial),emgLastTens(trial)];

end
