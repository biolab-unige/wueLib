function bool = wlb_checkDataConsistency(fnameMod1, fnameMod2)
%CHECKDATACONSISTENCY Description
%	BOOL = CHECKDATACONSISTENCY(FNAMEMOD1, FNAMEMOD2) Long description
%

fnameMod1 = {fnameMod1.name};
fnameMod2 = {fnameMod2.name};

[~,mod1,~] = cellfun(@fileparts,fnameMod1,'uni',false);
[~,mod2,~] = cellfun(@fileparts,fnameMod2,'uni',false);

nMod1Files = numel(mod1);
nMod2Files = numel(mod2);

nMatchingFiles = sum(ismember(mod1,mod2));

if nMod1Files == nMod2Files || nMatchingFiles == nMod1Files
    bool = true;
else
    bool = false;
end
end
