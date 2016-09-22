function fnames = wlb_filterFnames(fnames,pattern)
%FILTERFNAMES Description
%	FNAME = FILTERFNAMES(FNAMES,PATTERN) Long description
%
tmp = {fnames.name};
mask= zeros(numel(tmp),numel(pattern));

for el = 1:numel(tmp)
    mask(el,:) = ~cellfun(@isempty,regexp(tmp(el),pattern));
end

mask = logical(prod(mask,2));

fnames = fnames(mask);
end
