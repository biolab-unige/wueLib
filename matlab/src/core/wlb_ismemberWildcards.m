function mask = ismemberWildcards(stringsIn, patterns)
%ISMEMBERWILDCARDS Description
%	MASK = ISMEMBERWILDCARDS(STRINGSIN, PATTERNS) Long description
%
	
	nStrings = numel(stringsIn);
	nPatterns= numel(patterns);

	mask = cellfun(@(x) regexp(x,patterns),stringsIn,'Uni',false);
	mask = [mask{:}];
	mask = reshape(mask,[nPatterns,nStrings]);
	mask(cellfun(@isempty,mask)) = {0};

	mask = cell2mat(mask);

	mask = sum(mask) >= 1;
end
