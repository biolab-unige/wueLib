function varargout = bn_convertAllPLXtoEEG(varargin)
% bn_convertAllPLXtoEEG

% Edited 2015-04-02 by Gabriele Arnulfo <gabriele.arnulfo@gmail.com>
details = textscan(fopen('/mnt/data/BIOLAB_DATA/gabri/dbs/MER/plexon_files/LFPs_selected_plexon_file_code_list.csv')...
													,'%s %s\n','delimiter',',');

logFid = fopen('~/Desktop/Convert.log','w');

for files = 1:numel(details{1})
		warning('off','all');
	
  	path = '/mnt/data/BIOLAB_DATA/gabri/dbs/MER/plexon_files/LFPs';
		filename = strcat(details{1}{files},'.plx');

		if( strcmp(details{2}{files},'o') )
				disp(fullfile(path,filename));

				try 
						bn_convertPLXtoEEG(fullfile(path,filename));
				catch
						fprintf(logFid,strcat(filename,' : error while converting\n'));
						continue
				end

	
		end
end
 fclose(logFid);

end
