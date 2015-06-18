function [vert faces] = bn_read_flat_surf(filename)
% bn_read_flat_surf
%
% input, filename
% output, [ vert, faces] 

% file formats reverse engineered information
% #!ascii version of blabla [USELESS INFO]
% number_of_vertices number_of_faces

% Edited 2011-11-21 by Gabriele Arnulfo <gabriele.arnulfo@gmail.com>

[fid, message] = fopen(filename,'r');

if fid < 0
	error('MATLAB:bn_read_flat_surf:FOPENErr',strcat(filename,message));
end

% read dimensions
dim = fscanf(fid,'%*72c\n%d %d\n',2);

% allocate memory
vertices    = NaN(dim(1)*3,1);
faces		= NaN(dim(2)*4,1);

% read vertices
vertices    = fscanf(fid,'%d vno=%*d\n%f %f %*f\n',dim(1)*3);
vertices 	= reshape(vertices,3,dim(1))';

% reallocate vertices
vert		= NaN(max(vertices(:,1)),2);
vert(abs(vertices(:,1)),:) = vertices(:,2:end);

fgetl(fid);

% read faces
faces		= fscanf(fid,'%d\n %d %d %d\n',dim(2)*4);
faces		= reshape(faces,4,dim(2))'+ 1;

end
