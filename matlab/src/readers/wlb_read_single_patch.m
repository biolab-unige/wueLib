function vert = bn_read_single_patch(filename)
% bn_read_single_patch
%
% input: filename, the path to read the patch vertices from
% output: vert, 1 x N array of vertices contained in specified anatomical patch

% Edited 2013-08-09 by Gabriele Arnulfo <gabriele.arnulfo@gmail.com>

fid = fopen(filename);

fgetl(fid);
n_vert = fscanf(fid,'%d\n',1);

data = NaN(n_vert*5,1);

data = fscanf(fid,'%d %f %f %f %f\n',n_vert*5);
data = reshape(data,[5 n_vert]);

vert = data(1,:)+1';

fclose(fid);

end
