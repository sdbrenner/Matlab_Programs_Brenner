function M = Hmass_mat(p,t,Hp)
% Assemble FEM Mass Matrix
% p is the coordinates of the nodes
% t is the triangles
% p and t can be outputs from distmesh2d.m
% spd is a 3D matrix of wave speed c = c(x,y) where 
%   spd(:,:,1) are the x's,
%   spd(:,:,2) are the y's,
%   spd(:,:,3) are the c(x,y)'s,


N=size(p,1); % number of nodes
T=size(t,1); % number of triangles

M = sparse(N,N);

% Loop through triangles
for n = 1:T
  % get vertices for a given triangle:
    nodes = t(n,:); % row of t = node numbers of the 3 corners of triangle
    vertices = p(nodes,:);
    Hn = mean(Hp(nodes));
  % local contribution:
    Mn = Hn*loc_mass_mat(vertices);
  % add local contribution to global matrix:
    M(nodes,nodes) = M(nodes,nodes) + Mn;
end


end