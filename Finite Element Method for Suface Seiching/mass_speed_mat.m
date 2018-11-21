function MC = mass_speed_mat(p,t,cp)
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

MC = sparse(N,N);

% Loop through triangles
for n = 1:T
  % get vertices for a given triangle:
    nodes = t(n,:); % row of t = node numbers of the 3 corners of triangle
    vertices = p(nodes,:);
    Hn = mean(Hp(nodes));
  % local contribution:
    MCn = loc_mass_mat(vertices)*Hn;
  % add local contribution to global matrix:
    MC(nodes,nodes) = MC(nodes,nodes) + MCn;
end


end