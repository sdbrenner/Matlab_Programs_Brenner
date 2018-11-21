function M = mass_mat(p,t)
% Assemble FEM Mass Matrix
% p is the coordinates of the nodes
% t is the triangles
% p and t can be outputs from distmesh2d.m


N=size(p,1); % number of nodes
T=size(t,1); % number of triangles

M = sparse(N,N);

% Loop through triangles
for n = 1:T
  % get vertices for a given triangle:
    nodes = t(n,:); % row of t = node numbers of the 3 corners of triangle
    vertices = p(nodes,:);
  % local contribution:
    Mn = loc_mass_mat(vertices);
  % add local contribution to global matrix:
    M(nodes,nodes) = M(nodes,nodes) + Mn;
end


end