function K = Hstiff_mat(p,t,Hp)
% Assemble FEM Stiffness Matrix
% p is the coordinates of the nodes
% t is the triangles
% p and t can be outputs from distmesh2d.m


N=size(p,1); % number of nodes
T=size(t,1); % number of triangles

K = sparse(N,N);

% Loop through triangles
for n = 1:T
  % get vertices for a given triangle:
    nodes = t(n,:); % row of t = node numbers of the 3 corners of triangle
    vertices = p(nodes,:); % coordinates of the 3 corners
    Hn = mean(Hp(nodes));  % average depth of the triangle
  % local contribution:
    Kn = Hn*loc_stiff_mat(vertices);
  % add local contribution to global matrix:
    K(nodes,nodes) = K(nodes,nodes) + Kn;
end


end