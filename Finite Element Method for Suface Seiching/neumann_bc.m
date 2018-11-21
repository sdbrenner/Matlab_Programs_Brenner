function F = neumann_bc(p,e,g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS CODE IS INCOMPLETE AND WILL NOT CALCULATE THE EQUIVILANT FORCING  %
% FOR NEUMANN BC.                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Apply Neumann Boundary conditions (grad(u)= 0);
% p is the coordinates of the nodes
% e is the edges
% g is a function for the values of du/dn
% p can be the output from distmesh2d.m
% e can be the output from boundedges.m (if the entire domain has
% Neumann B.C.)
% Neumann B.C. are applied in FEM as additional force contributions on the
% right-hand-side of the equation


N=size(p,1); % number of nodes
E=size(e,1); % number of edges

F = sparse(N,N);

% Loop through edges
for n = 1:E;
  % get vertices for a given triangle:
    nodes = e(n,:); % row of e = node numbers of the 2 corners of edge
    vertices = p(nodes,:);
  % length of edge:
    E = norm(vertices);
  % local contribution:
    Fn = loc_stiff_mat(vertices);
  % add local contribution to global matrix:
    F(nodes,nodes) = F(nodes,nodes) + Fn;
end


end









