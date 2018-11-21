function fh = hnonuniform(p, fd)
%    We want small triangles near the boundaries of the domain, and
%    larger ones in the interior.  To achieve this, we assume that P
%    represents a large sampling of points in the region; we compute
%    the minimum and maximum distances of the points to the boundary,
%    and we assign mesh density values of HMIN to the closest points,
%    HMAX to the furthest ones, and linearly vary H between them. 
%
%    Note that the points inside the region have negative signed distance,
%    and those furthest from the boundary have the most negative value.
%    Thus, we take the absolute value of this distance to get the positive
%    distance we would prefer to work with.  The program does not expect to
%    receive input points which are actually outside the region.
%

%    Modified from code provided by John Burkardt (21 August 2008) under
%    the GNU LGPL license.
%    http://people.sc.fsu.edu/~jburkardt/m_src/distmesh/p18_nonuniform_fh.m


%    Input, real P(N,2), one or more points, where the mesh density function
%    is to be evaluated.
%
%    Output, real H(N), the value of the mesh density function H(P).
%
  hmax = 5;
  hmin = 1;
  d = abs( fd(p) );
  
  dmax = max(d);
  dmin = min(d);

  fh = (  (dmax - d)*hmin +(d - dmin)*hmax  )/ (dmax - dmin);
%   fh = ((hmax-hmin)/(dmax-dmin))*(d-dmin) + hmin;  
 

end