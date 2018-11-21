function Mi = loc_mass_mat(vertices)
% Local mass matrix
% vertices given as [x1 y1; x2 y2; x3 y3]; 

d = 2;   % dimension
T = det([ones(1,d+1);vertices'])/2; % area of 1 triangle 
%     (half of the area of the corresponding parallelogram).
Mi = (1/12)*T*(ones(3)+eye(3));

end