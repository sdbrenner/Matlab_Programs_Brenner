function Ki = loc_stiff_mat(vertices)
% Local stiffness matrix
% vertices given as [x1 y1; x2 y2; x3 y3]; 

d = 2;   % dimension
T = det([ones(1,d+1);vertices'])/2; % area of 1 triangle 
%     (half of the area of the corresponding parallelogram).
G = [ones(1,d+1);vertices'] \ [zeros(1,d);eye(d)];
Ki = T*G*G';

end