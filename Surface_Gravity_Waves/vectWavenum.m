function [k,err] = vectWavenum(omega,depth,tolerance,num_iterations)
% VECTWAVENUM Linear dispersion relation wavenumber solver
% 
%	k = vectWavenum(omega,depth) solves for the wave-numbers 'k' for a
%	given depth that satisfy the dispersion relation:
%
%   	omega^2 = g*k*tanh(k*depth)
%	
%   using the Newton-Raphson method.  The variables omega and depth can be
%	arrays.
% 
%	k = vectWavenum(omega,depth,tolerance) assigns a tolerance value for
%	this calculation. If 'tolerance' is omitted, or entered as [], then a
%	default value of 1e-6 is used.  The solver continues iterating until
%	all values of 'k' meet this tolerance level or until it reaches the
%	maximum number of iterations.
% 
%	k = vectWavenum(omega,depth,tolerance,num_iterations) computes 'k'
%   within a given tolerance using a maximum number of interations:
%   'num_iterations'. If 'num_iterations' is omitted, or entered as [],
%   then a default value of 10 is used.
% 
%	[k,err] = vectWavenum(omega,depth,...) returns the approximation error
%   associated with each value 'k': err = abs( g*k*tanh(k*d) - omega^2 )
% 
%	S.D.Brenner, April 2018

%% Error checking:
narginchk(2,4);
nargoutchk(1,2);

% omega and depth should be either scalars, or arrays of the same size:
if all(size(omega) ~= size(depth)) && length(omega(:)) ~= 1 && length(depth(:)) ~= 1
    error('''omega'' and ''depth'' must be arrays of the same size or be of length 1');
end

% If tolerance and num_iterations are included, they should be scalars:
if nargin >= 3 &&  length(tolerance) >1
   error('If ''tolerance'' is included, it should be length 1');
end
if nargin == 4 && length(num_iterations) >1
    error('If ''num_iterations'' is included, it should be length 1');
end

%% Set optional argument values
g = 9.81;
args = {1e-6,10}; % Default values
if nargin >= 3 && ~isempty(tolerance)
    args{1} = tolerance;
end
if nargin == 4 && ~isempty(num_iterations)
    args{2} = num_iterations;
end
tol = args{1};
num_iter = args{2};

% If one of omega or depth are scalars, they should be converted to vectors
% of equal size:
if length(omega(:)) == 1
    omega = omega* ones(size(depth));
end
if length(depth(:)) == 1
    depth = depth* ones(size(omega));
end


%% Initial guess:

% from CEM II-1-11:
% kn = (w.^2/g) .* ( tanh( (w.^2/g).*d ) ).^(-1/2)

% from Fenton and McKee, 1990:
k0 = omega.^2/g; % deep water wavenumber
kn = k0 .* ( tanh((k0.*depth).^(3/4)) ).^(-2/3);

% Error for initial guess
err = abs( omega.^2 - g*kn.*tanh(kn.*depth) ); 

%% Use Newton's method: k_(n+1) = k_n + f(k_n)/f'(k_n)

converged = 0; % flag to identify if the iteration has converged with tolerance
for n = 1:num_iter
    
    % on each loop, re-compute only the k's that don't meet the tolerance
    idx = find(err > tol); 
    
    fkn = g*kn(idx).*tanh(kn(idx).*depth(idx)) - omega(idx).^2;
    dfkn = g*tanh(kn(idx).*depth(idx)) + g*depth(idx).*kn(idx).*sech(kn(idx).*depth(idx)).^2;
    kn(idx) = kn(idx) - fkn./dfkn;
    
    % Estimation error
    err = abs( omega.^2 - g*kn.*tanh(kn.*depth) );
    
    % Tolerance check
    if max(err(:)) <= tol 
        converged = 1;
        break;
    end
end
if converged == 0
    warning('Method failed to a tolerance of %2.2d converge in %g iterations.  Maximum error of %2.2d recorded. Considering decreasing ''tolerance'' or increasing ''num_iterations''',tol, num_iter, max(err(:)));
end

% Return result:
k = kn;

% For any places with zero depth, set k to zero:
k(depth==0) = 0;


end
