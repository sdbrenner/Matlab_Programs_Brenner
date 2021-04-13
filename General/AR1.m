function x = AR1(N,phi)
% AR1 generates a zero-mean AR1 red noise signal of length N
%
%   x = AR1(N,phi) where phi is the persistance parameter. To produce a
%   process with lag τ, phi = 1-Δt/τ, where Δt is the signal timestep. 
%
%   The AR1 process has a nominal mean of 0 and nominal variance of 1/(1-phi^2).
%
%   S.D.Brenner, 2021.


%% Parse inputs

    if abs(phi)>=1
        warning('Input phi should be less than 1; results may be unstable');
    end
    if N<1 || floor(N)~=N
        error('Input N must be a positive integer');
    end

%% Generate noise

    x = NaN(1,N);   
    x(1) = randn;
    for n = 2:N
        x(n) = phi* x(n-1) + randn ;
    end
    
end