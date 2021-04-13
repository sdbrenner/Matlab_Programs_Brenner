function output = PWP( input )
% PWP
%
%
%
%
%

% Model Implementation (§4.3):
%   1. Absorb solar insolation (eqn. 6)
%   2. Calculate air-sea heat lost and freshwater flux (E-P) from top-most
%      grid
%   3. Calculate density profile 
%       (a) Free convective mixing: Adjust for static stablity (eqn. 8) by
%           mixing from surface downward
%   4. Absorb wind stress within the mixed layer
%   5. Step the momentum equation forward in time
%   6. Calculate mixed-layer (bulk) Richardson number, Rb
%       (a) Mixed-layer entrainment: If Rb is less than a critical value,
%           entrain deeper water to bring Rb above the critical value (eqn
%           9; Rb_crit = 0.65)
%   7. Calculate gradient Richardson number, Rg, over the stratified part
%      of the profile (not in the mixed-layer).
%       (a) Gradient mixing: If Rg is less a critical value, partially mix
%           field variables rho, T, S, and V (Rg_crit = 0.25, but a value
%           of Rg = 0.3 is used in PWP, 1986 to help with convergence).
%       (b) After partial mixing, re-calculate Rg and continuing iterating
%
% Steps 1-6 correspond to the DIM model (Price et al., 1978).


%% Define model constants

% depth and times-steps
dz = 0.5;
dt = 900;

% Solar insolation constants
I1 = 0.62; 
I2 = 1-0.62;
lambda1 = 0.6;
lambda2 = 20;

% Critical Richardson numbers
Rb_crit = 0.65;
Rg_crit = 0.25;

%% Load input data

%% Organize data


%% Model Loop

for n = 1:N
   allFields(n+1) = pwpEvolve( allFields(n) , allConsts );  
end


%% Save output


end
%% EMBEDDED FUNCTIONS %% ==================================================


function allFields = pwpEvolve( allFields , allConsts )

    % Extract variables
    [Io,L,E,P,SA,CT,sigma,u,v,mld] = extractFields( allFields );
    [] = extractConsts( allConsts );
    rho = sigma + 1000;
    
    % Absorb solar radiation and extract heat loss:
    Iz = insolation( Io, z, allConsts );
    Lz = zeros(size(z));
    Lz(1) = L;
    Fz = (Iz - Lz)/dz;
    CT = CT + dt*( 1./(rho*cp) .* Fz ); 
    % Check for below-freezing values (and replace with in-situ freezing)
    % [ note: this step was not mentioned or included as part of the
    %   original PWP model; however, given the temperature range I work in,
    %   it is necessary to include ]
    CTfr = gsw_CT_freezing(SA,z);
    CT( CT<CTfr ) = CTfr( CT<CTfr );
    % Extract freshwater flux:
    Ez = zeros(size(z));
    Ez(1) = SA(1) * (E-P) / dz;
    SA = SA - dt*( Ez );

    % Check for static stability
    [SA,CT,sigma,u,v] = mixConv(SA,CT,sigma,u,v);
    
    % Absorb wind and step forward momentum equations
    % (i.e. apply coriolis force by rotating u, v)
    u = u + dt*(-f*v + Gxz  );
    v = v + dt*( f*u + Gyz  );
    
    % Mixed-layer entrainment (Rb)
    [SA,CT,sigma,u,v] = mixBulk(SA,CT,sigma,u,v,Rb_crit);
    
    % Gradient mixing (Rg)
    [SA,CT,sigma,u,v] = mixGrad(SA,CT,sigma,u,v,Rg_crit);

    % package variables
    allFields = packageFields( );
end




function extractFields( allFields )
end


function allFields = packageFields( )
end


function Iz = insolation( Io, z, allConsts )
% Find solar insolation profile
    [] = extractConsts( allConsts );

    Iz = Io*( I1*exp(-z/lambda1) + I2*exp(-z/lambda2) );

end



function [SA,CT,sigma,u,v] = mixConv(SA,CT,sigma,u,v)
% Mix convectively from the surface downwards until a stable density
% profile is acheived
    while 1
        sigma = gsw_sigma0( SA, CT );
        % find index of first stable bin:
        convDepthInd = find( diff(sigma)>=0, 1, 'first');
        % If the surface bin is stable, end loop
        if convDepthInd == 1
            break;
        % Otherwise, mix fields within the convectively unstable plume
        else
            SA(1:convDepthInd) = mean( SA(1:convDepthInd) );
            CT(1:convDepthInd) = mean( CT(1:convDepthInd) );
             u(1:convDepthInd) = mean(  u(1:convDepthInd) );
             v(1:convDepthInd) = mean(  v(1:convDepthInd) );
        end
    end
end



function [SA,CT,sigma,u,v] = mixBulk(SA,CT,sigma,u,v,Rb_crit)
% Entrain the mixed layer deeper to increase the bulk Richardson number
% above its critical threshold
    while 1
        hInd = find( (z-h)>0, 1, 'first');
        dr = sigma(hInd+1) - sigma(hInd-1);
        du = u(hInd+1) - u(hInd-1);
        dv = v(hInd+1) - v(hInd-1);
        dV2 = du.^2 + dv.^2;
        Rb = g*dr*h/(rho0*dV2);
        if Rb > Rb_crit
            break;
        % Otherwise, mix fields 
        else
            SA(1:hInd+1) = mean( SA(1:hInd+1) );
            CT(1:hInd+1) = mean( CT(1:hInd+1) );
             u(1:hInd+1) = mean(  u(1:hInd+1) );
             v(1:hInd+1) = mean(  v(1:hInd+1) );
             sigma = gsw_sigma0( SA, CT );
        end
    end
end

function [SA,CT,sigma,u,v] = mixGrad(SA,CT,sigma,u,v,Rg_crit)
end
