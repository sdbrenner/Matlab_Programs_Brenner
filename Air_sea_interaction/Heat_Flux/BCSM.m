function [Rs,alf_sky] = BCSM(R_et,Z,P)
%% BCSM Bird Clear Sky Model
%
%   [Rs,alf_sky] = BCSM(R_et,Z,P) calculates the transmitted downwelling
%   solar radiation, Rs, and atmospheric albedo, alf_sky, for given
%   top-of-atmosphere solar radiation, R_et, solar zenith angle, Z
%   [degrees], and barometric pressure, P [mb]
%   The units of Rs will be the same as the units of R_et.
%
%   Surface transmitted downwelling solar radiation, nominally based on the
%   procedure of Bird (1994), but equations are copied from:  
%       A General Lake Model (GLM 2.4) for linking with high-frequency
%       sensor data from the Global Lake Ecological Observatory Network
%       (GLEON).
%       Hipsey et. al., 2017 
%       doi: 10.5194/gmd-2017-257
%
%   ** Note: The Bird (1994) model uses a spectral parameterization while
%   Hipsey et. al. (2017) convert that into a bulk form.  In some cases it
%   is difficult to see whether the two versions match. Further, there are
%   some cases where model parameters used in GLM were unclear, so some
%   assumptions have been made have been made about these values.



%% Direct Irradiance

% Rayleigh scattering
P0 = 1013;
M = ( cosd(Z) + 0.15*(93.885 - Z).^(-1.253) ).^(-1); % relative air mass
Mp = M.*P/P0;
T_rayleigh = exp(-Mp);

% Aerosol attenuation
tauA = 0.21; %?
T_aerosol = exp( (-tauA.^0.873).*(1+tauA-tauA.^0.7088).*(M.^0.9108) );
% Water vapor absorbtion
a = 0.08;
b = 2.01;
Td = 10; %???
W = exp( a *Td + b);
Wm = W*Mp;
T_wv = 1 - 2.4959*Wm./( 1 + (79.034*Wm).^0.6828 + 6.385*Wm);
% Ozone absorbtion
Oz_C = 0.3;
Oz_M = M*Oz_C;
T_oz = 1-(0.1611*Oz_M.*(1+39.48*Oz_M).^(-0.3035))...
        - 0.002715*Oz_M./(1 + 0.044*Oz_M + 0.0003*Oz_M.^2);
% Uniformly mixed gas absorbtion
T_mix = exp( -0.0127 * Mp.^0.26 );

% Direct transmitance
R_dir = R_et .* T_rayleigh .* T_aerosol .* T_wv .* T_oz .* T_mix;

%% Diffuse Irradiance

% Aerosol scattering
T_aa = 1- (0.1*( 1 - M + M.^1.06 ).*(1-T_aerosol) );

% Diffuse transmitance
R_dif = R_et .* T_oz .* T_mix .* T_wv .* T_aa;


%% Total:
Rs = 0.9662 * R_dir .* cosd(Z) + 0.79 * R_dif;
alf_sky = 0.068 + (1-0.84)*(1- T_aerosol./T_aa );
Rs = Rs.*(1-alf_sky);

end%function