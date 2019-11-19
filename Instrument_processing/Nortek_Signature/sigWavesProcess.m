function  [sigWaves,astWaves,leWaves,presWaves]=sigWavesProcess(Data,waterdepth,lon,lat)
% SIGWAVESPROCESS Create wave estimates from Signature-series burst data.
%
% NOTE: THIS FUNCTION IS CURRENTLY A WORK IN PROGRESS
%
%   sigWaves = sigWaves(Data,waterdepth,lat,lon) creates wave estimates
%   through a combination of altimeter (AST and LE), and pressure burst
%   data. Results are stored in a 'SWIFT-compatible" data structure. Works
%   on a single burst at a time; application to longer data sets
%	could be performed in a loop.
%
%   Actual wave processing is done using the UVZwaves subroutine from the
%   SWIFT codes (https://github.com/jthomson-apluw/SWIFT-codes). Program
%   adapted from an AWAC code written by J. Thomson:
%   'AWACwaves_readandprocess.m'
%
%   Notes:  
%   (1) This function is developed to operate on Data structures that are
%   output by converting raw .ad2cp data to .mat files using MIDAS
%   software.  Data converted with Signature Deployment software may not
%   have matching variable names.
%   (2) Assumes that beam velocities have been converted to cartesian
%   Earth-referenced velocities ENU
%
%   S.D.Brenner, 2019

%% Define function constants:

astQualityCutoff = 4500;
leQualityCutoff = 5500;
despike = false;                % apply phase-space despiking to raw data
extrapEquilibriumRange = false; % extrapolate the pressure spectra when beyond noise floor
declination = 0;             % deg (positive east)
finalScreening = true;          % optional final screening of results
maxWavePeriod = 16;             % max wave period allowed during final screening

%% Extract simple values from Data structure

time = Data.Burst_MatlabTimeStamp(1);
pres = Data.Burst_AltimeterPressure;
ast = Data.Burst_AltimeterAST;
le = Data.Burst_AltimeterLE;
astQual = Data.Burst_AltimeterQualityAST;
leQual = Data.Burst_AltimeterQualityLE;
% heading = Data.Burst_Heading;
% pitch = Data.Burst_Pitch;
% roll = Data.Burst_Roll;

depth = mean(pres);

% Calculate sampling frequency
ts = mean(diff(Data.Burst_MatlabTimeStamp)) * 86400; % sampling period  [sec]
fs = 1/ts;                     % sample frequency [Hz]

%% Quality control the altimeter data:
%  Check that there are enough good data left for a reasonable wave
%  estimate.  If so, replace low quality values with record mean
numCutoff = 1600;

% Acoustic surface tracking
astBadBool = astQual < astQualityCutoff;
ast( astBadBool ) = mean( ast(~astBadBool) );

% Leading edge
leBadBool = leQual < leQualityCutoff;
le( leBadBool ) = mean( le(~leBadBool) );

% If length is too short, discard
if length(ast(~astBadBool)) < numCutoff &&...
   length(le(~leBadBool)) < numCutoff 
    
    sigWaves.time = [NaN];
    sigWaves.lat = [NaN];
    sigWaves.lon = [NaN];
    sigWaves.sigwaveheight = [NaN];
    sigWaves.peakwaveperiod = [NaN];
    sigWaves.peakwavedirT = [NaN];
    sigWaves.wavespectra.energy = [NaN];
    sigWaves.wavespectra.freq = [NaN];
    sigWaves.wavespectra.a1 = [NaN];
    sigWaves.wavespectra.b1 = [NaN];
    sigWaves.wavespectra.a2 = [NaN];
    sigWaves.wavespectra.b2 = [NaN];
    sigWaves.wavespectra.check = [NaN];
    return;
end


%% Extract velocity from Data structure
% From the "Signature Principles of Operation" manual (p.27):
%
%   From the recorded profile one may select a level below the surface
%   where the measurements form an array projected from the Signature to
%   just below the surface. Managing the fact that orbital velocities
%   attenuate exponentially with depth means that the data used for wave
%   processing are the ones that are measured close to the surface, while
%   ensuring that there is no contamination from the surface either
%   directly from the cells touching the surface or indirectly from
%   sidelobes. This can be managed by adaptively positioning the cells
%   just below the surface by a fraction of the measured depth; 10% of the
%   depth has proven to provide a good signal response without contamination.

% Find appropriate depth bin:
 
% % slight dynamic adjustment:
% ranges = Data.Burst_Range;
% rDepths = pres - ranges;
% targetDepth = 0.1*pres;
% [~,minInd] = min( abs(rDepths-targetDepth), [], 2 );
% 
% % Loop through and extract appropriate velocity data
% % [I feel like this can been done better with vector indices...]
% u = NaN(1,length(minInd));
% v = NaN(1,length(minInd));
% w = NaN(1,length(minInd));
% for n = 1:length(minInd)
%     ind = minInd(n);
%     u(n) = Data.Burst_VelEast(n,ind);
%     v(n) = Data.Burst_VelNorth(n,ind);
%     w(n) = Data.Burst_VelUp(n,ind);
% end

% Fixed depth (more appropriate for single burst)
ranges = Data.Burst_Range;
rDepths = depth - ranges;
targetDepth = 0.1*depth;
[~,minInd] = min( abs(rDepths-targetDepth) );

u = Data.Burst_VelEast(:,minInd);
v = Data.Burst_VelNorth(:,minInd);
try 
    w = Data.Burst_VelUp(:,minInd);
catch
    w = Data.Burst_VelUp1(:,minInd);
end



%% Despike

% Despike if requested
if despike
    % func_despike_phasespace3d uses 'cubic' for interp method and matlab will
    % throw warnings that 'pchip' should be used instead.  Supress those
    % warnings:
    wId = 'MATLAB:interp1:UsePCHIP';
    warning('off',wId);

    [ u, v, w, spikeinds] = func_despike_phasespace3d_3var( u, v, w, 2 );
    
    if ~all(isnan(ast))
        [ ast spikeinds] = func_despike_phasespace3d( ast, 0, 2 );
    end
    if ~all(isnan(le))
        [ le spikeinds] = func_despike_phasespace3d( le, 0, 2 );
    end
end

%% Rotate velocities for declination
% east = u;
% north = v;
% 
% u = east.*cosd(declination) - north.*sind(declination);
% v = east.*sind(declination) + north.*cosd(declination);

%% Waves proceess with AST

[ Hs, Tp, Dp, E, f, a1, b1, a2, b2, check] = UVZwaves(u, v, ast, fs);

astWaves.time = time;
astWaves.lat = lat;
astWaves.lon = lon;
astWaves.sigwaveheight = Hs;
astWaves.peakwaveperiod = Tp;
astWaves.peakwavedirT = Dp;
astWaves.wavespectra.energy = E';
astWaves.wavespectra.freq = f';
astWaves.wavespectra.a1 = a1';
astWaves.wavespectra.b1 = b1';
astWaves.wavespectra.a2 = a2';
astWaves.wavespectra.b2 = b2';
astWaves.wavespectra.check = check';

%% Waves proceess with LE

[ Hs, Tp, Dp, E, f, a1, b1, a2, b2, check] = UVZwaves(u, v, le, fs);

leWaves.time = time;
leWaves.lat = lat;
leWaves.lon = lon;
leWaves.sigwaveheight = Hs;
leWaves.peakwaveperiod = Tp;
leWaves.peakwavedirT = Dp;
leWaves.wavespectra.energy = E';
leWaves.wavespectra.freq = f';
leWaves.wavespectra.a1 = a1';
leWaves.wavespectra.b1 = b1';
leWaves.wavespectra.a2 = a2';
leWaves.wavespectra.b2 = b2';
leWaves.wavespectra.check = check';

%% Waves process with pressure

[ Hs, Tp, Dp, E, f, a1, b1, a2, b2, check] = UVZwaves(u, v, pres, fs);

presWaves.time = time;
presWaves.lat = lat;
presWaves.lon = lon;
presWaves.sigwaveheight = Hs;
presWaves.peakwaveperiod = Tp;
presWaves.peakwavedirT = Dp;
presWaves.wavespectra.energy = E';
presWaves.wavespectra.freq = f';
presWaves.wavespectra.a1 = a1';
presWaves.wavespectra.b1 = b1';
presWaves.wavespectra.a2 = a2';
presWaves.wavespectra.b2 = b2';
presWaves.wavespectra.check = check';

%% Adjust pressure waves:

% correct pressure for Bernouli (velocity registering as additional pressure)
E = E./4;

% correct pressure spectra for depth attenutation

%   find wavenumber for each frequency
% for j=1:length(f)
%     k(j) = wavenumber( f(j), waterdepth );
% end
k = vect_wavenum( 2*pi*f, waterdepth);


% transfer function
attenuation = cosh( k .* waterdepth ) ./ cosh( k.*(waterdepth-depth) ) ;
attenuation = attenuation.^2; % square for energy

% Adjust transfer funtion:
% limit the size of the attenuation correction  (to not amplify noise)
noise = attenuation > 100 | isnan(attenuation); % limit the size of the attenuation correction
attenuation( noise ) = NaN; % cut it off when correction too big, don't amplify noise

% Apply transfer function
E = E.*attenuation;

% Extrapolate the equilibruim range (if requested)
if extrapEquilibriumRange 
    E( noise ) = ( E( min(noise) - 1 ) .* f( min(noise) - 1 ).^4 ) .* f(noise).^-4; % extrapolate equilibrium range
end

% Ajust signficant wave height and spectra in structure
Hs = 4 * sqrt( nansum(E) * median(diff(f)) );
presWaves.sigwaveheight = Hs;
presWaves.wavespectra.energy = E';


%% make a final product by merging pressure (better at low freq) and AST (better at high freq)

sigWaves = astWaves;

% Use spectra from pressure-based processing outside of noise range
sigWaves.wavespectra.energy = presWaves.wavespectra.energy;
noise = find( isnan( presWaves.wavespectra.energy ) );

% use minimun available spectral energy density at high frequencies  
sigWaves.wavespectra.energy(noise) =...
    min( [astWaves.wavespectra.energy(noise) ,...
          leWaves.wavespectra.energy(noise) ],[],2);
  
% Re-calculate wave properties from merged spectra
sigWaves.sigwaveheight = 4 * ( nansum( sigWaves.wavespectra.energy .* median(diff(sigWaves.wavespectra.freq)) ) )^.5;
[~,peakindex] = max( sigWaves.wavespectra.energy );
sigWaves.peakwaveperiod = f(peakindex).^-1;


%% final screening of results

if finalScreening
    if sigWaves.peakwaveperiod > maxWavePeriod || sigWaves.sigwaveheight < 0.1
        sigWaves.sigwaveheight = NaN;
        sigWaves.peakwaveperiod = NaN;
        sigWaves.peakwavedirT = NaN;
    end
end



end