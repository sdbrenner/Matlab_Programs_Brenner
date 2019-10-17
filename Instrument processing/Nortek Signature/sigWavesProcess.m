function  [sigWaves,astWaves,leWaves,presWaves]=sigWavesProcess(Data,waterdepth,lon,lat)
%   [sigWaves]=sigWaves(Data,waterdepth,lat,lon)
%
% Create wave estimates from signature-500 burst data.
% Assumes data conversion has been done in MIDAS, and velocities have
% already been converted to East-North-Up.
%
% Variable names modified to account for data-preprocessing steps that I've 
% performed on the raw data; specifically, I've separated data types so 
% instead of 'Burst.MatlabTimeStamp' (and similar variables), the 
% variable is now: 'Burst.MatlabTimeStamp'.  To return to default names,
% simply find-and-replace 'Burst.' with 'Burst.'.
%
% Works on a single burst at a time.  Application to longer data sets
% could be performed in a loop.
%
% Calls a sub routine from the SWIFT codes for the actual processing
% and stores results in a SWIFT compatible structure.
% Program adapted from an AWAC code written by J. Thomson:
% 'AWACwaves_readandprocess.m'
%
% S.Brenner, 2019

%% Define function constants:
astQualityCutoff = 4500; %120;
leQualityCutoff = 5500;
despike = false;                % apply phase-space despiking to raw data
extrapEquilibriumRange = false; % extrapolate the pressure spectra when beyond noise floor
declination = 15.6;             % deg (positive east)
finalScreening = true;          % optional final screening of results
maxWavePeriod = 16;             % max wave period allowed during final screening

%% Extract simple values from Data structure

time = Burst.MatlabTimeStamp(1);
pres = Burst.Pressure;
ast = Burst.AltimeterAST;
le = Burst.AltimeterLE;
astQual = Burst.AltimeterQualityAST;
leQual = Burst.AltimeterQualityLE;
% heading = Burst.Heading;
% pitch = Burst.Pitch;
% roll = Burst.Roll;

depth = mean(pres);

% Calculate sampling frequency
ts = mean(diff(Burst.MatlabTimeStamp)) * 86400; % sampling period  [sec]
fs = 1/ts;                     % sample frequency [Hz]


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
% ranges = Burst.Range;
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
%     u(n) = Burst.VelEast(n,ind);
%     v(n) = Burst.VelNorth(n,ind);
%     w(n) = Burst.VelUp(n,ind);
% end

% Fixed depth (more appropriate for single burst)
ranges = Burst.Range;
rDepths = depth - ranges;
targetDepth = 0.1*depth;
[~,minInd] = min( abs(rDepths-targetDepth) );

u = Burst.VelEast(:,minInd);
v = Burst.VelNorth(:,minInd);
try 
    w = Burst.VelUp(:,minInd);
catch
    w = Burst.VelUp1(:,minInd);
end

%% Quality control the altimeter data:
%  replace low quality values with record mean

% Acoustic surface tracking
badBool = astQual < astQualityCutoff;
ast( badBool ) = mean( ast(~badBool) );

% Leading edge
badBool = leQual < leQualityCutoff;
le( badBool ) = mean( le(~badBool) );

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