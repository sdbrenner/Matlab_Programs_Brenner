function [dlHours,sunDeclination] = daylightHours(lat,dateTime)
% DAYLIGHTHOURS calculates the number of daylight hours for a given
% latitude and time of year
%
%   [hrs] = daylightHours(lat,datetime)
%     
%   Without any inputs, daylightHours returns the daylight hours for the
%   current datetime at a latitude of 47.6062° N (corresponding to Seattle,
%   WA, USA).
%
%   This code adapted from the java function "SUNDIAL", by Chris Acky:
%       https://github.com/Kalyse/Sundial
%       Copyright (c) 2012, Chris Acky 
%
%
%   NOTES ON ACCURACY:
%
%   SUNDIAL claims to be accurate to 0.0001 minutes by taking into account
%   the axial tilt of the earth and the "equation of time". However, the 
%   result does not appear to depend on the equation of time.  Furthemore,
%   comparisons with NOAA's solar calculator show notable discepancies
%   For example, the java implentation of SUNDIAL indicates 16.15 hours of
%   daylight on May 21 2012 at 50°N, whereas the NOAA calculator shows 15.7
%   daylight hours for the same date and latitude.
%
%   Significant differences between the two calculators can be traced to
%   two main sources:
%       1. An error in the Julian date calculation in the java
%          implementation of SUNDIAL (resulting in a 31 day offset)
%       2. The inclusion of an estimated atmospheric refraction term in the
%          NOAA calculator (assumed 0.833°) that is not included in SUNDIAL
%   The former of these errors is corrected here, while the atmospheric
%   refraction term is still ignored.  Accounting for atmospheric
%   refraction requires some model or parameterization of those effects,
%   which is non-trivial.
%       For more information, see: 
%       https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
%
%
%   S.D.Brenner, 2019


% The SUNDIAL code (https://github.com/Kalyse/Sundial) includes the
% following credits and references:
%  
% http://lexikon.astronomie.info/zeitgleichung/   EOT 
% http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?bibcode=1989MNRAS.238.1529H&db_key=AST&page_ind=2&plate_select=NO&data_type=GIF&type=SCREEN_GIF&classic=YES
% http://code.google.com/p/eesim/source/browse/trunk/EnergySim/src/sim/_environment.py?spec=svn6&r=6
% http://mathforum.org/library/drmath/view/56478.html 
% http://www.jgiesen.de/elevaz/basics/meeus.htm
% http://www.ehow.com/how_8495097_calculate-sunrise-latitude.html
% http://www.jgiesen.de/javascript/Beispiele/TN_Applet/DayNight125d.java
% http://astro.unl.edu/classaction/animations/coordsmotion/daylighthoursexplorer.html
% http://www.neoprogrammics.com/nutations/Nutation_In_Longitude_And_RA.php



%% Parse inputs

% Set default values:
if nargin < 2 || isempty(dateTime)
    dateTime = datetime(datestr( now() ));   % need to convert to appropriate date system
end
if nargin < 1 || isempty(lat)
    lat = 47.6062;
end

% If 2D matrix inputs are given for datetime & lat, they must be "plaid".
% I will simply check if they have the same size
if all(size(lat)>1) || all(size(dateTime)>1)
    if ~isequal( size(lat),size(dateTime) )
        error('If lat and dateTime are 2D matrix inputs, they must be the same size and be "plaid"');
    end
end


%% Convert time to Julian days

julianDays = juliandate( dateTime );
t = (julianDays - 2451545.0) / 36525.0;


%% Make plaid matrices from lat and t

if isvector(lat) || isvector(t)
    [T,LAT] = meshgrid(t,lat);
else
    T = t;
    LAT = lat;
end

%% Get "Equation of time"
% The equation of time is the difference between apparent solar time and
% mean solar time. At any given instant, this difference will be the same
% for every observer on Earth.
% https://en.wikipedia.org/wiki/Equation_of_time
%
% Note: SUNDIAL seems to claim that increased accuracy is given by the
% inclusion of the equation of time in the calculation of daylight hours;
% however, this calculation is not anywhere else in the code.  I have
% included it in this function for posterity.

obliquity = getObliquity(T);      
[RA,sunDeclination] = getAscensionDeclination(T);
LS = getSunsMeanLongitude(T);
deltaPsi = getDeltaPSI(T);
E = LS - 0.0057183 - RA + deltaPsi.*cosd(obliquity);
if (E>5)
    E = E - 360.0;
end
E = E*4; % deg. to min
E = round(1000*E)/1000;



%% Get daylight hours

Nenner = cosd(LAT).*cosd(sunDeclination);
C = -sind(sunDeclination).*sind(LAT)./Nenner; 
C2 = C.^2;

dlHours = 90.0 - atand( C./sqrt(1 - C2) );
dlHours = 2.0*dlHours/15.0;
dlHours = round(dlHours*100)/100; 

dlHours( C > 1 ) = 0;
dlHours( C <-1 ) = 24;


end

%% EMBEDDED FUNCTIONS %% ==================================================

%% Get right ascension
function [rightAscension,sunDeclination] = getAscensionDeclination(T)
    
    L = getSunsMeanLongitude(T); % Calculate the mean longitude of the Sun		
    M = 357.52910 + 35999.05030*T - 0.0001559*T.^2 - 0.00000048*T.^3; % Mean anomoly of the Sun
    M = mod( M , 360 );		
    M( M<0 ) = M( M<0 ) + 360;
    
    C = (1.914600 - 0.004817*T - 0.000014*T.^2).*sind(M);
    C = C + (0.019993 - 0.000101*T).*sind(2*M);
    C = C + 0.000290*sind(3*M);	
    
    theta = L + C; % get true longitude of the Sun	
    
    obliquity = getObliquity(T);				
    obliquity = obliquity + 0.00256*cosd(125.04 - 1934.136*T);		
    lambda = theta - 0.00569 - 0.00478*sind(125.04 - 1934.136*T); % get apparent longitude of the Sun
    rightAscension = atan2d( cosd(obliquity).*sind(lambda), cosd(lambda) );				
    rightAscension( rightAscension<0 ) = rightAscension( rightAscension< 0 ) + 360;
   
    sunDeclination = asind( sind(obliquity) .* sind(lambda) );					
end

%% Get sun's mean longitude
function LS = getSunsMeanLongitude(T)
    LS = 280.46645 + 36000.76983*T + 0.0003032*T.^2;
    LS = mod( LS, 360 );
    LS(LS<0) = LS(LS<0)+360;  
end

%% Get "deltaPSI"
% Nutation in ecliptical longitude expressed in degrees.
function deltaPsi = getDeltaPSI(T)

    LS = getSunsMeanLongitude(T);
    LM = 218.3165 + 481267.8813*T;
    LM = mod( LM, 360);
    LM( LM<0 ) = LM( LM<0 ) + 360;

    % Longitude of ascending node of lunar orbit on the ecliptic as measured
    % from the mean equinox of date.
    omega = 125.04452 - 1934.136261*T + 0.0020708*T.^2 + T.^3/450000;
    deltaPsi = -17.2*sind(omega) - 1.32*sind(2*LS) - 0.23*sind(2*LM) + 0.21*sind(2*omega);
    deltaPsi = deltaPsi/3600.0;
end

%% Get Obliquity
% Calculate Earths Obliquity Nutation (axial tilt of earth)
% ( T = Time Factor Time factor in Julian centuries reckoned from J2000.0,
%   corresponding to JD )
function obliquity = getObliquity(T)

LS = getSunsMeanLongitude(T);
LM = 218.3165 + 481267.8813*T;	
eps0 =  23.0 + 26.0/60.0 + 21.448/3600.0 - (46.8150*T + 0.00059*T.^2 - 0.001813*T.^3)/3600;
omega = 125.04452 - 1934.136261*T + 0.0020708*T.^2 + T.^3/450000;	
deltaEps = (9.20*cosd(omega) + 0.57*cosd(2*LS) + 0.10*cosd(2*LM) - 0.09*cosd(2*omega))/3600;

obliquity = eps0 + deltaEps;	

end