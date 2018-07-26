function [Q,cos_ths] = insolation(lat,lon,mattimeUTC)
% Top of atmosphere solar radiation input (insolation),
% calculated from Global Physical Climatology (Hartmann, 1994)
% https://www.elsevier.com/books/global-physical-climatology/hartmann/978-0-12-328531-7
%
% [Q,cos_ths] = insolation(lat,lon,mattimeUTC)
%
% Samuel Brenner, 2018

% convert from UTC time to local time:
timezone = lon/360;
mattime = mattimeUTC + timezone;

% convert from matlab time vector to a day-of-year time vector
[year,~,~,~,~,~] = datevec(mattime);
dn = mattime - datenum(year,1,1)+1;
doy = floor(dn);

%% Declination angle - Fourier series approximation (eq. A.6)
theta_d = 2*pi*dn/365; % (eq. A.5)
% Fourier coefficients:
a = [0.006918, -0.399912, -0.006758, -0.002697];
b = [0.000000,  0.070257,  0.000907,  0.001480];
N = 4;
fs = zeros([size(dn),N]); % pre-allocate
% Fourier series (written so lat, mattime can be scalars or arrays) 
for m = 1:N
    n = m-1;
    fs(:,:,m) = a(m).*cos(n*theta_d) + b(m).*sin(n*theta_d);
end
% sum
delta = sum(fs,3);


%% Squared ratio of the mean Earth-sun distance to the actual distance
% - Fourier series approximation (eq. A.7)
% Fourier coefficients:
a = [1.000110,  0.034221,  0.000719];
b = [0.000000,  0.001280,  0.000077];

N = length(n);
fs = zeros([size(dn),N]); % pre-allocate
% Fourier series (written so lat, mattime can be scalars or arrays) 
for m = 1:N
    n = m-1;
    fs(:,:,m) = a(m).*cos(n*theta_d) + b(m).*sin(n*theta_d);
end
dd2 = sum(fs,3);


%% Hour angle:
time = dn - floor(dn);   % local time of day (in units of [days])
h = 2*pi*( (time-0.5) ); % local time of day (in units of [rad]), solar noon is zero
% hour angle at sunrise and sunset:
phi = deg2rad(lat);
h0 = real(acos( -tan(phi) .* tan(delta))); % (eq. 2.16)
night = abs(h)  < h0; % zero during night time and one during daytime

%% Zenith angle
cos_ths = sin(phi).*sin(delta) + cos(phi).*cos(delta).*cos(h); % (eq. 2.15)

%% Solar flux
So = 1367; % Solar constant
Q = So.*dd2.*cos_ths.*night; % Solar flux (eq. 2.14)
Q(Q<0) = 0; % during polar night, Q will evaluate to negative.  This is unphysical.


%% Daily average solar flux (not used):
Qday = (So/pi)*dd2.*( h0.*sin(phi).*sin(delta) + cos(phi).*cos(delta).*sin(h0) ); % (eq. 2.17)
Qday(Qday<0) = 0;


end%function
