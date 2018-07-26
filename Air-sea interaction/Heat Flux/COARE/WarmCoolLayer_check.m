% Program written to incorporate CF warm layer code in coare35vn bulk
% algorithm
% Ludovic Bariteau & Jim Edson
% 1st version: March 2013

clear all;
%close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grab the 10 minute file from DYNAMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd c:\data\cwf\dropbox\matlabstf\cwf\bulkalg\cor3_5;
load Revelle10minutesLeg3_r3; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs Required
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Yday           %decimal yearday
%Ur             %relative wind speed, m/s
%Tsea           %Sea temperature, C at ts_depth 
%Tair           %air temperature, C
%Pair           %air temperature, C
%RH             %relative humidity (%)
%Lat            %latitude, deg
%Lon;           %longitude, deg
%Solar          %downward solar flux, W/m^2
%IR             %downward IR flux, W/m^2
%Rainrate       %rainrate, mm/hr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zu=10;              %anemometer height
zt=10;              %air T height
zq=10;              %humidity height

Yday=yday;          %decimal yearday Jan 1 at noon = 1.5
Ur=Ur10;            %relative wind speed
Pair=Pair10;        %pressure, mb
Tair=T10;           %air temperature, C
RH=RH10;            %relative humidity, %
Tsea=Tsea;          %bulk sea temperature, C
ts_depth=0.05;      %bulk sea temperature sensor depth
Solar=-Solardn;     %downward solar flux, W/m^2
IR=-IRdn;           %downward IR flux, W/m^2
Rainrate=P;         %rainrate, mm/hr
Lat=Lat;            %latitude, deg
Lon=Lon;            %longitude, deg
zi=600;             %inversion ht

%************************************************************
%  This test program uses the sea snake temperature at
%  approximately 5 cm to compute the depth and amplitude of 
%  the warm layer.  The heating at this depth is added to 
%  sea snake and the cool-skin affect is removed in the 
%  coare35vn algorithm.  Therefore, fluxes in the returned 
%  B array have been computed based on a Tsea measurements 
%  that have been corrected for the warm layer and cool skin.
%************************************************************

B=coare35vnWarm(yday,Ur,zu,Tair,zt,RH,zq,Pair,Tsea,Solar,IR,Lat,Lon,zi,Rainrate,ts_depth);

%Outputs
%B(40)=dt_wrm;   %warming across entire warm layer deg.C
%B(41)=tk_pwp;   %warm layer thickness m
%B(42)=dsea;     %heating at selected depth
dt_wrm=B(:,40)';
tk_pwp=B(:,41)';
dsea=B(:,42)';

%************************************************************
%  The following shows how to use the warm layer variables to
%  compare the temperature measured at another depth, i.e.,
%  the TSG. 
%  If the sensor is beneath the warm layer depth, the total 
%  temperature change (the warm layer amplitude) is added to 
%  the sensor. 
%  Note that these are the values just below the surface 
%  without the cool skin correction.
%************************************************************

tsg_depth=4;                      %Approximate depth of TSG
dseag=dt_wrm.*tsg_depth./tk_pwp;  %Linear heating profile is assumed
i=find(tk_pwp<tsg_depth);         %Use the total warming if the sensor
dseag(i)=dt_wrm(i);               %is below the warm layer.

tsg=Tsea5+dseag;                  %Add warm layer correction to adjust
ts=Tsea+dsea;                     %temps to just below surface values.

figure(1);clf
plot(yday,Tsea,'r--',yday,Tsea5,'b--')
hold
plot(yday,ts,'r',yday,tsg,'b','linewidth',2)
axis([315 325 29 33])
legend('Sea Snake','ThermoSalinoGraph','Snk Corrected','TSG Corrected',1)
hold

