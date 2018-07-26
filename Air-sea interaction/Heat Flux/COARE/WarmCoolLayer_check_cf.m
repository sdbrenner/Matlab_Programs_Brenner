% Program written to incorporate CF warm layer code in coare35vn bulk
% algorithm
% Ludovic Bariteau & Jim Edson
% 1st version: March 2013

%clear all;
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
U=U10;              %true wind speed
Pair=Pair10;        %pressure, mb
Tair=T10;           %air temperature, C
RH=RH10;            %relative humidity, %
Tsea=Tsea;          %bulk sea temperature, C
ts_depth=.05;       %bulk sea temperature sensor depth
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

B=coare35vnWarm(yday,U,zu,Tair,zt,RH,zq,Pair,Tsea,Solar,IR,Lat,Lon,zi,Rainrate,ts_depth);

%Outputs
%B(40)=dt_wrm;   %warming across entire warm layer deg.C
%B(41)=tk_pwp;   %warm layer thickness m
%B(42)=dsea;     %heating at selected depth
dt_wrmx=B(:,40)';
tk_pwpx=B(:,41)';
dseax=B(:,42)';

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

tsg_depthx=2;                      %Approximate depth of TSG
dseagx=dt_wrmx.*tsg_depthx./tk_pwpx;  %Linear heating profile is assumed
i=find(tk_pwpx<tsg_depthx);         %Use the total warming if the sensor
dseagx(i)=dt_wrmx(i);               %is below the warm layer.
tsgx=Tsea5+dseagx;                  %Add warm layer correction to adjust
tsx=Tsea+dseax;                     %temps to just below surface values.
dtwgx=dt_wrmx-dseagx;

figure(1);
plot(yday,Tsea,'r--',yday,Tsea5,'b--')
hold
plot(yday,tsx,'r',yday,tsgx,'b','linewidth',2)
axis([315 325 29 33])
legend('Sea Snake','ThermoSalinoGraph','Snk Corrected','TSG Corrected',1)
hold

jdy=Yday;
zu_etl=zu*ones(length(jdy),1);
zt_etl=zt*ones(length(jdy),1);
zq_etl=zq*ones(length(jdy),1);
%Lon=lon;Lat=lat;
Ts=Tsea;   
tsg=Tsea5;tsg_depth=2;%Use sea snake for water T, depth of snake is 0.05 m
ta=Tair;
tsnk=Tsea;
[Press,Tseak,Tairk,qse,Qsat,qa,Rhoair,Rhodry]=scalarv(Pair,Tsea,Tair,RH);
U=U10;
rs=Solar;
rl=IR;
tdk=272.15;
zts=ts_depth*ones(length(jdy),1);
Rgas=287.1;
barpress=Pair;
rhoa=barpress*100./(Rgas*(ta+tdk).*(1+0.61.*qa/1000));
%Hsb=hsb;Hlb=hlb;Taub=taub;usb=sqrt(taub./rhoa);
reldir=cdir;
org=Rainrate;
warm_calc_dyn35;
dtwrmx=dt;

figure;
subplot(2,1,1);plot(jdy,tsnk-Tsea5,jdy,dtw,'.',jdy,B(:,40),'o');xlim([315 325]);ylabel('DT (C)');
legend('Tsnk-Tsg','DT_{old}','DT_{Edson}');
subplot(2,1,2);plot(jdy,-(dtw-B(:,40)),'.');xlim([315 325]);ylabel('DT_{Edsond}-DT_{old}');
xlabel('Yearday (2011)');

figure;
subplot(3,1,1);plot(jdl,tsnk-tsb,jdl,dtw,'.',jdl,dt_wrmx,'.');xlim([315 325]);ylabel('DT (C)');ylim([-.5 3.5])
legend('Tsnk-Tsb','DT_{old}','DT_{Edson}');
subplot(3,1,2);plot(jdl,tsg-tsb,jdl,dtgw,'.',jdl,dtwgx,'.');xlim([315 325]);ylabel('DT_{sg} (C)');ylim([-.5 1])
subplot(3,1,3);plot(jdl,U);xlim([315 325]);ylabel('U_{10} (m/s)');
xlabel('Local Yearday (2011)');
