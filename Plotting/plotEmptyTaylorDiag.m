function [pax] = plotEmptyTaylorDiag( rhoLim,pax )
% PLOTTAYLORDIAG plots a Taylor Diagram
% 
%     plotTaylorDiag( X,Y ) calculates statistics and plots the Taylor
%     Diagram for reference (observation) series X, and test (model) series
%     Y. X must be a ROW vector (size 1xN); Y can be a vector of the same
%     size as X or a matrix where each ROW represents a separate test
%     series (size MxN).
%
% S.D.Brenner, 2021

%% Parse Inputs

% Check size/type ?

%% Transform statistics and define contours, axes, etc;

% Transform Correlations and normalize STD
% theta = real(acos(R));
% rho = stdY/stdX;

% Create contours of N-RMSD
thConts = linspace(0,pi/2,45);
NRMS = (1:0.5:rhoLim).';
clear rhConts
rhConts = cos(thConts) + sqrt( NRMS.^2 - sin(thConts).^2  );

% Correlation and STD labels
rLab = fliplr( [0:0.1:1] );
thLab = real(acos( rLab )) ;
rhoLab = 0:0.5:rhoLim;


%% Plot

if nargin < 2
    pax = polaraxes();
else
    axes(pax);
end

hold on;
rmsHand = polarplot( thConts, real(rhConts), '-','color',grey(0.7) );


pax.LineWidth = 1.5;
pax.ThetaLim = [0,90];
pax.ThetaTick = rad2deg(thLab);
pax.ThetaTickLabel = string(rLab);
pax.RTick = rhoLab;
pax.RLim = [0,rhoLim];

end
