function [theta,rho] = taylorDiagStats( X,Y )
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

% Check sizes ?
% Remove NaNs ?


%% Calculate statistics

N = sum(~isnan(Y),2);

% Means and anomalies
Xbar = nanmean(X);
Ybar = nanmean(Y,2);
Xp = X-Xbar;
Yp = Y-Ybar;

% Standard deviations (normalized by N)
stdX = nanstd(X,1,2);
stdY = nanstd(Y,1,2);

% Correlations
R = 1./N .* nansum( Yp.*Xp,2 ) ./ (stdX*stdY);

% (centered) Root-mean-square deviation
E = sqrt( 1./N .* nansum((Yp-Xp).^2,2) ); 

% Check that statistics add up
tol = 1e-6;
err = (stdY./stdX).^2 + 1 - 2*(stdY./stdX).*R - (E./stdX).^2;
if err>tol
    warning('statistics mismatch -- results may be incorrect');
end


%% Transform statistics and define contours, axes, etc;

% Transform Correlations and normalize STD
theta = real(acos(R));
rho = stdY/stdX;


end
