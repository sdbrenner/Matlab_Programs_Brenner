function [P,bUnc,aUnc,yConf,r2] = linFitUnc(X,Y,alpha)
% LINFITUNC calculates a linear fit and associated slope uncertainty based on
% a two-sided T-test
%
%   [P,bUnc,aUnc,yConf,r2] = linFitUnc(X,Y,alpha)
%   alpha = 1-uncertainty level (i.e. for 95% uncertainty, alpha = 0.05)
%
%   S.D.Brenner, 2020
    

%% Calculate linear fit
[P] = polyfit(X,Y,1);
    


%% Calculate values for uncertainty/confidence tests
% (tests follow Bendat and Piersol, pp. 129-133 )

N = length(X);                          % number of predictors
YP = polyval(P,X);                      % predicted values of Y
XM = mean(X);                           % mean value of X
T = tinv(1-alpha/2, N-2);               % Student-t distribution statistics
syx = sqrt( ( 1/(N-2) ) * sum( (Y-YP).^2) );                        % eqn. 4.72

% Intercept uncertainty
aUnc = syx * T * sqrt( 1/N + XM^2./sum((X-XM).^2) );                % eqn. 4.69

% Slope uncertainty
bUnc = syx * T * sqrt( 1./sum((X-XM).^2) );                         % eqn. 4.70

% Define confidience interval function
yConf = @(X0) syx * T * sqrt( 1/N + (X0-XM).^2./sum((X-XM).^2) );   % eqn. 4.71

% Correlation coefficient
% (following Hartmann & Barnes; eqn. 4.37)
Xprime = X-mean(X);
Yprime = Y-mean(Y);
r2 = mean( Xprime.*Yprime )^2 ./ ( mean(Xprime.^2) * mean(Yprime.^2) );

%% Calculate slope uncertainty (old version, 
% following Hartmann and Barnes, pp. 43-44)
% 
% M = length(X);
% residuals = Y - polyval(P,X);             % residuals
% varRes =  (1/(M-2)) *sum( residuals.^2 ); % standard variance of residuals
% varX = 1/M*sum( (X-mean(X)).^2 );         % variance of X  
% stdB = sqrt( (1/M) * (varRes)/varX );     % standard deviation in slope
% bUnc = tinv(1-alpha/2, M-2) * stdB;       % uncertainty in slope (2-sided t-test)


end