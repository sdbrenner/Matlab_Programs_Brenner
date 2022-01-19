function [isIce,hfVar,lfVar,S,f] =  spectralIceIdentification(Data,fs,fc,rc)
% spectralIceIdentification    Identifies ice and open-water bursts based
% on spectral energy partion (following Scherbina et al, 2016).
%
%   Ice is classified based on the relative partion of energy between a low
%   frequency band and a high frequency band. Ice-cover typically has a red
%   spectrum, so will contain energy primarily in low frequencies. Energy
%   in the high-frequency band corresponds to surface gravity waves in open
%   water. This scheme may have trouble differentiating between smooth ice
%   or calm open-water surface conditions, which will both have a flat,
%   low-energy spectra.
%
%
%   isIce = spectralIceIdentification(altimeter,fs) returns ice
%   classification isIce for altimeter data sampled at a frequency fs
%   [Hz].  When altimeter is a vector, isIce is a scalar value. When
%   altimeter is a matrix, spectralIceIdentification is carried out for
%   each column of the matrix and isIce is a vector with
%   corrsponding values. Columns of altimeter containing NaN values will
%   return NaN values for isIce. Input altimeter should be despike prior to
%   use of this function.
%
%   isIce = spectralIceIdentification(altimeter,fs,fc,rc) allows
%   specification of the frequency cutoff (fc) and variance ratio cutoff
%   (rc) used for separating high- and low-frequnecy data and for
%   designating ice conditions.  Default values are fc = 0.1 Hz, rc = 5.
%
%   [isIce,hfVar,lfVar] =  spectralIceIdentification(...) returns the
%   variance within the high frequency band (hfVar) and low frequency band
%   (lfVar).
%
%   [isIce,hfVar,lfVar,S,f] =  spectralIceIdentification(...) also returns
%   the frequency spectrum S and correspoding frequencies, f.
%
%   S.D.Brenner, 2020 

    % Parse inputs
    fCutoffDefault = 1/10; %0.1;
    ratioCutoffDefault = 5;
    fN = fs/2;

    % Remove incomplete bursts
    nanInd = all(~isnan( Data),1);
    Z = Data(:,nanInd);
    
    % Calculate spectra
    [S,f] = pwelch( detrend(Z),[128],[],[],fs);
       

    % Identify high- and low-frequency ranges
    hfInd = find( f > fc  & f < 0.9*fN );
    lfInd = find( f <= fc & f>0 );

    % Calculate variance in each range
    hfVar = NaN(1,size(Data,1));
    lfVar = NaN(1,size(Data,1)); 
    hfVar(:,nanInd) = trapz(f(hfInd), S(hfInd,:));
    lfVar(:,nanInd) = trapz(f(lfInd), S(lfInd,:));

    % Identify ice based on variance ratio
    varianceRatio = hfVar./lfVar;
    isIce = varianceRatio < rc;
    
    % 
    Snan = S;
    S = NaN( length(f), size(Data,2) );
    S(:,nanInd) = Snan;
    
%     % Additional identification based on peak period:
%     % if the spectal peak is in the wave band (from fc to the Nyquist
%     % frequency), then don't assign it as ice
%     [~,ind] = max( S );
%     fp = f(ind);
%     isIceFP(:,nanInd) = ( (fp<fc) & (fp<fN) ).' ;
%     
%     
%     % Total criteria:
%     isIce = isIceVar & isIceFP;
    

end