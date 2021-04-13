function [E,f] = jimSpectra(z,winSize,overlap,nfft,fs)



if ~ismatrix(z) 
    error('z must be vector or matrix');
elseif min(size(z)) == 1
	[E,f] = jimSpectra1D(z,winSize,overlap,nfft,fs);
elseif min(size(z)) > 1
    for n = 1:size(z,2)
        zn = z(:,n);
        [En,f] = jimSpectra1D(zn,winSize,overlap,nfft,fs);
        E(:,n) = En;
    end
else
    error('something has gone wrong');
end



end


function [E,f] = jimSpectra1D(z,winSize,overlap,nfft,fs)


    %% Default inputs

    wsecs = 256;    % window length in seconds, should make 2^N samples
%     wsecs = 512;
    merge = 1;      % freq bands to merge, must be odd?
    maxf =  2*fs;   % Nyquist frequency

    % detrend first
    z = detrend(z);
    pts = length(z);    

    % break into windows (use 75 percent overlap)
    w = round(fs * wsecs);  % window length in data points
    if rem(w,2)~=0, w = w-1; else end  % make w an even number
    windows = floor( 4*(pts/w - 1)+1 );   % number of windows, the 4 comes from a 75% overlap
    dof = 2*windows*merge; % degrees of freedom
    % loop to create a matrix of time series, where COLUMN = WINDOW 
    zwindow = zeros(w,windows);
    for q=1:windows
        zwindow(:,q) = z(  (q-1)*(.25*w)+1  :  (q-1)*(.25*w)+w  );  
    end

    % detrend individual windows (full series already detrended)
    for q=1:windows
    zwindow(:,q) = detrend(zwindow(:,q));
    end

    % taper and rescale (to preserve variance)
    % form taper matrix (columns of taper coef)
    taper = sin ( (1:w) * pi/w )' * ones(1,windows); 
    % taper each window
    zwindowtaper = zwindow .* taper;
    % now find the correction factor (comparing old/new variance)
    factz = sqrt( var(zwindow) ./ var(zwindowtaper) );
    % and correct for the change in variance
    % (mult each window by it's variance ratio factor)
    zwindowready = (ones(w,1)*factz).* zwindowtaper;


    % FFT
    % calculate Fourier coefs
    Zwindow = fft(zwindowready);
    % second half of fft is redundant, so throw it out
    Zwindow( (w/2+1):w, : ) = [];
    % throw out the mean (first coef) and add a zero (to make it the right length)  
    Zwindow(1,:)=[]; 
    Zwindow(w/2,:)=0; 
    % POWER SPECTRA (auto-spectra)
    ZZwindow = real ( Zwindow .* conj(Zwindow) );



    % merge neighboring freq bands (number of bands to merge is a fixed parameter)
    % initialize
    ZZwindowmerged = zeros(floor(w/(2*merge)),windows);

    for mi = merge:merge:(w/2) 
        ZZwindowmerged(mi/merge,:) = mean( ZZwindow((mi-merge+1):mi , : ) );
    end
    % freq range and bandwidth
    n = (w/2) / merge;                         % number of f bands
    Nyquist = .5 * fs;                % highest spectral frequency 
    bandwidth = Nyquist/n ;                    % freq (Hz) bandwitdh
    % find middle of each freq band, ONLY WORKS WHEN MERGING ODD NUMBER OF BANDS!
    f = 1/(wsecs) + bandwidth/2 + bandwidth.*(0:(n-1)) ; 


    % ensemble average windows together
    % take the average of all windows at each freq-band
    % and divide by N*samplerate to get power spectral density
    % the two is b/c Matlab's fft output is the symmetric FFT, and we did not use the redundant half (so need to multiply the psd by 2)
    ZZ = mean( ZZwindowmerged.' ) / (w/2 * fs  );

    % auto and cross displacement spectra, assuming deep water
    E = ZZ;  %[m^2/Hz]

    
end