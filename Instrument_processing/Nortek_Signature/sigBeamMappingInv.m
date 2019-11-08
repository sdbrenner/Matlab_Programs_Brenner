function Data = sigBeamMappingInv(Data,Config,mode)
% SIGBEAMMAPINV performs inverse mapping of enu coordinates to beam
% coordinates for Nortek Signature500 ADCPs (including those with AHRS
% systems)
%
%   Data = sigBeamMappingInv(Data,Config,mode)
%
%   S.D.Brenner, 2019

%% Parse inputs

    % Assign default behaviour:
    if nargin < 2 || isempty(mode);     mode = {'avg'}; end

    % Parse 'mode' choice(s)
    mode = lower(mode);
    modeChoices = {'avg','ice','burst'};
    dataWordChoices = {'Average','AverageIce','Burst'};
    configWordChoices = {'avg','avg','burst'};
    [modeLog,modeInd] = ismember( mode , modeChoices );
    if ~modeLog
        error('The input variable ''mode'' must be one of: ''avg'', ''ice'', or ''burst''');
    elseif length(modeLog) > 1
        % If multiple mode words are entered, recursively run this script for
        % each of them individually (this may break something)
        for n = 1:length(lia)
            modeN = modeChoices{locb(n)};
            Data = sigBeamMap(Data,Config,modeN,output);
        end
        return;
    else
        ad2cpstr = '';
        dataModeWord = dataWordChoices{modeInd};
        configModeWord = configWordChoices{modeInd};
    end
    
    
    % Assign and check number of beams (assuming that beam mapping is the 
    % same for all measurements in a data file)
    numBeams = length( find(Data.( [dataModeWord '_Physicalbeam'] )( 1, : ) > 0) );
    if numBeams ~=4
        warning('This function developed for 4-beam configurations and may produce undesirable results (or may break)');
    end
    
    
%% ENU to XYZ transformation

    % Extract data from structure (into 3D array)
    [~, numBins] = size( Data.([dataModeWord,'_VelEast']) );
    if numBins > 1   
        velENUU3D(1,:,:) = permute( Data.([dataModeWord,'_VelEast']),  [3,1,2]);
        velENUU3D(2,:,:) = permute( Data.([dataModeWord,'_VelNorth']), [3,1,2]);
        velENUU3D(3,:,:) = permute( Data.([dataModeWord,'_VelUp1']),   [3,1,2]);
        velENUU3D(4,:,:) = permute( Data.([dataModeWord,'_VelUp2']),   [3,1,2]);
    else
        velENUU3D(1,:) = permute( Data.([dataModeWord,'_VelEast']),  [3,1,2]);
        velENUU3D(2,:) = permute( Data.([dataModeWord,'_VelNorth']), [3,1,2]);
        velENUU3D(3,:) = permute( Data.([dataModeWord,'_VelUp1']),   [3,1,2]);
        velENUU3D(4,:) = permute( Data.([dataModeWord,'_VelUp2']),   [3,1,2]);        
    end
        
        
    % Check if the Data structure includes an AHRS rotation matrix:
    ahrsEnabled = isfield( Data, [dataModeWord,'_AHRSRotationMatrix'] );

    % Pre-allocate XYZZ velocity array
    velXYZZ3D = NaN( size(velENUU3D) );
    
    % Loop through timestamps
    [numSamples,numBins] = size( Data.([dataModeWord,'_VelEast']) );
    for n = 1:numSamples
        % Generate transformation matrix for sample n (depends on whether the
        % instrument has AHRS or not)
        if ahrsEnabled
            T_xyzz2enuu_n = makeTransMatAHRS(Data,dataModeWord,n);
        else
            T_xyzz2enuu_n = makeTransMat(Data,dataModeWord,n);
        end    
        
        % Extract ENUU velocity data for sample n
        velENUU_n = squeeze( velENUU3D(:,n,:) );
        
        % Apply transformation matrix
        velXYZZ_n = T_xyzz2enuu_n\velENUU_n; 

        % Update 3D velocity matrix
        velXYZZ3D(:,n,:) = permute(velXYZZ_n,[1,3,2]);
    end        

    % Update data structure
    Data.([dataModeWord,'_VelX'])  = squeeze( velXYZZ3D(1,:,:) );
    Data.([dataModeWord,'_VelY'])  = squeeze( velXYZZ3D(2,:,:) );
    Data.([dataModeWord,'_VelZ1']) = squeeze( velXYZZ3D(3,:,:) );
    Data.([dataModeWord,'_VelZ2']) = squeeze( velXYZZ3D(4,:,:) );
    

%% XYZZ to BEAM transformation

    % Get transformation matrix:
    T_beam2xyzz = reshape(Config.([configModeWord,'_beam2xyz']),[4,4]).';
    
    velBeam3D = NaN( size(velENUU3D) );
    
    % Loop through timestamps
    for n = 1:numSamples
        % Extract XYZZ velocity data for sample n
        velXYZZ_n = squeeze( velXYZZ3D(:,n,:) );
        
        % Apply transformation matrix
        velBeam_n = T_beam2xyzz\velXYZZ_n;
        
        % Update 3D velocity matrix
        velBeam3D(:,n,:) = permute(velBeam_n,[1,3,2]);
    end

    % Update data structure
    Data.([dataModeWord,'_VelBeam1']) = squeeze( velBeam3D(1,:,:) );
    Data.([dataModeWord,'_VelBeam2']) = squeeze( velBeam3D(2,:,:) );
    Data.([dataModeWord,'_VelBeam3']) = squeeze( velBeam3D(3,:,:) );
    Data.([dataModeWord,'_VelBeam4']) = squeeze( velBeam3D(4,:,:) );
    
    
    
    
%% Superceded %% XYZZ to BEAM transformation    
% 
%     % Re-map from 3D beam-sample-bin to 2D beam-[samplebin] (i.e., re-map
%     % depth and time to a single vector for each beam)
%     velXYZZ = reshape( velXYZZ3D, numBeams,[] );
%     for n = 1:numBeams
%         velXYZZ(n,:) = reshape( velXYZZ3D(n,:,:), 1, []);
%     end
%     
%     % ( note: it should be possible to vectorize this, but I was having
%     % issues )
% %     velXYZZ = reshape( permute(velXYZZ3D,[]), numBeams,[] );
% 
%     % Get transformation matrix:
%     T_beam2xyzz = reshape(Config.([configModeWord,'_beam2xyz']),[4,4]).';
% 
%     % Apply transformation   
%     velBeam = T_beam2xyzz\velXYZZ;
%     
%     % Update Data structure
%     for n = 1:numBeams
%         Data.([dataModeWord,'_VelBeam',num2str(n)]) =...
%                              reshape( velXYZZ(n,:), [numSamples,numBins] );
%     end

end


%% EMBEDDED FUNCTIONS %% ==================================================

function T_xyzz2enuu = makeTransMatAHRS(Data,dataModeWord,n)
    % Get raw transformation matrix for sample n
    T_xyzz2enuu = transpose(reshape(Data.([dataModeWord  '_AHRSRotationMatrix'])(n,:),3,3));
    % Modify transfromation matrix to account for "twoZs"
    T_xyzz2enuu = modTransMat(T_xyzz2enuu);
end


function T_xyzz2enuu = makeTransMat(Data,dataModeWord,n)
    heading = Data.( [dataModeWord '_Heading'] )(n);
    pitch = Data.( [dataModeWord '_Pitch'] )(n);
    roll = Data.( [dataModeWord '_Roll'] )(n);
    
    % X-axis points to North; rotate so East is the first component:
    heading = heading -90;
    
    % Make heading matrix
    ch = cosd(heading);
    sh = sind(heading);
    H = [ ch, sh, 0;
         -sh, ch, 0;
           0,  0, 1 ];
    % Make tilt matrix
    cp = cosd(pitch);
    sp = sind(pitch);
    cr = cosd(roll); 
    sr = sind(roll);
    PR = [ cp, -sp*sr, -sp*cr; 
            0,     cr,    -sr;
           sp,  cp*sr,   cp*cr ];
    % Make transformation matrix
    T_xyzz2enuu = H * PR; % matrix multiplication
    % Modify transfromation matrix to account for "twoZs"
    T_xyzz2enuu = modTransMat(T_xyzz2enuu);
end

function T_xyzz2enuu = modTransMat(T_xyzz2enuu)
    % Modify transfromation matrix to account for "twoZs"
    T_xyzz2enuu(1,3) = T_xyzz2enuu(1,3)/2;
    T_xyzz2enuu(1,4) = T_xyzz2enuu(1,3);
    T_xyzz2enuu(2,3) = T_xyzz2enuu(2,3)/2;
    T_xyzz2enuu(2,4) = T_xyzz2enuu(2,3);
    
    T_xyzz2enuu(4,:) = T_xyzz2enuu(3,:);
    T_xyzz2enuu(3,4) = 0;
    T_xyzz2enuu(4,4) = T_xyzz2enuu(3,3);
    T_xyzz2enuu(4,3) = 0;
end