function Data = sigBeamMapping(Data,Config,mode,output)
% SIGBEAMMAP Maps beam coordinates to xyz, enu for Nortek Signature500
% ADCPs (including those with AHRS systems)
%   Data = sigBeamMap(Data,Config) 
%
%   Data = sigBeamMap(Data,Config,mode) 
%   Data = sigBeamMap(Data,Config,mode,output) 
%
%   Note: Nortek can provide a Matlab script that performs the same
%   functionality, and is more general.  This script is a adapted from
%   theirs, and slightly simplified/adjusted to suit personal processing
%   needs and practices and provide speed improvements.
%
%   S.D.Brenner, 2019


%% Parse inputs

    % Assign default behaviours:
    if nargin < 3 || isempty(output);   output = {'beam','xyz','enu'};  end
    if nargin < 2 || isempty(mode);     mode = {'avg'};                 end

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


    % Check 'output' choice(s)
    output = lower(output);
    outputChoices = {'beam','xyz','enu'};
    if ~ismember( output , outputChoices )
        error('The input variable ''output'' must be one of ''beam'', ''xyz'', or ''enu''');
    end

    % Assign and check number of beams (assuming that beam mapping is the 
    % same for all measurements in a data file)
    numBeams = length( find(Data.( [dataModeWord '_Physicalbeam'] )( 1, : ) > 0) );
    if numBeams ~=4
        warning('This function developed for 4-beam configurations and may produce undesirable results (or may break)');
    end

%% Beam to XYZ transformation:

    % Get transformation matrix:
    T_beam2xyzz = reshape(Config.([configModeWord,'_beam2xyz']),[4,4]).';
    
    % Velocity is essentially 3D: depth (bin#), time, beam; But the 
    % transformation matrix can only be applied in 2D (matrix multiplication).  
    % I'll map depth and time onto a single vector to create a 2D matrix:
    [numSamples,numBins] = size( Data.([dataModeWord,'_VelBeam1']) );
    velBeam = NaN(numBeams, (numSamples*numBins) );
    for n = 1:numBeams
        velBeam(n,:) =...
               reshape(Data.([dataModeWord,'_VelBeam',num2str(n)]), 1, []);
    end

    % Apply transformation        
    velXYZZ = T_beam2xyzz * velBeam ;

    % Update Data structure if output asks for 'xyz'
    if any(ismember(output,'xyz'))
        Data.([dataModeWord,'_VelX']) = reshape( velXYZZ(1,:), numSamples, numBins);
        Data.([dataModeWord,'_VelY']) = reshape( velXYZZ(2,:), numSamples, numBins);
        Data.([dataModeWord,'_VelZ1']) = reshape( velXYZZ(3,:), numSamples, numBins);
        Data.([dataModeWord,'_VelZ2']) = reshape( velXYZZ(4,:), numSamples, numBins);
    end

%% XYZ to ENU transformation
% The method to transform from XYZ to ENU coordinates depends on whether
% AHRS is installed in the system.  First, check if the Data structure
% contains an AHRS rotation matrix, then apply the appropriate
% transformation.

    if any(ismember(output,'enu'))

        % Check if the Data structure includes an AHRS rotation matrix:
        ahrsEnabled = isfield( Data, [dataModeWord,'_AHRSRotationMatrix'] );

        % Re-organize velXYZZ to 3D data structure in dim-time-depth:
        velXYZZ3D = reshape(velXYZZ,[numBeams,numSamples,numBins]);

        % Pre-allocate enuu velocity matrix
        velENUU3D = NaN( size(velXYZZ3D) );

        % XYZ to ENU data have separate transformation matrices for each point in
        % time to to time-variable headings.  To apply the tranformation, loop
        % through timestamps
        for n = 1:numSamples
            % Generate transformation matrix for sample n (depends on whether the
            % instrument has AHRS or not)
            if ahrsEnabled
                T_xyzz2enuu_n = makeTransMatAHRS(Data,dataModeWord,n);
            else
                T_xyzz2enuu_n = makeTransMat(Data,dataModeWord,n);
            end
            % Extract XYZZ velocity data for sample n
            velXYZZ_n = squeeze( velXYZZ3D(:,n,:) );
            % Apply transformation matrix
            velENUU_n = T_xyzz2enuu_n * velXYZZ_n;
            % Update 3D velocity matrix
            velENUU3D(:,n,:) = permute(velENUU_n,[1,3,2]);
        end


        % Update Data structure
        Data.([dataModeWord,'_VelEast']) =  squeeze( velENUU3D(1,:,:) );
        Data.([dataModeWord,'_VelNorth']) = squeeze( velENUU3D(2,:,:) );
        Data.([dataModeWord,'_VelUp1']) =   squeeze( velENUU3D(3,:,:) );
        Data.([dataModeWord,'_VelUp2']) =   squeeze( velENUU3D(4,:,:) );


    end

    %% Check if beam velocities are a requested output
    % Otherwise, remove them from the data structure

    if ~any( ismember( output,'beam') )
        Data = rmfield(Data, { [ dataModeWord '_VelBeam1' ],...
                               [ dataModeWord '_VelBeam2' ],...
                               [ dataModeWord '_VelBeam3' ],...
                               [ dataModeWord '_VelBeam4' ] } );
    end


end





%% SUPERCEDED
% 
% % Transformation of spherical to cartesian coordinates
% 
% for n = 1:numBeams
%     theta(n,1) = Config.(['beamConfiguration',num2str(n),'_theta']);
%     phi(n,1) = Config.(['beamConfiguration',num2str(n),'_phi']);
% end
% 
% 
% mat = [ 1, -1, 1, 0 ;...
%   1, -1, 0, 1 ;...
%   1, -1, 1, 0 ;...
%   1, -1, 0, 1 ];
% 
% sphere2cart = mat.*[ sind(theta).*cosd(phi),...
%                 sind(theta).*sind(phi),...
%                cosd(theta),...
%                 cosd(theta) ];
% T_beam2xyz = inv( sphere2cart ); 



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
    T_xyzz2enuu = H*PR; % matrix multiplication
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