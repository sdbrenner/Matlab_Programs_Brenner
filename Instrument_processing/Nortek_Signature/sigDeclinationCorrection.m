function Data = sigDeclinationCorrection(Data,mode,declinationOffset,declinationMattime)
% SIGDECLINATIONCORRECTION Declination offset heading correction for Nortek Signature-series
% acoustic doppler current profilers
%
%   Data = sigDeclinationCorrection(Data,declinationOffset)
%   Data = sigDeclinationCorrection(Data,declinationOffset,declinationMattime)


%% Parse inputs

if isempty(mode); mode = 'avg'; end



% Parse mode choice
%   ( Note, 'mode' options could have instead been the 'dataWordChoices'
%     values, but instead are 'modeChoices' to be consistent with other
%     Nortek and Signature codes)
modeChoices = {'avg','ice','burst'};
dataWordChoices = {'Average','AverageIce','Burst'};
[lia,locb] = ismember( lower(mode) , modeChoices );
if ~lia
    error('The input variable ''mode'' must be one of: ''avg'', ''ice'', or ''burst''');
elseif length(lia)>1
    % If multiple mode words are entered, recursively run this script for
    % each of the individually (this may break something)
    for n = 1:length(lia)
        modeN = modeChoices{locb(n)};
        Data = sigDeclinationCorrection(Data,modeN,declinationOffset,declinationMattime);
    end
    return;
else
    dataModeWord = dataWordChoices{locb};
end


%% 

if length(declinationOffset) > 1
    if nargin < 4 || isempty(declinationMattime)
        error('If ''declinationOffset'' is input as a vector, the corresponding variable ''declinationMattime'' must also be supplied');
    end
    declinationOffset = interp1(declinationMattime,declinationOffset,...
                                Data.([dataModeWord,'_MatlabTimeStamp']) );
end



%% Apply heading correction

    heading = Data.([dataModeWord,'_Heading']);
    trueHeading = mod(heading + declinationOffset,360);
    Data.([dataModeWord,'_Heading']) = trueHeading;
    
    % Check if the Data structure includes an AHRS rotation matrix:
    ahrsEnabled = isfield(Data,[dataModeWord  '_AHRSRotationMatrix']);

    % If the machine has an AHRS system, the heading correction needs to be
    % incorporated into the AHRS rotation matrix. I will follow a procedure 
    % similar to what is available in a Nortek-provided code for declination 
    % correction, but replace the declination angle with the correction angle.
    % This needs to be done separately for each timestamp.
    if ahrsEnabled 
        numTimes = length( Data.([dataModeWord,'_Heading']));
        for n = 1:numTimes
            % extract and organize AHRS rotation matrix
            rotMatAHRS = ...
     transpose(reshape(Data.([dataModeWord  '_AHRSRotationMatrix'])(n,:),3,3));
            % calculate correction angle
            corAng = mod( trueHeading(n) - heading(n), 360);
            % construct and appy correction rotation matrix
            corMat = [  cosd(corAng) , sind(corAng) , 0 ;
                       -sind(corAng) , cosd(corAng) , 0 ;
                             0,             0,        1 ];
            corRotMatAHRS = corMat *  rotMatAHRS;            
            % update AHRS rotation matrix
            Data.([dataModeWord  '_AHRSRotationMatrix'])(n,:) = ... 
                                         reshape(transpose(corRotMatAHRS),1,9);
        end
    end
