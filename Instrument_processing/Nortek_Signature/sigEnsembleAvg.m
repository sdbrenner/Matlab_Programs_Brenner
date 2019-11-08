function Data = sigEnsembleAvg(Data,mode,flds,func)
% SIGENSEMBLEAVG averages ensembles
%
%   avgData = sigEnsembleAvg(Data,mode)
%
%   avgData = sigEnsembleAvg(Data,mode,fields)
%   avgData = sigEnsembleAvg(Data,mode,fields,func)

%% Parse inputs

    % Set default values
    if nargin < 4 || isempty(func); func = @nanmean; end
    if nargin < 3 || isempty(flds); flds = fields(Data); end
    if nargin < 2 || isempty(mode); mode = 'avg'; end
    
    % Parse mode choice
    %   ( Note, 'mode' options could have instead been the 'dataWordChoices'
    %     values, but instead are 'modeChoices' to be consistent with other
    %     Nortek and Signature codes)
    modeChoices = {'avg','ice','burst'};
    dataWordChoices = {'Average','AverageIce','Burst'};
    [modeLog,modeInd] = ismember( lower(mode) , modeChoices );
    if ~modeLog
        error('The input variable ''mode'' must be one of: ''avg'', ''ice'', or ''burst''');
    elseif length(modeLog)>1
        % If multiple mode words are entered, recursively run this script for
        % each of the individually (this may break something)
        for n = 1:length(modeLog)
            modeN = modeChoices{modeInd(n)};
            Data = sigEnsembleAvg(Data,modeN,flds,func);
        end
        return;
    else
        dataModeWord = dataWordChoices{modeInd};
    end
    
%% Check that ensembles and data sizes are well behaved

    % Flag to indicate if ensembles and data sizes are good:
    % (assume yes initially)
    ensCheck = 1;
    
    % Extract EnsembleCount values
    ensembleCount = Data.([dataModeWord,'_EnsembleCount']);
    % Samples per ensemble:
    sampPerEns = max( ensembleCount );
    % Check size and behaviour:
    try
        % Check that total length of matrix is an even number of ensembles:
        % (if not, the reshape command will fail)
        ensembleMatrix = reshape(ensembleCount,sampPerEns,[]).';
        if ~isStatic( ensembleMatrix ) || ...    % if the matrix isn't static
            any( diff(ensembleMatrix(1,:)) < 0 ) % if rows are not monotonic
            ensCheck = 0;
        end
    catch
        ensCheck = 0;
    end
    
%% Perform averaging
% If the ensembles and data sizes are well behaved (ensCheck=1), then
% average each field using reshape and matrix averages.  If not
% (enseCheck=0), loop through ensembles (still requires a certain degree of
% "good behaviour", but is more robust to missing measurements).

    numFlds = length(flds);

    if ensCheck

        % Loop through fields
        for n = 1:numFlds
            fldName = flds{n};
            [numSamples,numBins] = size( Data.(fldName) );
            if numSamples > 1 % only average non-"static" fields
                var = reshape( Data.(fldName), sampPerEns, [] );
                if strcmpi( fldName, [dataModeWord,'_MatlabTimeStamp'] )
                    avgData.(fldName) = var(1,:).';
                else
                avgData.(fldName) = reshape( func(var), [], numBins);
                end
            else
                % If the variable is static (e.g. 'Ranges'), don't average,
                % but add it to the averaged structure (once only though)
                avgData.(fldName) = Data.(fldName);
            end
        end


    else
        warning('Ensembles are not "well behaved"; performing average in loop');

        % Find indices when the ensemble count stops increasing:
        deltaCount = diff(double( Data.([dataModeWord,'_EnsembleCount']) ));
        negInd = find( deltaCount<0 );
        negInd = [ 0; negInd; length( Data.([dataModeWord,'_EnsembleCount']) ) ];
        numEnsembles = length(negInd)-1; 

        % Initialize (pre-allocate) fields
        flds = fields(Data);
        numFlds = length(flds);
        for n = 1:numFlds
            fldName = flds{n};
            [numSamples,numBins] = size( Data.(fldName) );
            if numSamples > 1
                avgData.(fldName) = NaN( numEnsembles, numBins );
            else
                avgData.(fldName) = NaN( 1, numBins);
            end
        end

        % Loop through ensembles
        for l = 1:numEnsembles
            % Get indices to average
            avgInd =   (negInd(l)+1) : negInd(l+1) ;
            % Loop through fields and perform average
            for n = 1:numFlds
                fldName = flds{n};
                [numSamples,~] = size( Data.(fldName) );
                if numSamples > 1 % only average non-"static" fields
                    if strcmpi( fldName, [dataModeWord,'_MatlabTimeStamp'] )
                        avgData.(fldName)(l,:) = Data.(fldName)(avgInd(1),:);
                    else
                        avgData.(fldName)(l,:) = func( Data.(fldName)(avgInd,:) );
                    end
                elseif numSamples == 1 && l==1
                    % If the variable is static (e.g. 'Ranges'), don't average,
                    % but add it to the averaged structure (once only though)
                    avgData.(fldName) = Data.(fldName);
                end
            end
        end
    end

    % Replace Data structure with new, averaged structure
    Data = avgData;
    
end   

%% EMBEDDED FUNCTIONS %% ==================================================

function static = isStatic(var)
    % check if matrix columns are static (ignoring nans)
    ind = ~isnan( var(:,1) );
    static = all( var(ind,:) == var(1,:) ,'all');
end    