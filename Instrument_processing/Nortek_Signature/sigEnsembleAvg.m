function Data = sigEnsembleAvg(Data,mode,flds,func)
% SIGENSEMBLEAVG averages data ensembles
%
%   Data = sigEnsembleAvg(Data) peforms ensemble averaging for the Average
%   mode data in the structure 'Data'.
%
%   Data = sigEnsembleAvg(Data,mode) allows specification of the input data
%   mode as 'avg', 'burst', or 'ice' (corresponding to Average, Burst, or
%   AverageIce structure variables).  The function can act on multiple data
%   types by including different modes by including a cell array of modes:
%   e.g. {'avg','ice'}
%
%   Data = sigEnsembleAvg(Data,mode,fields) will only perform ensemble
%   averaging for only the fields specified, where 'fields' is a cell array
%   containing field names. If fields is not input, or is input as an empty
%   value, then averaging will be performed on all fields by default.
%
%   Data = sigEnsembleAvg(...,func) allows for specification of the
%   averaging function.  By default, 'func' is @nanmedian: so ensemble
%   averages are the medions of each ensemble, ignorining NaN values
%   (nanmedian requires the Statistics toolbox).  'func' must be a function
%   definition appropriate for application to vector data, or across rows
%   of matrix data.  While the argument 'func' allows for specification of
%   alternative averageing functions (e.g. @mean), it can be used to
%   perform other operations on ensembles as well (such as @std).  For
%   example, the number of non-NaN values in each ensemble can be found by
%   calling sigEnsembleAvg as follows:
%       nancount = @(x) sum( ~isnan(x) );
%       avgDataCnt = sigEnsembleAvg(averageData,'avg',[],nancount);
%
%   Notes:  
%   (1) This function is developed to operate on Data structures that are
%   output by converting raw .ad2cp data to .mat files using MIDAS
%   software.  Data converted with Signature Deployment software may not
%   have matching variable names.
%   (2) Ensemble identification is done using the variable
%   Data.*_EnsembleCount. If there are no missing samples and no ensembles
%   are cut off at the beginning or end, then the data is considered 'well
%   behaved', and averaging is performed efficiently by reorangizing the
%   data into matrices so all ensembles can be averaged simultaneously.  If
%   data are not well behaved, then ensemble average in done in a loop,
%   where ensembles are delineated by indices where the ensemble count
%   stops increasing.  This is a much slower process.
%
%   S.D.Brenner, 2019

%% Parse inputs

    % Set default values
    if nargin < 4 || isempty(func); func = @nanmean;        end
    if nargin < 3 || isempty(flds); flds = fields(Data);    end
    if nargin < 2 || isempty(mode); mode = 'avg';           end
    
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