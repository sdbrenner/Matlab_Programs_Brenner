function bg_patch(patch_start,patch_end,colour,varargin)
    % Function for creating a patch in the background of a figure
    % (e.g. highlighting periods of interest)
    %
    % 'patch_start' and 'patch_end' can be scalar or vectors.  If they
    % are vectors, bg_patch will create multiple patches with start and end
    % locations based on matched pairs of 'patch_starts' and 'patch_ends'.
    %
    % bg_patch(ps,pe); creates light grey background patches that start at
    %   ps and end at pe
    % bg_patch(ps,pe,colour); allows for the specification of the patch
    %   colour.  The 'colour' variable should be a 1x3 vector.
    % bg_patch(ps,pe,colour,'Name','Value'); allows for the input of
    %   additional optional properties to the patch command.  These must
    %   take the form of name,value pairs as they would be passed directly
    %   to the patch command. 
    %   For example, the function can be called as:
    %       bg_patch(ps,pe,colour,'EdgeColor','r');
    %   For details of the choices for Name-Value pair arguments, see
    %   documentation for the patch command.
    %   To specify these parameters without specifying a patch colour,
    %   input an empty value for the colour variable:
    %       bg_patch(ps,pe,[],'Name','Value');
    %
    % Note: bg_patch will set the vertical extent of the patches to the
    %   current y-limits of the axis.  If the y-limits change, bg_patch
    %   objects will not be updated.  Therefore, it is recommended to call
    %   this function after the y-limits of the plot have already been set.
    
    
    %% Initialization (Error checking and setting defaults)
    
    % Number of patch objects
    num_patch = length(patch_start);
    
    % Check that length of patch_start and patch_end vectors match. If they
    % don't, reduce the dimension of the larger of the two vectors.
    if length(patch_end) > num_patch
        patch_end = patch_end(1:num_patch);
        warning( 'patch_start and patch_end should be of the same size. Using decreased number of patch_end vectors.  Behaviour may not be as expected.');
    elseif length(patch_end) < num_patch
        num_patch = length(patch_end);
        patch_start = patch_start(1:num_patch);
        warning( 'patch_start and patch_end should be of the same size. Using decreased number of patch_start vectors.  Behaviour may not be as expected.');
    end
    
    % If no colour is specified, use a light grey:
    if nargin < 3 || any(isempty(colour)) || any(isnan(colour))
        colour =  0.95*[1,1,1];
    end
    
    % I should include error checking for the size of the colour vector and
    % also the optional patch arguements
    
    
    %% Create patches
    
    YL = get(gca,'ylim');
    hold on;
    % loop through number of patches
    for n = 1:num_patch
        patch_x = [patch_start(n),patch_end(n),patch_end(n),patch_start(n)];
        patch_y = YL([1,1,2,2]);
        patch( patch_x,patch_y,colour,varargin{:} );
    end
    hold off;
    set(gca,'children',flipud(get(gca,'children')),...
            'layer','top',...
            'ylim',YL);
end