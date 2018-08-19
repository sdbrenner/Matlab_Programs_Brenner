function data_cursor_datest(ax)
% Modify data cursor to use date format.
% Based on:
%   https://www.mathworks.com/matlabcentral/answers/117068
% and adjusted to use any given input axis.
% If not input is given in the function call it will apply the date format
% to the x-axis label by default.
%
% Samuel Brenner, July 2018

if nargin == 0
    axn = 1;
else
    axn = find(strcmp({'x','y','z'},ax));
end

dcm_obj = datacursormode(gcf);
set(dcm_obj, 'UpdateFcn',{@myupdatefcn,axn});

function output_txt = myupdatefcn(~,event_obj,axn)
% Display the position of the data cursor
% ~          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
disp_pos = sprintfc('%4.4g',pos);
disp_pos{axn} = datestr(pos(axn));
output_txt = {['X: ', disp_pos{1}],... 
              ['Y: ', disp_pos{2}] };        
          
% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',disp_pos{3}];
end%if

% If a DataIndex is available, show that also
info = getCursorInfo(dcm_obj);
if isfield(info,'DataIndex')
    DataIndex = [info.DataIndex];
    output_txt{end+1} = sprintf('Index: %d\n', DataIndex(1));
end%if


end%myfunction


end%data_cursor_datest