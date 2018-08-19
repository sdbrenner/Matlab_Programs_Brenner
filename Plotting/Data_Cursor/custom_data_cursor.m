function custom_data_cursor(fig)
% Modify data cursor to include index values (if available).
% Based on: 
%    https://stackoverflow.com/questions/25477401/
%
% Samuel Brenner, July 2018

if nargin<1
    fig = gcf;
end

dcm_obj = datacursormode(fig);
set(dcm_obj, 'UpdateFcn',@myupdatefcn);


function output_txt = myupdatefcn(~,event_obj)
% Display the position of the data cursor
% ~          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
disp_pos = sprintfc('%4.4g',pos);

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


end%custom_data_cursor