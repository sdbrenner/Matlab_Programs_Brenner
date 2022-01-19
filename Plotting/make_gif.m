function make_gif(filename,im,delaytime)
% MAKE_GIF is used to generate a .gif file from a set of plot images
%
%   make_gif(filename,im) generates a .gif file with a given filename from
%   a MATLAB cell array, 'im'. The input im is expected to be generated
%   using getframe in a for-loop that is updating a plot.
%   For example:
%
%
%   make_gif(filename,im,delaytime) allows specification of the delay time
%   between frames in the gif. If it is left blank, then a delay time of
%   1/8 s is used.

if nargin == 2
    delaytime = 1/8;
end

for gif_idx = 1:length(im)
    [A,map] = rgb2ind(im{gif_idx}.cdata,256);
    if gif_idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',delaytime);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',delaytime);
    end
end
disp('gif saved');

end