function make_gif(filename,im,delaytime)
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
disp('complete');

end