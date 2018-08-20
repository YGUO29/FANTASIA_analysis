function movF = MovGaussFilter(mov, G, subsample)
% filtered
% G = fspecial('gaussian',[4 4],2);
if subsample == 0
    movF = zeros(size(mov), 'uint8');
    for f=1:size(mov,3)
        movF(:,:,f) = uint8(imfilter(mov(:,:,f),G,'same'));
    end 
else 
% filtered subtime
    movF = zeros([size(mov, 1), size(mov, 2), floor(size(mov, 3)/subsample)], 'uint8');
    for f=1:floor(size(mov,3)/subsample)
        movF(:,:,f) = uint8(imfilter(mov(:,:,f*subsample),G,'same'));
    end 
end
