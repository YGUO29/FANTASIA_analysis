% toolboxes to use in this analysis
clear
addpath(genpath(cd))

filename = 'D:\=data=\80Z_imaging\img_2p\M80Z_180625\20180625T112517_registered.tif';
m = loadtiff(filename);
% filename = 'D:\=data=\80Z_imaging\img_2p\M80Z_180625\AVG_20180625T112517_registered.tif';
% m_avg = loadtiff(filename);


%% code suggestions from Spencer Smith
% It's just a way to generate some contrast. 
% The blurring step is useful for making the ROIs look nicer, but I only use it to generate ROIs.
% The data from the traces is extracted from the unblurred video.

% m = the_calcium_imaging_video; %// load the video somehow
G = fspecial('gaussian',[30 30],2); %// make a gaussian filter
%%
for f=1:size(m,3) %// for each frame in the video
    mF1(:,:,f) = imfilter(m(:,:,f),G,'same'); %// gaussian filter the video (it's inefficient to use a loop here, but it works for trying this out)
end
mKurt1 = kurtosis(mF1,0,3); %// compute the kurtosis of each pixel's time course
%% 
a = [10 10]; 
b = -0.03;
bw=adaptivethreshold(mKurt1,a,b,0); % numbers determined empirically, a = [100 100], b=-0.1 worked for some data;

figure,imshow(bw,[]),title(['a = [ ',num2str(a),' ], b = ',num2str(b)])

regions=bwlabeln(bw);
s=regionprops(regions);

%%
figure,imshow(ccimage,[])