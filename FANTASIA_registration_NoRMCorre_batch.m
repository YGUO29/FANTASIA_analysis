% modified by Yueqi 2018/8
clear
gcp;
addpath(genpath('D:\=code=\FANTASIA-NoRMCorre'));
imgfolder = 'D:\=data=\80Z_imaging\img_2p\M80Z_180503';
imgfolderlist = dir(fullfile(imgfolder,'*.tif'));
regfolder = fullfile(imgfolder,'registered_NoRMCorre');
if ~exist(regfolder)
    mkdir(regfolder)
end

%% register within trials

plotON = 0;
for i = 1:size(imgfolderlist,1)
    filename = imgfolderlist(i).name

    imglist = imfinfo(fullfile(imgfolder,filename)); % the length is the number of frames
%     data1 = imread(fullfile(imgfolder,filename),1); % read the first frame

    sframe = 11;    
    warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
    tic; Y = read_file(fullfile(imgfolder,filename),sframe); toc; % read the file (optional, you can also pass the path in the function instead of Y)
    Y = single(Y);                 % convert to single precision 
    T = size(Y,ndims(Y));
    Y = Y - min(Y(:));
    
    % now try non-rigid motion correction (also in parallel)
    % clear large data in memory, to avoid out of memory problems
    options_nonrigid = NoRMCorreSetParms(...
    'd1',size(Y,1),...
    'd2',size(Y,2),...
    'grid_size',[128, 128],... % used to be 128
    'mot_uf',4,...
    'bin_width',200,...
    'max_shift',60,...
    'max_dev',3,...
    'us_fac',50,...
    'init_batch',200,...
    'correct_bidir',0,...
    'shifts_method', 'cubic'); % added 9/2018, to avoid shift limit errors

    % perform motion correction
    tic; [M2,shifts2,template2,options_nonrigid] = normcorre_batch(Y,options_nonrigid); toc
    M2 = uint16(M2);
    
    filename_reg = [regfolder,'\',filename];
    if exist(filename_reg) 
    % popout dialogut window
    choice = questdlg('Tiff file already exist, overwrite or rename?', ...
	'Duplicated file', ...
	'Overwrite','Rename','Cancel','Cancel');
    % Handle response
        switch choice
            case 'Overwrite'
                delete(filename_reg);
                saveastiff(M2,filename_reg); 
            case 'Rename'
                prompt={'Enter an alternative file name'};
                name = 'New name for tiff file';
                defaultans = {filename};
                options.Interpreter = 'tex';
                answer = inputdlg(prompt,name,[1 40],defaultans,options);
                tiffnametemp = [regfolder,'\', answer{1}, '.tif'];
            case 'Cancel'
                return
        end
    else 
        saveastiff(M2,filename_reg);
    end
    
    if plotON
        % compute metrics (non-rigid)
        nnY = quantile(Y(:),0.005);
        mmY = quantile(Y(:),0.995);
        [cY,mY,vY] = motion_metrics(Y,10);
        [cM2,mM2,vM2] = motion_metrics(M2,10);
        T = length(cY);
        % plot metrics

        figure;
        ax1 = subplot(2,4,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')
        ax2 = subplot(2,4,2); imagesc(mM2,[nnY,mmY]); axis equal; axis tight; axis off; title('mean non-rigid corrected','fontsize',14,'fontweight','bold')
        ax3 = subplot(2,4,[3,4]); plot(1:T,cY,'DisplayName','raw data'), hold on
            plot(1:T,cM2,'DisplayName','non-rigid'); 
            legend(ax3); title('correlation coefficients','fontsize',14,'fontweight','bold')
        % plot shift
        % shifts_r = squeeze(cat(3,shifts1(:).shifts));
        shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
        shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
        shifts_x = squeeze(shifts_nr(:,1,:))';
        shifts_y = squeeze(shifts_nr(:,2,:))';
        patch_id = 1:size(shifts_x,2);
        str = strtrim(cellstr(int2str(patch_id.')));
        str = cellfun(@(x) ['patch # ',x],str,'un',0);

        ax4 = subplot(2,4,[5,6]); plot(shifts_x); hold on; title('displacements along x','fontsize',14,'fontweight','bold')
                set(gca,'Xtick',[])
        ax5 = subplot(2,4,[7,8]); plot(shifts_y); hold on; title('displacements along y','fontsize',14,'fontweight','bold')
                xlabel('timestep','fontsize',14,'fontweight','bold')
        linkaxes([ax3, ax4, ax5],'x')

    end
end
%% register across trials using the template from 1st trial 
% Note: set grid_size/max_size > 2. Otherwise the max shift cannot be achieved
clear
gcp;
addpath(genpath('D:\=code=\FANTASIA-NoRMCorre'));
imgfolder = 'D:\=data=\80Z_imaging\img_2p\M80Z_180905';
imgfolderlist = dir(fullfile(imgfolder,'*.tif'));
regfolder = fullfile(imgfolder,'registered_NoRMCorre2');
if ~exist(regfolder)
    mkdir(regfolder)
end

plotON = 0;
NR = 1; % 1: non-rigid, 0: rigid
for i = 1:size(imgfolderlist,1)
    filename = imgfolderlist(i).name
    imglist = imfinfo(fullfile(imgfolder,filename)); % the length is the number of frames
    sframe = 11;    
    warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
    tic; Y = read_file(fullfile(imgfolder,filename),sframe); toc; % read the file (optional, you can also pass the path in the function instead of Y)
    Y = single(Y);                 % convert to single precision 
    T = size(Y,ndims(Y));
    Y = Y - min(Y(:));
    
    if NR == 1   % now try non-rigid motion correction (also in parallel)
        % clear large data in memory, to avoid out of memory problems
        options_nonrigid = NoRMCorreSetParms(...
        'd1',size(Y,1),...
        'd2',size(Y,2),...
        'grid_size',[128, 128],... % used to be 128
        'mot_uf',4,...
        'bin_width',200,...
        'max_shift',60,...
        'max_dev',3,...
        'us_fac',50,...
        'init_batch',200,...
        'correct_bidir',0,...
        'shifts_method', 'cubic'); % added 9/2018, to avoid shift limit errors
        % perform motion correction
        if i == 1
            tic; [M2,shifts2,template,options_nonrigid] = normcorre_batch(Y,options_nonrigid); toc
        else 
            tic; [M2,shifts2,template2,options_nonrigid] = normcorre_batch(Y,options_nonrigid,template); toc
        end
    else
        % rigid motion correction
        options_rigid = NoRMCorreSetParms(...
        'd1',           size(Y,1),...
        'd2',           size(Y,2),...
        'bin_width',    10,...
        'max_shift',    50,...
        'us_fac',       20,...
        'init_batch',   200,...
        'correct_bidir',0);
        % perform motion correction
        if i == 1
            tic; [M2,shifts2,template,options_rigid] = normcorre(Y,options_rigid); toc
        else 
            tic; [M2,shifts2,template2,options_rigid] = normcorre(Y,options_rigid,template); toc
        end
    end
    
    
    
    M2 = uint16(M2);
    filename_reg = [regfolder,'\',filename];
    if exist(filename_reg) 
    % popout dialogut window
    choice = questdlg('Tiff file already exist, overwrite or rename?', ...
	'Duplicated file', ...
	'Overwrite','Rename','Cancel','Cancel');
    % Handle response
        switch choice
            case 'Overwrite'
                delete(filename_reg);
                saveastiff(M2,filename_reg); 
            case 'Rename'
                prompt={'Enter an alternative file name'};
                name = 'New name for tiff file';
                defaultans = {filename};
                options.Interpreter = 'tex';
                answer = inputdlg(prompt,name,[1 40],defaultans,options);
                tiffnametemp = [regfolder,'\', answer{1}, '.tif'];
            case 'Cancel'
                return
        end
    else 
        saveastiff(M2,filename_reg);
    end
    
    if plotON
        % compute metrics (non-rigid)
        nnY = quantile(Y(:),0.005);
        mmY = quantile(Y(:),0.995);
        [cY,mY,vY] = motion_metrics(Y,10);
        [cM2,mM2,vM2] = motion_metrics(M2,10);
        T = length(cY);
        % plot metrics

        figure;
        ax1 = subplot(2,4,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')
        ax2 = subplot(2,4,2); imagesc(mM2,[nnY,mmY]); axis equal; axis tight; axis off; title('mean non-rigid corrected','fontsize',14,'fontweight','bold')
        ax3 = subplot(2,4,[3,4]); plot(1:T,cY,'DisplayName','raw data'), hold on
            plot(1:T,cM2,'DisplayName','non-rigid'); 
            legend(ax3); title('correlation coefficients','fontsize',14,'fontweight','bold')
        % plot shift
        % shifts_r = squeeze(cat(3,shifts1(:).shifts));
        shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
        shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
        shifts_x = squeeze(shifts_nr(:,1,:))';
        shifts_y = squeeze(shifts_nr(:,2,:))';
        patch_id = 1:size(shifts_x,2);
        str = strtrim(cellstr(int2str(patch_id.')));
        str = cellfun(@(x) ['patch # ',x],str,'un',0);

        ax4 = subplot(2,4,[5,6]); plot(shifts_x); hold on; title('displacements along x','fontsize',14,'fontweight','bold')
                set(gca,'Xtick',[])
        ax5 = subplot(2,4,[7,8]); plot(shifts_y); hold on; title('displacements along y','fontsize',14,'fontweight','bold')
                xlabel('timestep','fontsize',14,'fontweight','bold')
        linkaxes([ax3, ax4, ax5],'x')

    end
end


%% set parameters (first try out rigid motion correction)
options_rigid = NoRMCorreSetParms(...
    'd1',           size(Y,1),...
    'd2',           size(Y,2),...
    'bin_width',    10,...
    'max_shift',    25,...
    'us_fac',       20,...
    'init_batch',   200,...
    'correct_bidir',0); % do not perform bi-directional scanning correction (avoid zigzag artifact)

% perform motion correction
tic; [M1,shifts1,template1,options_rigid] = normcorre(Y,options_rigid); toc
% save file as tiff
M1 = uint16(M1);
saveastiff(M1,'D:\=data=\80Z_imaging\img_2p\CaImAn_Results\20180503T121108_nobidir.tif');
% saveastiff(template1,'D:\=data=\80Z_imaging\img_2p\CaImAn_Results\20180503T121108_template.tif');


%% plot a movie with the results

figure;
for t = 1:1:T
    subplot(121);imagesc(Y(:,:,t),[nnY,mmY]); xlabel('raw data','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    subplot(122);imagesc(M2(:,:,t),[nnY,mmY]); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    set(gca,'XTick',[],'YTick',[]);
    drawnow;
    pause(0.02);
end

%% create a stack of average images
imgfolder = 'D:\=data=\80Z_imaging\img_2p\M80Z_180419';
imgfolderlist = dir(fullfile(imgfolder,'*.tif'));

Y_mean = zeros(653,652,length(imgfolderlist));
for i = 1:size(imgfolderlist,1)
    filename = imgfolderlist(i).name
    imglist = imfinfo(fullfile(imgfolder,filename)); % the length is the number of frames
    tic; Y = read_file(fullfile(imgfolder,filename)); toc; % read the file (optional, you can also pass the path in the function instead of Y)
    Y_mean(:,:,i) = mean(Y,3);
end

saveastiff(uint16(Y_mean),fullfile(imgfolder,'mean_images.tif'));


