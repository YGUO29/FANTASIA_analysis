% registration modified from suite2p
% by Yueqi Guo 2018/8
% based on code from Spencer Smith lab (by Yiyi Yu)
%%
% folder with exported tiff images
addpath(genpath('D:\=code=\Suite2P'))
folder = 'D:\=data=\80Z_imaging\img_2p\M80Z_180625';
imgfolderlist = dir(fullfile(folder,'*.tif'));
metadata.expdate = folder(end-5:end); % experiment date
metadata.stim = ' '; % visual stimulus
metadata.Regsavepath = 'D:\=data=\80Z_imaging\img_2p\M80Z_180625\registered'; % path for saving registration results
RegPara.useGPU                 = 0; % if you can use a GPU in matlab this accelerate registration approx 3 times
RegPara.NimgFirstRegistration  = 200;
RegPara.PhaseCorrelation       = 1; % set to 0 for non-whitened cross-correlation
RegPara.SubPixel               = Inf; % 2 is alignment by 0.5 pixel, Inf is the exact number from phase correlation
RegPara.NiterPrealign          = 5; % (was 20) number of iterations for the target computation (iterative re-alignment of subset of frames)
RegPara.RawPrecision           = 'uint8';
RegPara.kriging                = 1; % ==== added by Yueqi ====

for i = 1:size(imgfolderlist,1)
        filename = imgfolderlist(i).name
        imglist = imfinfo(fullfile(folder,filename)); % the length is the number of frames
        data1 = imread(fullfile(folder,filename),1); % read the first frame
        metadata.imgsize = size(data1);
        metadata.numframe = length(imglist);
        metadata.name = strtrim(filename);
        % initial alignment (to get target image RegPara.mimg, and other parameters) 
        clear data
        f1 = round(rand(1)*length(imglist)/2); % randomly chose a starting frame, located at the first half of the stack
        for j = f1+1:f1+RegPara.NimgFirstRegistration
            data(:,:,j-f1) = imread(fullfile(folder,filename),j);
        end
        data = single(data);
        RegPara = alignIterative_modified(data, RegPara);           
        clear data
        
        % register whole movie
        RegPara.RegFile = fullfile(metadata.Regsavepath, [metadata.name(1:end-4),'_reg.bin']);
        regdir = fileparts(RegPara.RegFile);
        if ~exist(regdir, 'dir')
            mkdir(regdir);
        end
        
        tic
        ix0 = 0; % index of starting frame
        Nbatch = 500; % how many frames to process in batch
        % init array for registered movie
        dreg = zeros([size(RegPara.mimg), length(imglist)],class(data1)); % data format same as first frame of unregistered data
        
        dsall = []; CorrAll = [];
        while ix0<length(imglist)
            clear mov            
            indxr = ix0 + (1:Nbatch); % take Nbatch frames            
            indxr(indxr>length(imglist)) = [];
            mov = zeros([size(data1), length(indxr)],class(data1));
            for j = indxr(1):indxr(end)
                mov(:,:,j-ix0) =  imread(fullfile(folder,filename),j);
            end
            
            [ds, Corr]  = registration_offsets(mov, RegPara, 0);
            CorrAll = cat(1,CorrAll, Corr);
            dsall = cat(1,dsall, ds);

            dreg(:, :, indxr) = rigidRegFrames(mov, RegPara, ds);
            ix0 = ix0 + Nbatch;
            indxr(end)
        end
        toc
        
        RegPara.dsall = dsall;
        RegPara.CorrAll = CorrAll;
%         RegPara = effectiveRange(RegPara);
        % save result
        clear fid
        fid = fopen(RegPara.RegFile, 'w');
        fwrite(fid, dreg, class(data1));
        fclose(fid);            
        % save reg sample (downsampled temporally 15x)
        sampledreg = dreg(:,:,1:15:end);
        clear mov
        imwrite(sampledreg(:,:,1), fullfile(metadata.Regsavepath, [metadata.name,'_RegSample.tif']))
        for k = 2:size(sampledreg,3)
            imwrite(sampledreg(:,:,k), fullfile(metadata.Regsavepath, [metadata.name(1:end-4),'_RegSample.tif']), 'writemode', 'append');
        end
        clear dreg; clear sampledreg
        save(fullfile(metadata.Regsavepath, [metadata.name(1:end-4),'_regobj.mat']),'metadata','RegPara');
toc
end

%% read bin data
fid = fopen(RegPara.RegFile, 'r');
data = fread(fid, metadata.imgsize(2)*metadata.imgsize(1)*metadata.numframe, '*uint16');
mov = reshape(data, metadata.imgsize(1),metadata.imgsize(2),metadata.numframe);