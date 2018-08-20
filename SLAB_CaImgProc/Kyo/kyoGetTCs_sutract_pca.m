function [tcraw, tc, tcraw_pca, tc_pcasub] = kyoGetTCs_sutract_pca(mov, rois, subtractPara)

r1 = subtractPara(1);
r2 = subtractPara(2); 

[d1,d2,T] = size(mov);

tic
lowSignal=0;
rg=1; roiall = sum(rois,3);
tcraw = zeros(size(mov,3),size(rois,3));
tc = zeros(size(mov,3),size(rois,3));
tc_pcasub = zeros(size(mov,3),size(rois,3));
tcraw_pca = zeros(size(mov,3),size(rois,3));

for r=1:size(rois,3);
    roi=rois(:,:,r);
    if r1>0
        dila1 = bwmorph(roi,'dilate', r1);
    else 
        dila1 = roi;
    end
    dila2 = bwmorph(roi,'dilate', r2);
    subbg1 = dila2 - dila1; 
    subbg = and(subbg1, ~roiall);
    
    [ii, jj] = find(subbg==1);
    [ii1, jj1] = find(roi==1);
    siglist = []; bglist = [];
    for i = 1:length(ii)
        bglist(:,i) = squeeze(mov(ii(i),jj(i),:));
    end
    for i = 1:length(ii1)
        siglist(:,i) = squeeze(mov(ii1(i),jj1(i),:));
    end
%     [U,S,~,~,~,mu] = pca(double(bglist));
    [U,S,~,~,~,mu] = pca(double(bglist),'Algorithm', 'SVD','NumComponents', 1);
    
    recontrace = S(:,1)*U(:,1)';
    df = mean(siglist,2)- mean(recontrace,2);
    df1 = mean(siglist,2);
    
    tcraw(:,rg) = df1;
    f = quantile(df1,0.4); % This seems to work nicely for most signals, but it may need to be changed for highly active cells
        if f==0 % in case of very weak signal, use the mean instead
            lowSignal=lowSignal+1;
            tc(:,rg) = (df1./mean(df1))-mean(df1);
        else
            tc(:,rg) = (df1./f)-1;
        end
        
        tcraw_pca(:, rg)=df;
        f = quantile(df,0.4); % This seems to work nicely for most signals, but it may need to be changed for highly active cells
        if f==0 % in case of very weak signal, use the mean instead
            lowSignal=lowSignal+1;
            tc_pcasub(:,rg) = (df./mean(df))-mean(df);
        else
            tc_pcasub(:,rg) = (df./f)-1;
        end
        rg = rg+1;
end
toc
if lowSignal>0
    fprintf('%i ROI(s) had weak signals.',lowSignal);
end
