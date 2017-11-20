% modified from make dF/F 030716, to correct the bug from spike2 recordings with the initial protocol 05/07/16 

clear; clc;


% cd 'E:\Lab\Data\wholeBrain\fMRI\Data_from_Eve\2017-07-27'

addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/piotr_toolbox'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/wholeBrainDX'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/sigTOOL'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/CalciumDX'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/bfmatlab'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/chatAnalysis'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/piotr_toolbox'))



filelist = readtext('files_pre.txt', ' ');
fnms = filelist(:, 1);
blueInitial = filelist(:, 2);
mask_fnms = filelist(:, 3);
angle = filelist(:, 4);
no_movies = length(fnms);
downSampleRatio = 0.5;
fixUV = 1;
hat = 300; % window size for rolling hat algorithm

rotateAngle = angle{1};
% pxlCoords = [200, 60; 200, 400; 110, 160; 110, 300; 260, 120; 260, 340; ...
%     350, 180; 350, 260; 390, 120; 390, 330];

% for n = 1
for n = 1:no_movies
    
    clear mask_id filtered1 filtered2 A_sliced B_sliced A_filtered B_filtered diff 
    
    filename = fnms{n};
    concatList = dir(fullfile([filename(1 : end-4), '*.tif']));
    
    A = [];
    B = [];
    for c = 1 : length(concatList)
        imgall = openMovie(concatList(c).name);
        szall = size(imgall);
        if mod(c, 2) == 1
            blueFrames = blueInitial{n} : 2: szall(3);
            uvFrames = setdiff(1:szall(3), blueFrames);
        else
            uvFrames = blueInitial{n} : 2: szall(3);
            blueFrames = setdiff(1:szall(3), uvFrames);
        end
            
        A = cat(3, A, imresize(imgall(:, :, blueFrames), downSampleRatio, 'bilinear'));
        B = cat(3, B, imresize(imgall(:, :, uvFrames), downSampleRatio, 'bilinear'));
        wrongId = 200:200:size(B, 3);
        for id = wrongId(1:end-1)
            A(:, :, id) = (A(:, :, id-1) + A(:, :, id+1))/2;
            B(:, :, id) = (B(:, :, id-1) + B(:, :, id+1))/2;
        end
        
        clear imgall
    end
    
    A = imrotate(A, rotateAngle);
    B = imrotate(B, rotateAngle);


    if size(A, 3) >= size(B, 3)
        A = A(:, :, 1:size(B, 3));
    else
        B = B(:, :, 1:size(A, 3));
    end
    A = A(:, :, 2:end-1);      
    B = B(:, :, 2:end-1);



    sz = size(A); szZ=sz(3);
    npix = prod(sz(1:2));
    A = reshape(A, npix, szZ); %reshape 3D array into space-time matrix
    B = reshape(B, npix, szZ);


%     if fixUV == 1
%         sumUV = sum(B);
%         dI = sumUV(2:end) - sumUV(1:end-1);
%         shortFrameId = find(dI < mean(dI) - 3*std(dI)) + 1;
%         if ~isempty(shortFrameId)
%             if shortFrameId(end) == szZ 
%                 tmpId = shortFrameId(1 : end-1);
%                 B(:, tmpId) = (B(:, tmpId - 1) + B(:, tmpId + 1)) / 2;
%                 B(:, shortFrameId(end)) = B(:, shortFrameId(end - 1));
%             else
%                 B(:, shortFrameId) = (B(:, shortFrameId - 1) + B(:, shortFrameId + 1)) / 2;
%             end
%         end
%     end



    ROI = ReadImageJROI(mask_fnms{n});

    mask = poly2mask(ROI.mnCoordinates(:, 1), ROI.mnCoordinates(:, 2), ...
        szall(1), szall(2));
    mask = imresize(mask, downSampleRatio, 'bilinear');
    mask = imrotate(mask, rotateAngle);

    [mask_r, mask_c] = find(mask > 0);
    mask_id = find(mask > 0);


    % remove slow drifting trend
        hat = 300;
        se = strel('line', hat, 0);
        A_imgDiff = A(mask_id, 1) - A(mask_id, 300);
        B_imgDiff = B(mask_id, 1) - B(mask_id, 300);
        aa = flipud(A(mask_id, 1:300)')';
        bb = flipud(A(mask_id, 1:300)')';
        A_sliced = [-aa + 2 * A(mask_id, 1), A(mask_id, :)];
        B_sliced = [-bb + 2 * B(mask_id, 1), B(mask_id, :)];
%         A_sliced = A(mask_id, :);
%         B_sliced = B(mask_id, :);
        
        parfor p = 1:length(mask_id)   
        % parfor p = 1:1000 
            filtered1(p, :) = imtophat(A_sliced(p, :), se);
%             filtered2(p, :) = imtophat(B_sliced(p, :), se);
        end
%         A_filtered = zeros(size(A, 1), size(A, 2) + 300);
%         B_filtered = zeros(size(B, 1), size(B, 2) + 300);
        A_filtered = zeros(size(A));
%         B_filtered = zeros(size(B));
        A_filtered(mask_id, :) = filtered1(:, 301:end);
%         B_filtered(mask_id, :) = filtered2(:, 301:end);
        baseline_A = A - A_filtered;
%         baseline_B = B - B_filtered;
        mean_A = mean(baseline_A, 2);
%         mean_B = mean(baseline_B, 2);
        clear baseline_A baseline_B


%         %%%%%%%%%%%%%%%%%
%         % this part is to compare the traces before and after removing drifting
%         % trend
%         for p = 1
% %         for p = 1 : size(pxlCoords, 1)
%             r = pxlCoords(p, 1)* downSampleRatio - 2 : pxlCoords(p, 1)* downSampleRatio + 2;
%             c = pxlCoords(p, 2)* downSampleRatio - 2 : pxlCoords(p, 2)* downSampleRatio + 2;
%             [R, C] = meshgrid(r, c);
% 
%             p_list = sub2ind(sz(1:2), R, C);
%             p_list = p_list(:);
% 
%             before_avg = mean(B_sliced(p_list, :));
%             after_avg = mean(B_filtered(p_list, :));
% 
%             figure; 
%             subplot(2, 1, 1)
%             plot(before_avg(301:800))
%             subplot(2, 1, 2)
%             plot(after_avg)
%         end 



    frStart = 1;
    frEnd = size(A, 2);

    F = mean(A(mask_id, 300:end), 2);

%   original df/f
    dA = zeros(size(A));
    dA(mask_id, :) = A_filtered(mask_id, :) ./ repmat(F, 1, szZ);
    % dA(mask_id, :) = A(mask_id, :)./repmat(mean(A(mask_id, :), 2), 1, szZ) - 1;

    tmp = dA(mask_id, :);
    s0 = std(tmp(:));
    m0 = mean(tmp(:));
    movieRange = [-2*s0+m0, 5*s0+m0];
    I0 = mat2gray(dA, movieRange);   %scale the whole array so that min = 0, max = 1
    I0 = reshape(I0, sz(1), sz(2), size(A, 2));


    [Iarr0, ~] = gray2ind(I0, 256);
%         Iarr0_small = imresize(Iarr0, 0.5, 'bilinear');
    [~, fn, ~] = fileparts(filename);
%         Iarr2avi(Iarr0_small, frStart, frEnd, ['dff_', fn])
%     Iarr2avi(Iarr0, frStart, frEnd, ['dff_', fn])


    save([fnms{n}(1:end-4), '_preprocessed_dA.mat'], 'dA', 'movieRange', '-v7.3');
        
end

