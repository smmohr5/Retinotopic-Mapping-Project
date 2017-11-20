% modified from make dF/F 030716, to correct the bug from spike2 recordings with the initial protocol 05/07/16 

clear; clc;

% cd 'D:\Lab\Data\wholeBrain\160920_p15'
% cd 'D:\Lab\Data\ChAT\Tra2b\170311_tra2bSnap25_p6_4.0g_female_102x'
% cd 'E:\Lab\Data\wholeBrain\fMRI\170912_visual_paw_stim'

addpath(genpath('/ysm-gpfs/home/xg77/from_louise/CalciumImgCode/piotr_toolbox'))
addpath(genpath('/ysm-gpfs/home/xg77/from_louise/CalciumImgCode/wholeBrainDX'))
addpath(genpath('/ysm-gpfs/home/xg77/from_louise/CalciumImgCode/sigTOOL'))
addpath(genpath('/ysm-gpfs/home/xg77/from_louise/CalciumImgCode/CalciumDX'))
addpath(genpath('/ysm-gpfs/home/xg77/from_louise/CalciumImgCode/bfmatlab'))
addpath(genpath('/ysm-gpfs/home/xg77/from_louise/CalciumImgCode/chatAnalysis'))
addpath(genpath('/ysm-gpfs/home/xg77/from_louise/CalciumImgCode/piotr_toolbox'))


filelist = readtext('files_pre.txt', ' ');
fnms = filelist(:, 1);
blueInitial = filelist(:, 2);
mask_fnms = filelist(:, 3);
no_movies = length(fnms);
downSampleRatio = 0.5;
fixUV = 1;
hat = 300; % window size for rolling hat algorithm



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
        
        
        if fixUV == 1
            sumUV = sum(B);
            dI = sumUV(2:end) - sumUV(1:end-1);
            shortFrameId = find(dI < mean(dI) - 3*std(dI)) + 1;
            if ~isempty(shortFrameId)
                if shortFrameId(end) == szZ 
                    tmpId = shortFrameId(1 : end-1);
                    B(:, tmpId) = (B(:, tmpId - 1) + B(:, tmpId + 1)) / 2;
                    B(:, shortFrameId(end)) = B(:, shortFrameId(end - 1));
                else
                    B(:, shortFrameId) = (B(:, shortFrameId - 1) + B(:, shortFrameId + 1)) / 2;
                end
            end
        end
        
        
        
        ROI = ReadImageJROI(mask_fnms{n});
        
        mask = poly2mask(ROI.mnCoordinates(:, 1), ROI.mnCoordinates(:, 2), sz(1)/downSampleRatio, sz(2)/downSampleRatio);
        mask = imresize(mask, downSampleRatio, 'bilinear');

        [mask_r, mask_c] = find(mask > 0);
        mask_id = find(mask > 0);


        % remove slow drifting trend
        hat = 300;
        se = strel('line', hat, 0);
        A_sliced = A(mask_id, :);
        B_sliced = B(mask_id, :);
        parfor p = 1:length(mask_id)   
        % parfor p = 1:1000 
            filtered1(p, :) = imtophat(A_sliced(p, :), se);
            filtered2(p, :) = imtophat(B_sliced(p, :), se);
        end
        A_filtered = zeros(size(A));
        B_filtered = zeros(size(B));
        A_filtered(mask_id, :) = filtered1;
        B_filtered(mask_id, :) = filtered2;
        baseline_A = A - A_filtered;
        baseline_B = B - B_filtered;
        mean_A = mean(baseline_A, 2);
        mean_B = mean(baseline_B, 2);
        clear baseline_B

        

        frStart = 1;
        frEnd = size(A, 2);

        % use regression to compute gain
        small_A = A_filtered(mask_id, :) + repmat(mean_A(mask_id, :), 1, szZ);
        small_B = B_filtered(mask_id, :) + repmat(mean_B(mask_id, :), 1, szZ);
        parfor i = 1:length(mask_id)
            g_pixel(i) = regress(small_A(i, :)', small_B(i, :)');
        end
        diff = zeros(size(A_filtered));
        diff(mask_id, :) = A_filtered(mask_id, :) - repmat(g_pixel', 1, szZ) .* B_filtered(mask_id, :);


        tmp = diff(mask_id, :);
        s = std(tmp(:));
        m = mean(tmp(:));  
        I=mat2gray(diff, [-1*s+m, 2.5*s+m]);   %scale the whole array so that min = 0, max = 1
        I = reshape(I, sz(1), sz(2), size(A, 2));

        [Iarr, ~] = gray2ind(I, 256);
%         Iarr_small = imresize(Iarr, 0.5, 'bilinear');

%         moviefn = ['regress_', filename];
        moviefn2 = ['regress_dff_', filename];




        F = mean(A(mask_id, :), 2);
        
        % convert difference movie to dff
        ddiff = zeros(size(A));
        ddiff(mask_id, :) = diff(mask_id, :) ./ repmat(F, 1, szZ);
        tmp = ddiff(mask_id, :);
        m = mean(tmp(:));       
        s = std(tmp(:));
        Idff = mat2gray(ddiff, [-2*s + m, 5*s + m]);   %scale the whole array so that min = 0, max = 1
        Idff = reshape(Idff, sz(1), sz(2), size(A, 2));

        [Iarrdff, ~] = gray2ind(Idff, 256); 
        Iarr2avi(Iarrdff, frStart, frEnd, moviefn2)


% %         original df/f
%         dA = zeros(size(A));
%         dA(mask_id, :) = A_filtered(mask_id, :) ./ repmat(F, 1, szZ);
%         % dA(mask_id, :) = A(mask_id, :)./repmat(mean(A(mask_id, :), 2), 1, szZ) - 1;
% 
%         tmp = dA(mask_id, :);
%         s0 = std(tmp(:));
%         m0 = mean(tmp(:));
%         I0=mat2gray(dA, [-2*s0+m0, 5*s0+m0]);   %scale the whole array so that min = 0, max = 1
%         I0 = reshape(I0, sz(1), sz(2), size(A, 2));
% 
% 
%         [Iarr0, ~] = gray2ind(I0, 256);
% %         Iarr0_small = imresize(Iarr0, 0.5, 'bilinear');
%         [~, fn, ~] = fileparts(filename);
% %         Iarr2avi(Iarr0_small, frStart, frEnd, ['dff_', fn])
%         Iarr2avi(Iarr0, frStart, frEnd, ['dff_', fn])


        save([fnms{n}(1:end-4), '_preprocessed.mat'], 'ddiff', '-v7.3');
        
end

