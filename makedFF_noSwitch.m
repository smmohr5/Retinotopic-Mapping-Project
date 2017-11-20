% preprocess (mostly just downsample and save into preprocessed.mat data)
% for movies with no switching
% 02/27/17

% 09/13/17 removed -1 when computing dA


clear; clc;

% cd 'E:\Lab\Data\wholeBrain\fMRI\170914_emxG6_p33_male_21.2g_40.5x'


addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/piotr_toolbox'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/wholeBrainDX'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/sigTOOL'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/CalciumDX'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/bfmatlab'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/chatAnalysis'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/piotr_toolbox'))




filelist = readtext('files_pre2.txt', ' ');
fnms = filelist(:, 1);
mask_fnms = filelist(:, 2);
no_movies = length(fnms);
downSampleRatio = 0.5;
hat = 300; % window size for rolling hat algorithm
% hat = 1200; % longer hat for faster switching

highSampling = 0;

for n = 1:no_movies
    
    clear mask_id filtered1 filtered2 A_sliced B_sliced A_filtered B_filtered diff 
    
    filename = fnms{n};
    concatList = dir(fullfile([filename(1 : end-4), '*.tif']));
    
        
    A = [];
    for c = 1:min(6, length(concatList))
        imgall = openMovie(concatList(c).name);
        szall = size(imgall);
        initial = 2;
        if highSampling
            if mod(c, 2) == 1
                FramesId = initial : 2: szall(3);
            else
                FramesId = (3 - initial) : 2: szall(3);
            end
            
            A = cat(3, A, imresize(imgall(:, :, FramesId), downSampleRatio, 'bilinear'));
        else
            A = cat(3, A, imresize(imgall, downSampleRatio, 'bilinear'));
        end
        
        clear imgall
    end



    sz = size(A); szZ=sz(3);
    npix = prod(sz(1:2));
    A = reshape(A, npix, szZ); %reshape 3D array into space-time matrix             


    ROI = ReadImageJROI(mask_fnms{n});

    mask = poly2mask(ROI.mnCoordinates(:, 1), ROI.mnCoordinates(:, 2), sz(1)/downSampleRatio, sz(2)/downSampleRatio);
    mask = imresize(mask, downSampleRatio, 'bilinear');

    [mask_r, mask_c] = find(mask > 0);
    mask_id = find(mask > 0);


    % remove slow drifting trend
    hat = 300;
    se = strel('line', hat, 0);
    A_sliced = A(mask_id, :);
    parfor p = 1:length(mask_id)   
    % parfor p = 1:1000 
        filtered1(p, :) = imtophat(A_sliced(p, :), se);
    end
    A_filtered = zeros(size(A));
    A_filtered(mask_id, :) = filtered1;
    baseline_A = A - A_filtered;
    mean_A = mean(baseline_A, 2);
    clear baseline_A 


    frStart = 1;
    frEnd = size(A, 2);


    F = mean(A(mask_id, :), 2);




%   df/f
    dA = zeros(size(A));
    dA(mask_id, :) = A_filtered(mask_id, :) ./ repmat(F, 1, szZ);
    % dA(mask_id, :) = A(mask_id, :)./repmat(mean(A(mask_id, :), 2), 1, szZ) - 1;

    tmp = dA(mask_id, :);
    s0 = std(tmp(:));
    m0 = mean(tmp(:));
    I0=mat2gray(dA, [-2*s0+m0, 5*s0+m0]);   %scale the whole array so that min = 0, max = 1
    I0 = reshape(I0, sz(1), sz(2), size(A, 2));


    [Iarr0, ~] = gray2ind(I0, 256);
    [~, fn, ~] = fileparts(filename);
    Iarr2avi(Iarr0, frStart, frEnd, ['dff_', fn])

    
    movieRange = [-2*s0+m0, 5*s0+m0];
    save([fnms{n}(1:end-4), '_preprocessed.mat'], 'dA', 'movieRange', '-v7.3');
        
end

