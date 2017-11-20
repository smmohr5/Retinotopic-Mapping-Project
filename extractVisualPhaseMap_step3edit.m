clear; clc;

% cd Z:\Sarah\171030_Rx_p17_male_7.0g
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/piotr_toolbox'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/wholeBrainDX'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/sigTOOL'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/CalciumDX'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/bfmatlab'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/chatAnalysis'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/piotr_toolbox'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/Utilities'))

% cd E:\Lab\Data\withSarah\data\171010
% addpath(genpath('/ysm-gpfs/home/xg77/from_louise/CalciumImgCode/piotr_toolbox'))
% addpath(genpath('/ysm-gpfs/home/xg77/from_louise/CalciumImgCode/wholeBrainDX'))
% addpath(genpath('/ysm-gpfs/home/xg77/from_louise/CalciumImgCode/sigTOOL'))
% addpath(genpath('/ysm-gpfs/home/xg77/from_louise/CalciumImgCode/CalciumDX'))
% addpath(genpath('/ysm-gpfs/home/xg77/from_louise/CalciumImgCode/bfmatlab'))
% addpath(genpath('/ysm-gpfs/home/xg77/from_louise/CalciumImgCode/chatAnalysis'))
% addpath(genpath('/ysm-gpfs/home/xg77/from_louise/CalciumImgCode/piotr_toolbox'))

w = 52;  % width of screen, in cm
h = 33;  % height of screen, in cm
z_dist = 10; % distance from the eye to screen

mask_fn = '171023_mask.roi';
% mask_fn = '171015_mask2.roi';
% mask_fn = 'ROI_retin_171011.roi';
% mask_fn = '171030_viscor.roi';
% mask_fn = '171030_cran.roi';
file_list = dir(fullfile('*_resp.mat'));

downSampleRatio = .5; % down sample ratio for preprocessed movies
order = 3; % butterworth filter order
fs = 10; % imaging sampled at 10hz
sz = [256 250];

ROI = ReadImageJROI(mask_fn);     


for n = 1:4 %length(file_list)
    
    load(file_list(n).name)
    mask = poly2mask(ROI.mnCoordinates(:, 1), ROI.mnCoordinates(:, 2), sz(1), sz(2));
    mask_id = find(mask > 0);
    
    
    for d = 1 : length(avgResponseM)
        
        sz_d = size(responseM{d}); 
        concat_movie{d} = reshape(responseM{d}, sz_d(1), sz_d(2) * sz_d(3)); 
        fr_n =  sz_d(2) * sz_d(3);
        k = sz_d(3); % number of repetitions in one movie
        
        
        stimFrame = size(avgResponseM{d}, 3);
        phaseMap{d} = zeros(size(avgResponseM{d}(:, :, 1)));
        Kmap{d} = zeros(size(avgResponseM{d}(:, :, 1)));
        

        nfft = fr_n;
        f = fs / 2 * linspace(0, 1, nfft / 2 + 1); % choosing correct frequency axes
        tmp = f - fs/sz_d(2);
        freq_id = find(abs(tmp) == min(abs(tmp)));
        phase = zeros(size(mask_id));
        for p = 1 : length(mask_id)   
            y = concat_movie{d}(mask_id(p), :);
            if mean(y) ~= 0
                
                res = sum(y.*exp(-1i * 2 * pi * k * (1:nfft)/nfft));
                k_value{d}(p) = res;
                
                include_id(p) = 1;
            end
        end     
        
        Kmap{d} = zeros(sz);
        Kmap{d}(mask_id(include_id > 0)) = k_value{d}(include_id > 0);
        
       
        A = angle(Kmap{d})*180/pi;
            B = find(A < -90);
            A(B) = A(B)+ 360;
                
     %with interstime frames phase angles vert:(-90, 135) horz:(-90, 171.82)    
        if d == 1 || d==2
            A = (A/1.6314)-13.793;
        else
            A = (A/2.23)-18.343;
        end
        
%         if no interstim frames (-90 - 270 degrees)
%         if d == 1 || d==2
%             A = (A/2.6102)-34.48;
%         else
%             A = (A/3.0664)-29.35;
%         end
%         
        
        
        h = figure; imagesc((A)); colormap HSV; axis image; colorbar %test
                                    
        
        title(['original kmap d', num2str(d)])
        saveas(h, [file_list(n).name(1:end-4), '_phase_d', num2str(d), '.png'])
        
    end

    
    save([file_list(n).name(1:end-4), '_Kmap.mat'], 'Kmap')
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % combine kmaps for horizontal vs. vertical, smoothing for kmaps
    for lowpass = [0.5, 0.75, 1]
        [kmap_hor kmap_vert delay_hor delay_vert magS] = ...
            Gprocesskret_batch(Kmap, mask, lowpass, []); % refer to the original function for input arguments
        h = figure; 
        set(h, 'Position', [0 0 1000 800])
        subplot(1, 2, 1);imagesc(kmap_hor); colormap HSV; axis image; colorbar %test
        title(['kmap_hor lowpass', num2str(lowpass)])
        subplot(1, 2, 2);imagesc(kmap_vert); colormap HSV; axis image; colorbar %test
        title(['kmap_vert lowpass', num2str(lowpass)])
        
        saveas(h, [file_list(n).name(1:end-4), '_lowpas', num2str(lowpass), '.png'])
    end  
        
end

close all
    
