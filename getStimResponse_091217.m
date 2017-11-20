% 09/12/2017 
% 1. Add imgall stimulus red sign on dF/F movies for movies with
% stimulations
% 2. Extract averaged response movie

% files_pre.txt with movie names, spike2 files

% Output: dF/F movie with stim marker, averaged response dF/F (both
% downsampled by 0.5), response matrix



clear; clc;

% cd E:\Lab\Data\wholeBrain\fMRI\170912_visual_paw_stim

addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/piotr_toolbox'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/wholeBrainDX'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/sigTOOL'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/CalciumDX'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/bfmatlab'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/chatAnalysis'))
addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/piotr_toolbox'))


filelist = readtext('files_stim.txt', ' ');
fnms = filelist(:, 1);
stim_fnm = filelist(:, 2);
no_movies = length(fnms);
downSampleRatio = 0.5;

sz = [256 250];
colorRange1 = [0, 0.05];
colorRange2 = [-0.02, 0.1];

no_switchingId = 6 : 10;

for n = 2:no_movies

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Import tiff movie and get dFofF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fnm = fnms{n};
    load(fnm)
    
    if ismember(n, no_switchingId)
        isSwitching = 0;
    else
        isSwitching = 1;
    end
    
%     if exist('dA', 'var')
%         imgall = dA; 
%         clear dA
%         isSwitching = 0;
%     else
%         imgall = ddiff;
%         clear ddiff
%         isSwitching = 1;
%     end
    imgall = dA; clear dA

    sz(3) = size(imgall, 2);
    imgall = reshape(imgall, sz(1) * sz(2), sz(3));
   


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Import stimulus data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fHeader = fopen(stim_fnm{n});
    channels = [2 4]; % channel 1: camera frame, channel 2: stim

    th = 1;
    stimGap = 2; % considered to be two separate stimuli if their time difference is larger than 1s
    sampleFreq = 1000; % Hz
    if isSwitching
        cameraFreq = 20;
    else
        cameraFreq = 10;
    end
    

    for ChannelIndex = 1:length(channels)
        [data_raw{ChannelIndex}, my_header{ChannelIndex}] = SONGetChannel(fHeader,channels(ChannelIndex),'scale');
        ChannelIndex
    end

    isStim = (data_raw{2} > th);

    stimOn = find (isStim(2:end) - isStim(1:end-1) == 1) + 1;
    stimOff = find(isStim(2:end) - isStim(1:end-1) == -1) + 1;
    if stimOff(1) < stimOn(1)
        stimOff = stimOff(2:end);
    end
    
    stimOn = stimOn(1:length(stimOff));

    Onsets = [stimOn(1); stimOn(find(stimOn(2:end) - stimOn(1:end-1) > (sampleFreq *stimGap)) + 1)];
    Offsets = [stimOff(stimOn(2:end) - stimOn(1:end-1) > (sampleFreq *stimGap)); stimOff(end)];

    Offsets = Offsets(Offsets < sampleFreq * 630);
    Onsets = Onsets(1:length(Offsets));
%     Offsets = Offsets(2:end);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add stim marker and get average response movie
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    frameStep = data_raw{1} * sampleFreq;
    if isSwitching
        frameStep = frameStep(1:2:end);
    end
    frameOn = -20;
    frameOff = 180;
    
    [stimMarker_x, stimMarker_y] = meshgrid(1 : 30, sz(2)-29 : sz(2));
    stimMarker_id = sub2ind(sz(1:2), stimMarker_x, stimMarker_y);
    for i = 1 : length(Onsets)
        OnFrames = find(frameStep >= Onsets(i) - sampleFreq/cameraFreq & frameStep <= Offsets(i)); 
        if isempty(OnFrames)
            break
        end
        frameOnset(i) = min(OnFrames);
        if isempty(find(frameStep > Offsets(i)))
            break
        else
            frameOffset(i) = find(frameStep > Offsets(i), 1);
        end
        imgall(stimMarker_id, frameOnset(i) : frameOffset(i)) = 10;
        
        if (i + 1 < length(Onsets))
            avgResponseM(:, :, i) = imgall(:, frameOnset(i) + frameOn - 1 : frameOnset(i) + frameOff - 1);
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write movies
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stimLength = frameOffset(1) - frameOnset(1);
    meanResp = mean(avgResponseM, 3);
    meanResp = reshape(meanResp, sz(1), sz(2), size(meanResp, 2));        
    meanResp2 = mat2gray(meanResp(:, :, 1:180), colorRange1);
    [I_resp, ~] = gray2ind(meanResp2, 256);
    

    imgall = reshape(imgall, sz(1), sz(2), sz(3));        
    imgall2 = mat2gray(imgall, colorRange2);
    [I, ~] = gray2ind(imgall2, 256);
    
    rgbColors = jet(256);
   
    frStart = 1;
    frEnd = sz(3);
    moviefn = fnms{n}(1:end-17);
    Iarr2avi(I, frStart, frEnd, [moviefn, '_addStim.avi']); % dF/F with stim marker
    Iarr2avi(I_resp, frStart, size(meanResp2, 3), [moviefn, '_meanResponse.avi']) % averaged response
    
    save([moviefn, '_resp.mat'], 'avgResponseM', 'meanResp', 'colorRange1', 'sz')
end