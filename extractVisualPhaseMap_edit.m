% 09/12/2017 
% 1. Add imgall stimulus red sign on dF/F movies for movies with
% stimulations
% 2. Extract averaged response movie

% files_pre.txt with movie names, spike2 files

% Output: dF/F movie with stim marker, averaged response dF/F (both
% downsampled by 0.5), response matrix



clear; clc;

%cd E:\Lab\Data\withSarah\data\171018
cd Z:\Sarah\testfolder
addpath(genpath('Z:\Sarah\toolbox'))


%addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/piotr_toolbox'))
%addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/wholeBrainDX'))
%addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/sigTOOL'))
%addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/CalciumDX'))
%addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/bfmatlab'))
%addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/chatAnalysis'))
%addpath(genpath('/ysm-gpfs/project/sm2736/toolbox/toolbox/piotr_toolbox'))

filelist = readtext('files_stim.txt', ' ');
fnms = filelist(:, 1);
spike2_fnm = filelist(:, 2);
retin_fnm = filelist(:, 3);
diode_th = filelist(:, 4);
no_movies = length(fnms);
downSampleRatio = 0.5;

stimGap = 6; % min gap should be larger than 6s
sampleFreq = 5000; % Hz

sz = [256 250];
colorRange1 = [-0.02, 0.04];

no_switchingId = 1:2;

% for n = 1
for n = 1:no_movies
    
    clear avgResponseM

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
    
    
    if exist('dA', 'var')
        imgall = dA; 
        clear dA
    else
        imgall = ddiff;
        clear ddiff
    end

    sz(3) = size(imgall, 2);
    imgall = reshape(imgall, sz(1) * sz(2), sz(3));
   


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Import stimulus data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fHeader = fopen(spike2_fnm{n});
    channels = [2 4]; % channel 2: camera frame, channel 4: photoDiode

    
    if isSwitching
        cameraFreq = 20;
    else
        cameraFreq = 10;
    end
    

    for ChannelIndex = 1:length(channels)
        [data_raw{ChannelIndex}, my_header{ChannelIndex}] = SONGetChannel(fHeader,channels(ChannelIndex),'scale');
        ChannelIndex
    end

    isStim = (data_raw{2} < diode_th{n}); % white has lower value for photo diode

    stimOn = find (isStim(2:end) - isStim(1:end-1) == 1) + 1;
    stimOff = find(isStim(2:end) - isStim(1:end-1) == -1) + 1;
    if stimOff(1) < stimOn(1)
        stimOff = stimOff(2:end);
    end
    
    stimOn = stimOn(1:length(stimOff));

    Onsets = [stimOn(1); stimOn(find(stimOn(2:end) - stimOn(1:end-1) > (sampleFreq *stimGap)) + 1)];
    Offsets = [stimOff(stimOn(2:end) - stimOn(1:end-1) > (sampleFreq *stimGap)); stimOff(end)];


    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get average response movie for four directions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(retin_fnm{n});

    frameStep = data_raw{1} * sampleFreq;
    if isSwitching
        frameStep = frameStep(1:2:end);
    end
  
    

    responseM = cell(1, 4);
    t = 1;
    for i = 1 : length(Onsets)
        if t<= length(sweepDir_id)
            OnFrames = find(frameStep >= Onsets(i) & frameStep <= Offsets(i)); 
            if ~isempty(OnFrames) && ~isempty(find(frameStep > Offsets(i)))
                frameOnset(t) = min(OnFrames);
            end
            t=t+1;
        end
    end
    
    for ii = 1 : length(frameOnset)
        
        if isempty(responseM{sweepDir_id(ii)})
        [ii sweepDir_id(ii)]
        responseM{sweepDir_id(t)} = imgall(:, frameOnset(ii) : frameOnset(ii+1)); %replace frameOffset with frameOnset(i+1)
        else    
        [ii sweepDir_id(ii)]
        responseM{sweepDir_id(ii)} = cat(3, responseM{sweepDir_id(ii)}, ...
        imgall(:, frameOnset(ii) : size(responseM{sweepDir_id(ii)}, 2) + frameOnset(ii) - 1));
        end
                   
    end
       



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write movies
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for d = 1 : 4
        avgResponseM{d} = mean(responseM{d}, 3);
        avgResponseM{d} = reshape(avgResponseM{d}, sz(1), sz(2), size(avgResponseM{d}, 2));        
        meanResp_gray{d} = mat2gray(avgResponseM{d}, colorRange1);

        frStart = 1;
        frEnd = size(meanResp_gray{d}, 3);
        [I_resp{d}, ~] = gray2ind(meanResp_gray{d}, 256);
        moviefn = fnm(1 : 9);
        Iarr2avi(I_resp{d}, frStart, frEnd, [moviefn, '_avgD', num2str(d), '.avi']); % dF/F with stim marker
    end


    save([moviefn, '_resp.mat'], 'avgResponseM', 'responseM', 'colorRange1', '-v7.3')
end