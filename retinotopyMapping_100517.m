

sca;
close all;
clearvars;
PsychDefaultSetup(2);

initialize_Param; % prepare parameters




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
screens = Screen('Screens');
screenNumber = max(screens);

white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);

gray = round((white+black)/2);

if gray == white
    gray=white / 2;
end

% Open window
w = PsychImaging('OpenWindow',screenNumber, gray);
rect = Screen('Rect', w);
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

HideCursor;

% Inter frame interval
ifi = Screen('GetFlipInterval', w);


% Compute each frame of the movie and convert the those frames for 4
% sweeping directions separately
tic
f = 1;
for directSweep = 1 : 4
    sweep = singleSweep (directSweep, x_width, y_width, barWidth, numberFrames,...
    numberFramesVert, sizeSq, flickFrame); % generate a single sweep matrix

    for i = 1 : size(sweep, 3)
        img = convertDistortion_largeScreen(sweep(:, :, i), x_width, y_width); % warping
        img(1 : sqrSz, x_width-sqrSz+1 : x_width) = 1; % white square indicating stimulus ON
        tex(f) = Screen('MakeTexture', w, img); 
        f = f + 1; %?
    end
end

% break image, with a black square indicating stimulus off
img = .5 * ones(size(img));
img(1 : sqrSz, x_width-sqrSz+1 : x_width) = 0;
tex(f) = Screen('MakeTexture', w, img);
toc





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Present animation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vbl = Screen('Flip', w); % Sync and get a time stamp

% Use realtime priority for better timing precision:
priorityLevel = MaxPriority(w);
Priority(priorityLevel);


img_rect = [1 1 x_width y_width]; % size of image to show on the screen

% Animation loop:
for sweepId = sweepDir_id
    
    % determine sweep direction
    switch sweepId
        case 1
            f_start = 1;
            f_end = 480;
        case 2
            f_start = 481;
            f_end = 960;
        case 3
            f_start = 961;
            f_end = 1260;
        case 4
            f_start = 1261;
            f_end = 1560;
    end
    
    
    % present images
    tic
    for t = f_start : f_end
        Screen('DrawTexture', w, tex(t), img_rect, rect);
        vbl  = Screen('Flip', w, vbl + waitframes * ifi); % Flip to the screen
    end
    toc
    
    Screen('DrawTexture', w, tex(end), img_rect, rect);
    vbl  = Screen('Flip', w);
    
    WaitSecs(6);
    
end

Priority(0);

Screen('Close');
sca;


