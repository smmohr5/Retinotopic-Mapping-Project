% initial parameters for retinotopic mapping 

fd_name = uigetdir('select saving path');
fn = input('input session name: ', 's');
numSweepPerDir = input('input number of repetitions per direction: ');


ratio = 4; % the downsample ratio to make sweep movies (in order to save time and memory)
x_width = 1920/ratio; % x of small images
y_width = 1200/ratio; % y of small images
barWidth = round(130/ratio); % size of sweeping band, 20 degrees visual angle
numberFrames = 480; % number of frames for horizontal sweep, at 30fps
numberFramesVert = 300; % number of frames for vertical sweep, at 30fps
sizeSq = 164/ratio; % size of checkerboard squares, 25 degrees visual angle
flickFrame = 5; % to get checkerboard flickering every 6 frames
waitframes = 1; % 30fps presentation, 60fps if waitframes = 0
sqrSz = 25; % size of the squre on the top-right corner


% generate sequence for sweeping directions, 10 repetition per direction
rng('shuffle')
sweepDir_id = [ones(1, numSweepPerDir), 2 * ones(1, numSweepPerDir), ...
    3 * ones(1, numSweepPerDir), 4 * ones(1, numSweepPerDir)];
order = randperm(4 * numSweepPerDir);
sweepDir_id = sweepDir_id(order);

save([fd_name, '\', fn, '.mat'], 'sweepDir_id', 'ratio', 'x_width', 'y_width', 'barWidth', ...
    'numberFrames', 'numberFramesVert', 'sizeSq', 'flickFrame', 'waitframes', 'sqrSz', 'numSweepPerDir');


