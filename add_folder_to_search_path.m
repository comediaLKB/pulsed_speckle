%% add folder and subfolders to search path

mfilepath = mfilename('fullpath');
folder = fileparts(mfilepath);
subfolders = genpath(folder);
addpath(subfolders)

clear mfilepath folder subfolders