function AddPaths()
%%%%%%%%%%%%%%%%%%%%%
curDir = fileparts(mfilename('fullpath'));

addpath([curDir '\mex']);

addpath([curDir '\matlab']);
addpath([curDir '\matlab\helpers']);
addpath([curDir '\matlab\main_code']);
addpath([curDir '\matlab\visualization']);


addpath([curDir '\3rdParty']);

% addpath(genpath([curDir '\shapes']));


