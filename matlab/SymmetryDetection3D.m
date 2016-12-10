function [detectedSymmetries,statisticsArray] = ...
    SymmetryDetection3D(shapeIn,shapeName,params,saveWSdir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath(genpath('../shapes/'));

% close all
% clear all
% clc
dbstop if error
set(0,'DefaultFigureWindowStyle','docked')

rand( 'seed',0) %#ok<RAND>
randn('seed',0) %#ok<RAND>

rootDir = fileparts(fileparts(mfilename('fullpath')));
if ~exist('FindBestTransformation3D.m','file')
    run(fullfile(rootDir,'AddPaths.m'))
end

if ~exist('params','var') || isempty(params) || ~isstruct(params)
    params = struct();
end
params = GetParams_Sym(params);


if ~exist('saveWSdir','var')
    saveWSdir = [rootDir '\DefaultSaveFolder\'];
end
if ~isdir(saveWSdir)
    mkdir(saveWSdir)
end

if ~exist('shapeName','var')
    shapeName = 'DefaultShapeName';
end
%% reading/creating the shape Iin

if (params.defineYourOwnShape)
    close all
        shapeName = 'icosahedron';
%     shapeName = 'dodecahedron';

    params.breakSearchByThreshold = true;
    params.expandDetectedSymmetry = false;
    params.autoDetectK = 0;
    params.K = 20;
    params.epsilon = 0.13;

    [I1,I1_radius] = GetPredifinedVolume(shapeName);
    shapeIn = I1;
else
    if (~exist('shapeIn','var'))
        % load a prepared instance
        l = load(fullfile(rootDir,'shapes\COSEG\SELECTED','27.vox.mat'));
        Iin= l.G;
    elseif isstruct(shapeIn)
        Iin = shapeIn.G;
    else
        Iin = shapeIn;
    end
    [I1,I1_radius] = BoxifyShape(Iin);

end
% don't know the GTmaxDim = max(size(I1));
optVec = []; optMat = [];

%% create (truncated) DT volume
if (params.autoDetectK)
    kSTR = 'auto';
else
    K = params.K;
    kSTR = num2str(K);
end
fileNamePrefix = [saveWSdir shapeName '_k_' kSTR '_'];

% check if a file with prefix exist, if not create a tmep one (will be deleted later)
if params.returnIfResultFileExist && numel(dir([fileNamePrefix '*.mat'])) > 0 && ~params.defineYourOwnShape
    [detectedSymmetries,statisticsArray] = deal([]);
    return
else
    save([fileNamePrefix '.mat'],'params')
end

if (params.autoDetectK)
    KvsATVstats = load([rootDir '\KvsATVstatsLinear.mat'],'Ks','avgATVs','stdATVs');
    [Img,roiMask,bounds,steps,ATV,K] = AutoDetectKandGetVolumeMaskAndAV(I1,params,I1_radius,KvsATVstats);
else
    [Img,roiMask,bounds,steps,ATV] = GetVolumeMaskAndAV(K,I1,params,I1_radius);
end

% get the binary version for later calculations
[binImg,binroiMask,binbounds,binsteps,binATV] = GetVolumeMaskAndAV(0,I1,params,I1_radius);
onlyGridSize = 1;
[~,bingridSize] = CreateListOfConfigs3D(binbounds,binsteps,params,onlyGridSize);



% debug
if params.overrideWithBinaryImage
    Img     = binImg;
    roiMask = binroiMask;
end

%params.allowReflections = 0;
[~,gridSize] = CreateListOfConfigs3D(bounds,steps,params,onlyGridSize);

fprintf('K=%3d, ATV=%.3f, gridSize=%d\n',K,ATV,gridSize)

filename = sprintf('%sKvsATV__%s__K_%d_ATV_%.4f_gridSize_%d.mat',saveWSdir,shapeName,K,ATV,gridSize);
save(filename,'shapeName','K','ATV','gridSize');


%% run the search!
photometricInvariance = 0;


% config form: [tx,ty,tz,s,latitude,longitude,roll]
% identityConfig = [0,0,0,1,0,0,0];

detectedSymmetries = []; % identityConfig;
nDetected = 0;
goodSymmetryThreshold = params.goodSymmetryThreshold;

statisticsArray = {};
allErrors = [];

for i = 1 : 1000
    if (params.constantNumberOfBestSymmetries)
        goodSymmetryThreshold = 1;
        if (nDetected >= params.maxNumberOfBestSymmetries)
            break
        end
    end

    % main call - "FindBestTransformation3D"
    params.detectedSymmetries = detectedSymmetries;
    [optError,bestConfig,bestA,optMat,fullError,sampledError,overlapError,statistics] = ...
        FindBestTransformation3D(shapeName,bounds,steps,Img,Img,roiMask,params);


    % collect statistics
    statisticsArray{end+1} = statistics;
    if params.calcFullError
        curError = fullError;
    else
        curError = sampledError;
    end

    if (params.breakSearchByThreshold && curError > goodSymmetryThreshold)
        break
    end

    %     nonNegErrors = allErrors(allErrors>0.01);
    %     if ((curError > goodSymmetryThreshold) || ...
    %         (length(nonNegErrors)>1&&(curError > 3*mean(allErrors(1:end-1)))) || ...
    %         (length(nonNegErrors)==1&&(curError > 5*mean(nonNegErrors))) )
    %         break
    %     end

    allErrors(end+1) = curError;

    [binfullError,binsampledError]  = ...
        CalcFullAndSampledError3D(I1,I1,bestConfig,binroiMask,'both',params.epsilon);

    % add naively to detected list
    currentSymmetry.config = bestConfig;
    currentSymmetry.isReflection = bestConfig(4) < 0; % otherwise 'rotation'
    currentSymmetry.rotationFold = -1;
    currentSymmetry.binFullError = binfullError;
    currentSymmetry.binSampled = binsampledError;
    currentSymmetry.fullError = fullError;
    currentSymmetry.sampledError = sampledError;

    detectedSymmetries = ManageNewSymmetry(detectedSymmetries,currentSymmetry,Img,roiMask,params);

    lastError = detectedSymmetries(end).fullError;
    if (params.breakSearchByThreshold && lastError > goodSymmetryThreshold)
        detectedSymmetries = detectedSymmetries(1:end-1);
        nDetected = numel(detectedSymmetries);
        break
    end

    nDetected = numel(detectedSymmetries);

    fprintf('Detected Symmetry #%d:',i)
    curConf = 180*detectedSymmetries(end).config/pi;
    fprintf('reflection=%d, lat=%6.1f, lon=%6.1f, roll=%6.1f, fold:%d, fullBinErr=%5.3f\n',...
        curConf(4)<0,curConf(5:7),detectedSymmetries(end).rotationFold,detectedSymmetries(end).binFullError)



end

%% visualization (and image saving) loop
if params.doVisualization && params.plot_symmetries
    plot_symmetries(shapeIn,shapeName,detectedSymmetries,fileNamePrefix,params)
end

%%
debug = 0;
if (debug)
    sumTime = zeros(1,nDetected);
    for j = 1 : nDetected
        sumTime(j) = statisticsArray{j}.totTime;
    end
    sumTime
    sum(sumTime)
    mean(sumTime)
    save('here.mat')

    g = figure; plot(1:length(sumTime),cumsum(sumTime),'-*');
    xlim([0.1,49])
    xlabel('detected symmetry number','fontweight','bold','fontsize',11');
    ylabel('total (accumulated) time in seconds','fontweight','bold','fontsize',11');
    set(gca,'fontweight','bold');

    % inset
    img = imread([fileNamePrefix 'dode.png']);
    [h,w,d] = size(img);
    img = img(h/5:4*h/5,w/3.1:2.2*w/3.1,:);% imshow2(img)
    axes('position',[.47 .13 .4 .45]);
    image(img);axis off

    paperDir = '..\..\latex\paper_v1\SK\images\dodec_runtimes\';
    mkdir(paperDir);
    saveas(g,[paperDir 'dodec.eps'],'epsc');
    saveas(g,[paperDir 'dodec.png']);

    % and the dodec figure
    fv1 = isosurface(I1,0.5);
    dod_h = figure; hold on; axis equal;
    patch(fv1,'facecolor','r','edgecolor','none');
    view(3); camlight
    axis off
    set(gcf,'color','w');
    saveas(dod_h,[fileNamePrefix 'dode.png'])
end
%% saving
close all
statisticsArray = [statisticsArray{:}];
save([fileNamePrefix '_info_.mat'],'statisticsArray','K','ATV','binATV','gridSize','bingridSize',...
    'detectedSymmetries','params');


% delete temp file
if exist([fileNamePrefix '.mat'],'file')
    delete([fileNamePrefix '.mat'])
end


function plot_symmetries(shapeIn,shapeName,detectedSymmetries,fileNamePrefix,params)
%%%%%%%%%%%%%%%%%%%%%%%%%
close all
nDetected = numel(detectedSymmetries);
for symmInd = 1 : nDetected

    bestConfig = detectedSymmetries(symmInd).config;
    bestA = CreateAffineTransformation3D(bestConfig);
    sampledError = detectedSymmetries(symmInd).sampledError;
    fullError = detectedSymmetries(symmInd).fullError;

    numPoints = []; xs = []; ys = []; zs = []; optA = [];
    VisualizeSamplingGrayscale3D(strrep(shapeName,'_',' '),shapeIn,shapeIn,numPoints,...
        bestA,bestA(1,4),bestA(2,4),bestA(3,4),...
        xs,ys,zs,bestConfig,optA,sampledError,fullError,detectedSymmetries(symmInd).rotationFold);


    %% saving
    axis off
    saveas(gcf,[fileNamePrefix '_symmetry_number_' num2str(symmInd) '.png']);
    saveas(gcf,[fileNamePrefix '_symmetry_number_' num2str(symmInd) '.fig']);

    if params.createEpsFile
        saveas(gcf,[fileNamePrefix '_symmetry_number_' num2str(symmInd) '.eps'],'epsc');
    end

    if params.createBlenderFile
        axis off;title('');
        FixAndExportToBlender(gcf,[fileNamePrefix '_symmetry_number_' num2str(symmInd) '.wrl'])
    end

end
close all







