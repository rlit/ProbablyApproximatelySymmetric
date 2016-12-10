function [optError,bestConfig,bestA,optA,fullError,sampledError,overlapError,statistics] = ...
    FindBestTransformation3D(~,bounds,steps,I1,I2,roiMask,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('params','var') || isempty(params) || ~isstruct(params)
    params = struct();
end
params = GetParams_Sym(params);

% ----------- extract params from struct
quiet  = params.quiet;
silent = params.silent;
%debug  = params.debug;
doVisualization = params.doVisualization;
if ~silent
verboseFun = @(varargin)fprintf(varargin{:});
else
verboseFun = @(varargin)[];
end

epsilon = params.epsilon;
delta   = params.delta;

photometricInvariance = params.photometricInvariance;
minAngleBetweenSymmetries = params.minAngleBetweenSymmetries;
% allowReflections = params.allowReflections;
% onlyReflections = params.onlyReflections;
detectedSymmetries = params.detectedSymmetries;
bnbMaxLevel = params.bnbMaxLevel;
allowRestart = params.allowRestart;

optA = params.optA;

if ~params.allowReflections
    params.onlyReflections = false;
end
%%


[h1,w1,d1] = size(I1); r1x = 0.5*(w1-1); r1y = 0.5*(h1-1); r1z = 0.5*(d1-1);
[h2,w2,d2] = size(I2); r2x = 0.5*(w2-1); r2y = 0.5*(h2-1); r2z = 0.5*(d2-1);
% maxDim = max(size(I1));
DTsampling = 0;
EQUALsampling = 0;

%% generate the Net
[configs,gridSize] = CreateListOfConfigs3D(bounds,steps,params);

verboseFun('>> gridSize is: %d\n',gridSize);

assert(size(configs,1) < 40000000,'more than 40 million configs!');


% initualize statistics
statistics.simulatedGridSizes = [];
statistics.bestDists = [];
statistics.numInsiderConfigs = [];
statistics.numConfigs = [];
statistics.deltas = [];
statistics.goodConfigs = [];
statistics.expandedConfigs = [];
statistics.bnbConfigs = {};
statistics.bnbTimes = [];
statistics.carveTimes = [];
statistics.totTime = -1;
statistics.visTime = -1;



%% main loop

deltaFact = 2;%1.511;
level = 0;
bestDists = [];
bestGridVec = [];
newDelta = delta;
bnbTimes = [];
carveTimes = [];
while (1)
    level = level + 1;
    currFactor = deltaFact^(level-1);


    % 0] generate Theta(1/eps^2) random points (and fresh ones each iteration later on)
    numPoints = round(3/epsilon^2);
    roiIdxs = find(roiMask)';
    idxs = randsample(roiIdxs,numPoints);
    [ys,xs,zs] = ind2sub([h1,w1,d1],idxs);

    % visualize sample
    if (doVisualization && ~(DTsampling || EQUALsampling))
        sampleInds = sub2ind(size(I1),ys,xs,zs);
        insideInds = sampleInds(I1(sampleInds)==1);
        [inside_ys,inside_xs,inside_zs] = ind2sub(size(I1),insideInds);
        outsideInds = sampleInds(I1(sampleInds)~=1);
        [outside_ys,outside_xs,outside_zs] = ind2sub(size(I1),outsideInds);
        if (~quiet&&~silent)
            figure(66); clf;
            hold on; axis equal;
            fv1 = isosurface(I1,0.5);
            patch(fv1,'facecolor','none','edgecolor','c');
            plot3(inside_xs,inside_ys,inside_zs,'.m');
            plot3(outside_xs,outside_ys,outside_zs,'.r');
            view(3); camlight
            ttl = sprintf('number of pixels: %d (%dx%dx%d), inside: %d (%.1f%%)\nnumber of samples: %d, inside: %d (%.1f%%)',...
                numel(I1),h1,w1,d1,sum2(I1==1),100*sum2(I1==1)/numel(I1),...
                numel(sampleInds),numel(insideInds),100*numel(insideInds)/numel(sampleInds));
            title(ttl);
        end
    end


    % 1] translate config vectors to matrix form
    Configs2AffineMEX = tic;
    verboseFun('----- Configs2Affine, with %d configs -----\n',size(configs,1));
    [matrixConfigs_mex, insiders] = ...
        Configs2Affine_mex_3D(configs',int32(h1),int32(w1),int32(d1),int32(h2),int32(w2),int32(d2),...
        int32(r1x),int32(r1y),int32(r1z),int32(r2x),int32(r2y),int32(r2z));
    Configs2Affine_mex_time = toc(Configs2AffineMEX);

    % 1b] remove configs
    CloseConfigsDetection = tic;
    if (exist('minAngleBetweenSymmetries','var'))
        %  too close to previously detected symmetries
        idx_to_keep = FindConfigsFarFromOnesInList(configs,detectedSymmetries,minAngleBetweenSymmetries);
    else
        idx_to_keep = find(insiders);
    end
    CloseConfigsDetection_time = toc(CloseConfigsDetection);
    matrixConfigs_mex = matrixConfigs_mex(:,idx_to_keep);
    origNumConfigs = size(configs,1);

    configs = configs(idx_to_keep,:);


    % 2] evaluate all configurations
    EvaluateConfigsMEX = tic;

    verboseFun('----- Evaluate Configs, with %d configs ----- using %d random points\n',size(configs,1),length(zs));
    distances = EvaluateConfigs_mex_3D(I1,I2,matrixConfigs_mex,int32(xs),int32(ys),int32(zs),int32(photometricInvariance));
    EvaluateConfigs_mex_time = toc(EvaluateConfigsMEX);

    bnbTimes(level) = Configs2Affine_mex_time + EvaluateConfigs_mex_time; %#ok<AGROW> 
    carveTimes(level) = CloseConfigsDetection_time;%#ok<AGROW>

    [bestDist,ind] = min(distances);
    bestConfig = configs(ind,:);
    bestA = CreateAffineTransformation3D(configs(ind,:));
    statistics.bnbConfigs{level} = bestConfig;

    bestDists(level) = bestDist; %#ok<AGROW>


    % for statistics
    [~,simulatedGridSize] = CreateListOfConfigs3D(bounds,ReduceBy(steps,currFactor),params,1); % estimate the size of the complete grid
    statistics.simulatedGridSizes(end+1) = simulatedGridSize;
    statistics.bestDists(end+1) = bestDist;
    statistics.numInsiderConfigs(end+1) = size(matrixConfigs_mex,2);
    statistics.numConfigs(end+1) = origNumConfigs;
    statistics.deltas(end+1) = newDelta;
    if (exist('goodConfigs','var'))
        statistics.goodConfigs(end+1) = size(goodConfigs,1);
    end
    if (exist('expandedConfigs','var'))
        statistics.expandedConfigs(end+1) = size(expandedConfigs,1);
    end

    verboseFun('$$$ bestDist = %.3f\n',bestDist);
    %     if ( (bestDist < 0.005) || ((level > 2) && (bestDist < 0.01)) || (level >= 20) || ...
    %             ((level > 3) && (bestDist > mean(bestDists(level-3:level-1))*0.97)) || ...
    %             ((level > 8) && (statistics.goodConfigs(end) > 20 ) ) )
    %         break
    %     end
    if ( (bestDist < 0.001) || (level >= 20) || ...
            ((level > 1) && statistics.goodConfigs(end) < 2) || ...
            ((level > 6) && (bestDist > mean(bestDists(level-3:level-1))*0.97)) || ...
            ((level > 8) && (statistics.goodConfigs(end) > 20 ) ) || ...
            ((level > 2) && (statistics.numConfigs(end) > 4*statistics.numConfigs(end-1))) || ...
            ((level > 7) && (statistics.goodConfigs(end) > 0.8*statistics.goodConfigs(end-2))) )

        %         maxFullTests = 1000;
        %         nFullTests = min(size(configs,1),maxFullTests);
        %
        %         [~,best4Full] = sort(distances);
        %         best4Full = best4Full(1:nFullTests);
        %         fullErrors = CalcFullError3D(I1,I2,configs(best4Full,:),roiMask);
        %         [~,ind] = min(fullErrors);
        %
        %         bestConfig = configs(best4Full(ind),:);
        %         bestA = CreateAffineTransformation3D(configs(ind,:));
        break
    end
    %     bestConfig = configs(ind,:);
    %     bestA = CreateAffineTransformation3D(bestConfig);


    % break

    % 5] choose the 'surviving' configs and delta for next round
    [~,goodConfigs,repeatOrigGrid] = GetGoodConfigsByDistance3D(configs,bestDist,newDelta,distances,bestGridVec);

    if (isempty(goodConfigs))
        keyboard
    end
    numConfigs = size(goodConfigs,1);

    
    % 7] expand 'surviving' configs for next round
    if (repeatOrigGrid && (bestDist > 0.05) &&(level==1) && (origNumConfigs < 7.5*10^6))
        % don't keep restarting if don't see improvement
         if (exist('ppBestDist','var') && bestDist > 0.9*ppBestDist)
       doRestart = 0;
         else
             doRestart = 1;
         end
    else
        doRestart = 0;
    end

    if (doRestart && allowRestart)
        fact = 0.9;
        verboseFun('##### RESTARTING!!! changing from delta: %.3f, to delta: %.3f\n', newDelta, newDelta*fact);
        if (exist('pBestDist','var'))
            ppBestDist = pBestDist;
        end
        pBestDist = bestDist;
        newDelta = newDelta*fact;
        level = 0;
        steps.tx = fact*steps.tx;
        steps.ty = fact*steps.ty;
        steps.tz = fact*steps.tz;
        steps.r = fact*steps.r;
        steps.s = fact*steps.s;
        [configs,~] = CreateListOfConfigs3D(bounds,steps,params);

    elseif level >= bnbMaxLevel
        verboseFun('##### MAX B&B LEVEL REACHED!!! finishing...\n');
        break

    else
        prevDelta = newDelta;
        newDelta = newDelta/deltaFact;
        verboseFun('##### CONTINUING!!! prevDelta = %.3f,  newDelta = %.3f \n',prevDelta,newDelta);

        % expand the good configs
        expandType = 'fullExpansion'; % 'randomExpansion'; %    'deltaGrid'; %
        switch expandType
            case 'randomExpansion'
                expandedConfigs = ExpandConfigsRandom3D(goodConfigs,steps,bounds,level,100,deltaFact);
            case 'fullExpansion'
                expandedConfigs = ExpandConfigsFull3D(goodConfigs,steps,bounds,level,deltaFact);
        end

        configs = [goodConfigs ; expandedConfigs];
        if (any(isnan(expandedConfigs(:))))
            error('there are NaNs here!')
        end
    end


    %     configs = unique(configs,'rows'); % REMOVED THIS - IT IS WORTHWHILE

    verboseFun('***\n');
    verboseFun('*** level %d:|goodConfigs| = %d, |expandedConfigs| = %d\n',level,numConfigs,size(configs,1));
    verboseFun('***\n');


end

%% for debug

optError = 0;
overlapError = [];

if params.calcFullError
    [fullError,~] = CalcFullAndSampledError3D(I1,I2,bestConfig,roiMask,'full',epsilon);
else
    fullError = NaN;
end

%% for output
sampledError = bestDist;
statistics.level = level;


statistics.bnbTimes = bnbTimes;
statistics.carveTimes = carveTimes;
statistics.totTime = sum(bnbTimes)+sum(carveTimes);


return



function reduced = ReduceBy(steps,currFactor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reduced.tx = steps.tx/currFactor;
reduced.ty = steps.ty/currFactor;
reduced.tz = steps.tz/currFactor;
reduced.s = steps.s/currFactor;
reduced.r = steps.r/currFactor;



function [res,i] = IsMemberApprox(A,row,err)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res = 0;
for i = 1 : size(A,1)
    if (norm(A(i,:)-row) < err)
        res = 1;
        return
    end
end


function [thresh,goodConfigs,repeatOrigGrid] = GetGoodConfigsByDistance3D(configs,bestDist,newDelta,distances,bestGridVec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
repeatOrigGrid = 0;
% targetNum = 20000;
% thresh = bestDist + newDelta/3;
thresh = bestDist + GetThreshPerDelta3D(newDelta);
goodConfigs = configs(distances <= thresh, :); % bestDist + levelPrecision,:);
numConfigs = size(goodConfigs,1);
while (numConfigs > 23000)
    thresh = thresh * 0.99;
    disp('CHANGING THRESHOLD - TOO MANY CONFIGS!');
    goodConfigs = configs(distances <= thresh, :); % bestDist + levelPrecision,:);
    numConfigs = size(goodConfigs,1);
end

if (thresh < min(distances))
    thresh = min(distances);
    goodConfigs = configs(distances <= thresh, :); % bestDist + levelPrecision,:);
    numConfigs = size(goodConfigs,1);
    if (numConfigs > 10000)
        goodConfigs = goodConfigs(randsample(numConfigs,10000),:);
        numConfigs = 10000;
    end
end

percentage = numConfigs/size(configs,1);
while (percentage < 1/500)
    thresh = thresh / 0.99;
    goodConfigs = configs(distances <= thresh, :); % bestDist + levelPrecision,:);
    numConfigs = size(goodConfigs,1);
    percentage = numConfigs/size(configs,1);
end


if (percentage > 0.022)
    repeatOrigGrid = 1;
end

if (~isempty(bestGridVec))
    [exists,~] = IsMemberApprox(goodConfigs,bestGridVec,1000*eps);
    if (~exists)
        warning('problem with configs');
    end
end


