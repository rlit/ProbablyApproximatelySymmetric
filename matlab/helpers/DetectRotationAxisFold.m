function [nFold,representativeConfig,representativeScore] = ...
    DetectRotationAxisFold(lat,lon,Img,roiMask,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tol = 1e7;
maxRotations = ceil(2 * pi / params.minAngleBetweenSymmetries);
configTemplate = [0 0 0 1 lat lon 0];

%%

rollsAll = [];
nRot     = [];

for iRot = 2:maxRotations
    step = 2*pi/iRot;
    vals = step:step:(2*pi-params.minAngleBetweenSymmetries/2);

    rollsAll = [rollsAll vals]; %#ok<AGROW>
    nRot     = [nRot     iRot*ones(1,numel(vals))]; %#ok<AGROW>
end

rollsAll = round(rollsAll*tol)/tol;
[rollsUnq,~,I] = unique(rollsAll','first');

configsUnq = repmat(configTemplate,numel(rollsUnq),1);
configsUnq(:,7) = rollsUnq;
good_epsilon = 0.05;
[~,scores] = CalcFullAndSampledError3D(Img,Img,configsUnq,roiMask,'sampled',good_epsilon);

scoresAll = scores(I);
meanScores = [0 arrayfun(@(i)mean(scoresAll(nRot==i)),2:maxRotations)];

% figure;stem(meanScores)

%% revolutional detection
% veryGoodThresh = 0.03;
% nRolls = find()
% if (all(meanScores < veryGoodThresh) || all(meanScores < 1.5*origScoreFull))
%     configsEx = configTemplate;
%     configsEx(4) = 0; % mark the revolution
%     scoresExFull = origScoreFull;
%     scoresExSampled = origScoreSampled;
%     return
% end

%% another possible detection scheme:
% pcts = prctile(meanScores,[40 90]);
% flag = meanScores>pcts(1) & meanScores<pcts(2);
% std40 = std( meanScores(flag) );
% avg40 = mean(meanScores(flag));
% (meanScores - avg40)/std40

% %%
% figure; plot(meanScores,'*'); hold on;
% plot([1,length(meanScores)],median(meanScores)*[1,1],'m-');
% plot([1,length(meanScores)],median(meanScores)-2*iqr(meanScores)*[1,1],'g-');
% 
%%
% nRolls = find(meanScores < median(meanScores) - 2*iqr(meanScores),1,'last');
%%
bestFoldSoFar = -1;
bestScoreSoFar = meanScores(2);
epsilon = 0.003;

if (mean(meanScores(2:end))<epsilon)
    bestFoldSoFar = -1;
else
    for i = 2 : length(meanScores)
        if (meanScores(i) < median(meanScores) - 2*iqr(meanScores))
            if (meanScores(i)<1.5*bestScoreSoFar+epsilon)
                bestFoldSoFar = i;
            end
            if (meanScores(i)<bestScoreSoFar)
                bestScoreSoFar=meanScores(i);
            end
        end
%         figure(77); subplot(4,5,i);  hold on; plot(meanScores,'*');
%         plot([1,length(meanScores)],median(relevantMeanScores)*[1,1],'m-');
%         plot([1,length(meanScores)],(median(relevantMeanScores)-1*std(relevantMeanScores))*[1,1],'g-');
%         plot([1,length(meanScores)],(median(meanScores)-2*iqr(meanScores))*[1,1],'k-');
%         plot(i,meanScores(i),'sr');
        
    end
end

nFold = bestFoldSoFar;
%% 

representativeConfig = configTemplate;

if (nFold == -1)
% revolutionary
representativeScore = mean(meanScores);
representativeConfig(7) = -2*pi/7;
    nFold = 0;
else
    representativeScore = meanScores(nFold);
    representativeConfig(7) = -2*pi/nFold;
end


