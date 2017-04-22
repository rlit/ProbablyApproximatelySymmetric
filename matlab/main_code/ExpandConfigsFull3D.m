function expandedConfigs = ExpandConfigsFull3D(configs,steps,bounds,level,deltaFact)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stepFactor = deltaFact^level;

expandedConfigs = ExpandSingleConfig(configs,steps,bounds,stepFactor);

% % --------- Alternative implementation
% configCell = num2cell(configs,2);
% configCellEx = cellfun(@(c)ExpandSingleConfig(c,steps,bounds,stepFactor),configCell, 'UniformOutput',0);
% configEx = vertcat(configCellEx{:});
% isequal(sortrows(configEx),sortrows(expandedConfigs))

%% bound expansions to valid limits
expandedConfigs(:,1) = BoundBy(expandedConfigs(:,1),bounds.tx(1),bounds.tx(2));
expandedConfigs(:,2) = BoundBy(expandedConfigs(:,2),bounds.ty(1),bounds.ty(2));
expandedConfigs(:,3) = BoundBy(expandedConfigs(:,3),bounds.tz(1),bounds.tz(2));
% expandedConfigs(:,4) = BoundBy(expandedConfigs(:,4),bounds.s);
expandedConfigs(:,5) = BoundBy(expandedConfigs(:,5),0,pi);%lat
expandedConfigs(:,6) = BoundBy(expandedConfigs(:,6),0,2*pi);%lon
expandedConfigs(:,7) = BoundBy(expandedConfigs(:,7),-pi,0);%roll

% disp('done fully expanding the good configs');


function configsEx = ExpandSingleConfig(configIn,steps,bounds,stepFactor)
configsEx = configIn;


% expand X+Y+Z
if range(bounds.tx)>0
    configsEx = ExpandSingleCoordinate(configsEx,steps.tx/stepFactor,1);
end
if range(bounds.ty)>0
    configsEx = ExpandSingleCoordinate(configsEx,steps.ty/stepFactor,2);
end
if range(bounds.tz)>0
    configsEx = ExpandSingleCoordinate(configsEx,steps.tz/stepFactor,3);
end

if range(bounds.r)>0
    % latitude
    configsEx = ExpandSingleCoordinate(configsEx,steps.r/stepFactor,5);

    % longitude
    origLat = configsEx(:,5);
    totPatchAreaAtLat = 2*pi*(cos(origLat-steps.r) - cos(origLat));
    numPatchesAtInclinations = max(1,round(totPatchAreaAtLat/(steps.r)^2)); % prevent illegal inclinations
    halfThetaSteps = 2*pi./numPatchesAtInclinations;
    configsEx = ExpandSingleCoordinate(configsEx,halfThetaSteps,6);

    % roll
    isReflection  = configsEx(:,4) < 0;
    refConfigs    = configsEx(isReflection,:);
    nonRefConfigs = configsEx(~isReflection,:);
    nonRefConfigsEx = ExpandSingleCoordinate(nonRefConfigs,steps.r/stepFactor,7);

    configsEx = [refConfigs ; nonRefConfigsEx];

end


function configsEx = ExpandSingleCoordinate(configIn,step,coordIdx)
configWithStep1 = configIn;
configWithStep1(:,coordIdx) = configWithStep1(:,coordIdx) + step;

configWithStep2 = configIn;
configWithStep2(:,coordIdx) = configWithStep2(:,coordIdx) - step;

configsEx = [configWithStep1 ; configWithStep2 ; configIn];


