function farIdxs = FindConfigsFarFromOnesInList(configs,listToBeFarFrom,angleTh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% config form: [tx,ty,tz,s,latitude,longitude,roll]
colScale = 4; % we actually use this for reflection (denoted by -1)
colLat   = 5;
colLon   = 6;
colRoll  = 7;
isRollHigh = abs(configs(:,colRoll)) > angleTh;
isReflected = (configs(:,colScale) < 0);
isFar = isRollHigh | isReflected; % small roll is too close to the identity transform
for ii = 1:numel(listToBeFarFrom)
    curConfig = listToBeFarFrom(ii).config;

    sphereDist = greatcircledist(configs(:,colLat),configs(:,colLon),curConfig(colLat),curConfig(colLon));
    sphereDist = abs(sphereDist);
    % sphereDist in [0,pi]

    if listToBeFarFrom(ii).isReflection
        isSameReflection = isReflected;
    else
        isSameReflection = ~isReflected;
    end

    % handle close pairs
    isNearOnSameSide  = (sphereDist <    angleTh) & isSameReflection;
    isNearOnOtherSide = (sphereDist > pi-angleTh) & isSameReflection;

    isFar = isFar & ...
        ~(isNearOnSameSide | isNearOnOtherSide);

    if ~listToBeFarFrom(ii).isReflection && listToBeFarFrom(ii).rotationFold < 2
        % this is a continous srotational symmetry - remove reflections that are orthogonal to this axis
        isOrthogonal = abs(sphereDist - pi/2) < angleTh & isReflected;
        isFar = isFar & ~isOrthogonal;
    end

end

% remainingConfigs = configs(isFar,:);
% removedConfigs = configs(~isFar,:);
farIdxs = find(isFar);


