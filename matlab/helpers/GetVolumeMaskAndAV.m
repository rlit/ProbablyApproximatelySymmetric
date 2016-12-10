function [Img,roiMask,bounds,steps,ATV,times] = GetVolumeMaskAndAV(K,I1,params,I1_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
workOnDTs = K>0;
tStart = tic;
if workOnDTs
    % create full DT volume
    I1dt_full = double(bwdist(bwperim(I1)));
    I1dt_full(I1==0) = -I1dt_full(I1==0);
    I1dt_full = round(I1dt_full);

    % truncate at k
    I1dt = I1dt_full;
    clear I1dt_full
    I1dt(abs(I1dt)>K) = K * sign(I1dt(abs(I1dt)>K));

    % normalize values to [0,1]
    I1dt = I1dt - min(I1dt(:));
    I1dt = I1dt / max(I1dt(:));
    Img = I1dt;
    clear I1dt
else
    Img = I1;
end
times.tsdf = toc(tStart);
%% roi mask
switch params.roiMaskType
    case 'shpere'
        idx = 1:numel(I1);
        [Ys,Xs,Zs] = ind2sub(size(I1),idx');
        maxDim = size(I1,1);
        rad = (maxDim-1)/2;
        cntr = [rad rad rad]+1;
        diffs = bsxfun(@minus,[Ys,Xs,Zs],cntr)';
        dists = sqrt(sum(diffs.^2));
        dists = reshape(dists,size(I1));
        roiMask = dists < rad;
    case 'interior'
        roiMask = Img>0;
    case 'interiorAndSphereIntersect'
        idx = 1:numel(I1);
        [Ys,Xs,Zs] = ind2sub(size(I1),idx');
        clear idx
        maxDim = size(I1,1);
        rad = (maxDim-1)/2;
        cntr = [rad rad rad]+1;
        diffs = bsxfun(@minus,[Ys,Xs,Zs],cntr)';
        dists = sqrt(sum(diffs.^2));
        dists = reshape(dists,size(I1));
        roiMask = dists < rad;
        roiMask = roiMask & (Img>0); % the intersection
    otherwise
        error('bad roiMakType');
end

%% determine search limits and initial grid
tStart = tic;
[bounds,steps,ATV] = DetermineGridBoundsAndSteps(Img,I1_radius,params.delta,roiMask);
times.tv = toc(tStart);

return
