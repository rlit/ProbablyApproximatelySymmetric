function [bounds,steps,averageVariation] = DetermineGridBoundsAndSteps(I1,I1_radius,delta,roiMask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,w,d] = size(I1);

r2x = 0.5*(w-1);
r2y = 0.5*(h-1);
r2z = 0.5*(d-1);
r1x = 0.5*(w-1);
r1y = 0.5*(h-1);
r1z = 0.5*(d-1);


%% search limits - add reflections!!!
c = 1; % i.e. no scale

minTx = -(r2x-r1x/c);
maxTx = r2x-r1x/c;
minTy = -(r2y-r1y/c);
maxTy = r2y-r1y/c;
minTz = -(r2z-r1z/c);
maxTz = r2z-r1z/c;
minR = -pi;
maxR = pi;
minS = 1/c;
maxS = c;

%% parametrize the initial grid
% roiType = 'sphere';

% by morphology
Imin = imerode(I1,ones([3,3,3])); % the mininum in N
Imax = imdilate(I1,ones([3,3,3])); % the mininum in N
variation = max(I1-Imin,Imax-I1);


roiTV = variation(roiMask);
averageVariation = sum(roiTV(:))/nnz(roiMask);
roiRadius = I1_radius;


[bounds,steps] = GenerateGrid3D(delta,minTx,maxTx,minTy,maxTy,minTz,maxTz,minR,maxR,minS,maxS,...
    averageVariation,roiRadius);

if (c==1)
    steps.s = 0;
else
    error('TODO - verify scale step');
end

return