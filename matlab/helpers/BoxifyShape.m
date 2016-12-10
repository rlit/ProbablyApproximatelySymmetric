function [Iout,radius] = BoxifyShape(Iin)

% clipping
idx = find(Iin);
[Ys,Xs,Zs] = ind2sub(size(Iin),idx);
I_clip = Iin(min(Ys):max(Ys) , min(Xs):max(Xs) , min(Zs):max(Zs));

% COG % radius
idx = find(I_clip);
[Ys,Xs,Zs] = ind2sub(size(I_clip),idx);
COG = round(mean([Ys,Xs,Zs]));
diffs = bsxfun(@minus,[Ys,Xs,Zs],COG)';
dists = sqrt(sum(diffs.^2));
radius = ceil(max(dists));

sz = size(I_clip);
Iout = zeros([2*radius+1,2*radius+1,2*radius+1]);
rngFun = @(i)(radius+1-(COG(i)-1)) : (radius+1-(COG(i)-1) + sz(i)-1);
yRng = rngFun(1);
xRng = rngFun(2);
zRng = rngFun(3);
Iout(yRng,xRng,zRng) = I_clip;
