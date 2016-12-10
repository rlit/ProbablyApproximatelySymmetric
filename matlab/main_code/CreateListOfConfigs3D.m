function [configs,gridSize] = CreateListOfConfigs3D(bounds,steps,params,onlyGridSize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allowReflections = params.allowReflections;
onlyReflections = params.onlyReflections;

if ~allowReflections
    onlyReflections = false;
end


tx_steps = bounds.tx(1) : steps.tx : bounds.tx(2);
ty_steps = bounds.ty(1) : steps.ty : bounds.ty(2);
tz_steps = bounds.tz(1) : steps.tz : bounds.tz(2);
r_steps = bounds.r(1) : steps.r : bounds.r(2);
s_steps = bounds.s(1) : steps.s : bounds.s(2);
if (steps.s == 0)
    s_steps = bounds.s(1);
end
if (bounds.tx(2)-tx_steps(end) > 0.5*steps.tx)
    tx_steps = [tx_steps tx_steps(end)+steps.tx];
end
if (bounds.ty(2)-ty_steps(end) > 0.5*steps.ty)
    ty_steps = [ty_steps ty_steps(end)+steps.ty];
end
if (bounds.tz(2)-tz_steps(end) > 0.5*steps.tz)
    tz_steps = [tz_steps tz_steps(end)+steps.tz];
end
if (bounds.r(2)-r_steps(end) > steps.r) % this one is cyclic (so don't take *0.5)
    r_steps = [r_steps r_steps(end)+steps.r];
end
if (bounds.s(2)-s_steps(end) > 0.5*steps.s)
    s_steps = [s_steps s_steps(end)+steps.s];
end

ntx_steps = length(tx_steps);
nty_steps = length(ty_steps);
ntz_steps = length(tz_steps);
ns_steps = length(s_steps);
nr_steps = length(r_steps);

% NR2_steps = nr_steps;
% if (bounds.r(2)-bounds.r(1) == 2*pi)
%         quartile1_r_steps = r_steps(r_steps < -pi/2 + steps.r/2);
%         NR2_steps = length(quartile1_r_steps);
% end

% grid size calculation
L = steps.r; % rotation step size
n = ceil(pi/L); % number of patches around half the equator
N_spherePoints = ceil((2*n)^2/pi); % total number of patches on sphere
% latitudes =  (0:n)*L;
if (allowReflections && ~onlyReflections)
    numRolls = round(nr_steps/2) + 1;
elseif (allowReflections && onlyReflections)
    numRolls = 1;
else
    numRolls = round(nr_steps/2);
end

gridSize = ntx_steps*nty_steps*ntz_steps*ns_steps*N_spherePoints*numRolls;
% keyboard
if (exist('onlyGridSize','var') && onlyGridSize)
    configs = [];
    return
end

if (gridSize > 26000000)
    error('oversized gridSize');
end




%%

configs = CreateList(...
    ntx_steps,nty_steps,ntz_steps,...
    ns_steps,nr_steps,...
    tx_steps,ty_steps,tz_steps,...
    r_steps,s_steps,L,gridSize,...
    allowReflections,onlyReflections);




% Matlab version
function configs = CreateList(...
    ntx_steps,nty_steps,ntz_steps,...
    ns_steps,nr_steps,...
    tx_steps,ty_steps,tz_steps,...
    r_steps,s_steps,L,gridSize,...
    allowReflections,onlyReflections)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method follows that Paper "Novel Techniques for Creating Nearly Uniform Star Catalog" by Mortari and others

n = ceil(pi/L); % number of patches around half the equator
% N_spherePoints = ceil(n^2/pi); % total number of patches on sphere
latitudes =  (0:n)*L;

configs = zeros(gridSize,7);

gridInd = 0;
for tx_ind = 1 : ntx_steps
    tx = tx_steps(tx_ind);
    %fprintf('step %d out of %d\n',tx_ind,ntx_steps);
    for ty_ind = 1 : nty_steps
        ty = ty_steps(ty_ind);
        for tz_ind = 1 : ntz_steps
            tz = tz_steps(tz_ind);
            for s_ind =  1 : ns_steps
                s = s_steps(s_ind);
                %%  computing equidistance points on a sphere
                for latIndex = 1 : n % loop on latitudes
                    if (latIndex==1) % north pole
                        numPatchesAtInclination = 1;
                    else
                        totPatchAreaAtInclination = 2*pi*(cos(latitudes(latIndex-1)) - cos(latitudes(latIndex)));
                        numPatchesAtInclination = round(totPatchAreaAtInclination/L^2);
                    end
                    for j = 1 : numPatchesAtInclination % loop on longitudes
                        longitude = 2*pi*(j-1)/numPatchesAtInclination;

                        if ~onlyReflections
                            for r_ind = 1 : n % rolls in [-pi,0]
                                roll = r_steps(r_ind);
                                gridInd = gridInd + 1;
                                configs(gridInd,:) = [tx,ty,tz,s,latitudes(latIndex),longitude,roll];
                            end
                        end
                        %                                                 px = sin(latitudes(i))*cos(longitude);
                        %                                                 py = sin(latitudes(i))*sin(longitude);
                        %                                                 pz = cos(latitudes(i));
                        if(allowReflections && (latitudes(latIndex) < pi/2))
                            gridInd = gridInd + 1;
                            
                            configs(gridInd,:) = [tx,ty,tz,-s,latitudes(latIndex),longitude,0];
                        end
                    end
                end
                
            end
            
        end
    end
end

% keyboard
configs = configs(1:gridInd,:);

return


