function p = GetParams_Sym(p)

if nargin<1 || isempty(p) || ~isstruct(p)
    p = struct();
end

p = AddField(p,'returnIfResultFileExist',false);

p = AddField(p,'plot_symmetries',true);
p = AddField(p,'createBlenderFile',false);
p = AddField(p,'createEpsFile',false);

p = AddField(p,'defineYourOwnShape',true);
p = AddField(p,'calcFullError',true);
p = AddField(p,'findRollMultiplicities',false);

% precisions
p = AddField(p,'epsilon',0.12);
p = AddField(p,'delta',0.26);

% verbosity
p = AddField(p,'silent',true);
p = AddField(p,'quiet',true);
p = AddField(p,'debug',false);
p = AddField(p,'doVisualization',true);

p = AddField(p,'roiMaskType','interiorAndSphereIntersect'); % not 'interior' any more
% p = AddField(p,'roiMaskType','interior');
% p = AddField(p,'roiMaskType','shpere');

p = AddField(p,'K',0);% use 0 for Binary
p = AddField(p,'maxNumberOfBestSymmetries',4);
p = AddField(p,'constantNumberOfBestSymmetries',false);


p = AddField(p,'overrideWithBinaryImage',false);
p = AddField(p,'allowReflections',true);
p = AddField(p,'onlyReflections',false);
p = AddField(p,'minAngleBetweenSymmetries',18*pi/180);
p = AddField(p,'goodSymmetryThreshold',0.05);
p = AddField(p,'autoDetectK',true);
p = AddField(p,'detectedSymmetries',[]);

p = AddField(p,'breakSearchByThreshold',false);
p = AddField(p,'expandDetectedSymmetry',true);

% B&B
p = AddField(p,'bnbMaxLevel',inf);
p = AddField(p,'allowRestart',true);


% legacy parameters
p = AddField(p,'optA',[]); 
p = AddField(p,'transVec',[]); 
p = AddField(p,'photometricInvariance',false); 
p = AddField(p,'forPaper',false); 
p = AddField(p,'saveWSdir',''); 
p = AddField(p,'OutputDir',''); 


function p = AddField(p,name,val)
if ~isfield(p,name)
    p.(name) = val;
end
