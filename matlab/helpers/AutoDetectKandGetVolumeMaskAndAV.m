function [Img,roiMask,bounds,steps,ATV,targetK,times] = AutoDetectKandGetVolumeMaskAndAV(I1,params,I1_radius,KvsATVstats)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

targetAV = 3/I1_radius;

Ks      = KvsATVstats.Ks;
avgATVs = KvsATVstats.avgATVs;
stdATVs = KvsATVstats.stdATVs;


% find a good representative query K
[~,queryIndex] = min(abs(targetAV-avgATVs));
targetK = Ks(queryIndex)*I1_radius/100;

% make the query
[Img,roiMask,bounds,steps,queryATV,times] = GetVolumeMaskAndAV(targetK,I1,params,I1_radius);



if (abs(queryATV-targetAV)>0.05)
    prevQueryATV = queryATV;
    % find percentile w.r.t. learned distribution
    stdsAway = (queryATV-avgATVs(queryIndex))/stdATVs(queryIndex);

    beginInd = max(1,queryIndex-2);
    endInd = min(length(Ks),queryIndex+2);
    range = beginInd:endInd;
    standardizedGraph = avgATVs(range) + stdsAway*stdATVs(range);

    x = Ks(range);
    y = standardizedGraph;
    p = polyfit(x,y,2);
    xs = Ks(range(1)):Ks(range(end));
    f = polyval(p,xs);
    % figure; plot(x,y,xs,f,'o-');

    [minVal,minInd] = min(abs(targetAV-f));
    K = xs(minInd);

    targetK = K*I1_radius/100; % Ks(queryIndex-1) + minInd-1;
    estimatedATV = f(minInd);

    % make the query
    [Img,roiMask,bounds,steps,queryATV] = GetVolumeMaskAndAV(targetK,I1,params,I1_radius);

    %% debug figure - a figure that tells the whole story
    if (0)
        expandedRange = max(1,queryIndex-2):min(queryIndex+2,length(Ks));
        ex_x = Ks(expandedRange);
        figure; hold on;
        errorbar(ex_x,avgATVs(expandedRange),stdATVs(expandedRange));
        plot([ex_x(1),ex_x(end)],targetAV*[1,1],':k');
        plot(x,y,'dm');
        plot(xs,f,'-.g');
        plot(xs(minInd),f(minInd),'*m')
        plot(targetK,queryATV,'*r');
    end
end

ATV = queryATV;

% standardised
targetK = targetK*100/I1_radius; % minInd*I1_radius/100; % Ks(queryIndex-1) + minInd-1;


return
