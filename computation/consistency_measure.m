%simple measures of consistency for the three datasets

%%

load('cocovals') %these are extracted from the xls available from cocomac.g-node.org (loaded using readCoCoMatrix.m)

MINMEANEXP=1;%all exps for mean (1), or only those used for var (2)?
keepZeros=0;%0 if discarding zeroes and working with log values
belowQuantile=.5;%1 if all connections. lower for only compute this over lower values
%% transform values to what they would look like for tracing
strCoco=coco;
strCoco(strCoco==101.23)=nan;%cant use 'exists' for strength analysis
if ~keepZeros
    strCoco(strCoco==100)=nan;%cant use non-existence for strength analysis
else
    strCoco=10.^(strCoco-100);
end

noStr=squeeze(sum(~isnan(strCoco),3));
% tabulate(strCoco(:))
vars=squeeze(nanvar(permute(strCoco,[3,1,2])));
means=squeeze(nanmean(permute(strCoco,[3,1,2])));

%remove the one-value ones
vars(noStr<2)=nan;
means(noStr<MINMEANEXP)=nan;

if keepZeros
    %     remove the ones with only zeroes
    vars(sum(permute(strCoco,[3,1,2])>1,1)<1)=nan;
    means(sum(permute(strCoco,[3,1,2])>1,1)<1)=nan;
end
hist((vars(:)))
title([nanmedian(vars(:)),nanmean(vars(:))])
nanmean(vars(:))
nanvar(means(:))
V_Coco=nanmean(vars(:))/nanvar(means(:))%this is invariant to any scaling of strCoco above
imagesc(vars/nanvar(means(:)))

[InK, percOutK, FLN, source, unS, target,cocoNames, data, PNASdata, total]=loadKennedyData();
OutK=bsxfun(@times,percOutK,total);
%%

ids=find(sum(InK)>0);
varsK=[];
meansK=[];
assocMeansK=[];
for i=1:length(ids)
    id=ids(i);
    exps=find(InK(:,id)==1);
    if ~keepZeros
        out=log10(OutK(exps,:));
        out(out==-Inf)=nan;
    else
        out = OutK(exps,:);
    end
    
    toUse=(sum(~isnan(out),1)>1)&sum(InK(exps,:))==0;
    if keepZeros
        %     remove the ones with only zeroes
        toUse(sum(out>0)==0)=0;
    end
    if sum(toUse)>0
    	varsK=[varsK,(nanvar(out(:,toUse)))];
        assocMeansK=[assocMeansK,(nanmean(out(:,toUse)))];
    end
    toUse=(sum(~isnan(out),1)>=MINMEANEXP)&sum(InK(exps,:),1)==0;
    if keepZeros
        %     remove the ones with only zeroes
        toUse(sum(out>0)==0)=0;
    end
    if sum(toUse)>0
        meansK=[meansK,nanmean(out(:,toUse),1)];
    end
    
end
if sum(varsK(:)==0)>0
    fout
end

     size(varsK)
     size(meansK)
     cu=quantile(meansK,belowQuantile);
    varsK=varsK(assocMeansK<=cu);
%     meansK=meansK(meansK<=cu);
         size(varsK)
     size(meansK)
nanmean(varsK)
nanvar(meansK)%how spread are the estimates?
V_Kenn=nanmean(varsK)/nanvar(meansK)%V=measurement var over estimates var (lower is better)

%% in the same, simple way, look at Allen data. 
load('cortInOutCCO50')
%%
percOut=bsxfun(@times,Out,1./sum(Out,2));
relOut=bsxfun(@times,Out,1./sum(In,2));
cus=[0.01:.01:1];
meas=nan(length(cus),1);
    cutoff=.5;
    ids=find(sum(normIn>cutoff)>0);
    varsA=[];
    assocMeansA=[];%the means related to the vars above
    meansA=[];
    for i=1:length(ids)
        id=ids(i);
        exps=find(normIn(:,id)>cutoff);
        if ~keepZeros
            out=log10(percOut(exps,:));%use percOut to make comparable to Kennedy macaque data
        else
            out=percOut(exps,:);
        end
        toUse=sum(~isnan(out),1)>1&max(normIn(exps,:),[],1)<cutoff;
        if keepZeros
            %     remove the ones with only zeroes
            toUse(sum(out>0)==0)=0;
        end
        if sum(toUse)>0
        	varsA=[varsA,(nanvar(out(:,toUse)))];
        	assocMeansA=[assocMeansA,(nanmean(out(:,toUse)))];
        end
        toUse=sum(~isnan(out),1)>=MINMEANEXP&max(normIn(exps,:),[],1)<cutoff;
        if keepZeros
            %     remove the ones with only zeroes
            toUse(sum(out>0)==0)=0;
        end
        if sum(toUse)>0
        meansA=[meansA,nanmean(out(:,toUse),1)];
        end
    end
    cu=quantile(meansA,cutoff);
    varsA=varsA(assocMeansA<=cu);
%     meansA=meansA(meansA<=cu);
    %with no identical values, var=0 means there was only one measurement
    % meansA(varsA==0)=nan;%one measurement
    % varsA(varsA==0)=nan;
    %should be no var=0
    if sum(varsA(:)==0)>0
        fout
    end
    % nanmean(varsA)
    % nanvar(meansA)%how spread are the estimates?
    [numel(varsA) numel(meansA)]
    V_Allen=nanmean(varsA)/nanvar(meansA)

