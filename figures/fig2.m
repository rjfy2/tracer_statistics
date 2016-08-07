%investigate consistency of injections;
PERC=true;%look at percentages or ncd?
%Kennedy
[macIn, macOut, FLN, source, unS, target, macNames, data, PNASdata,total] = loadKennedyData();

if ~PERC
    macOut=bsxfun(@times,macOut,total);
end
%% plot the repeat measurements for Kenn, and look a bit at consistency
imagesc(log10(macOut([1:10,34,37:end],:)))
expSets={1:5,6:8,9:10,[34,37:39]};
for i=1:4
    %density
    density=mean(sum(macOut(expSets{i},:))>0);
    %density per experiment
    dpe=mean(mean(macOut(expSets{i},:)>0,2));
    %density if we take out those singly found
    string=mean(sum(macOut(expSets{i},:)>0)>1);
    %density if we take only those always found
    string2=mean(sum(macOut(expSets{i},:)>0)==length(expSets{i}));
    [density,(dpe-density)/density,(string-density)/density,(string2-density)/density]
    
    %distribution of no of times found
    tabulate(sum(macOut(expSets{i},:)>0))
end

%%
load('atlasData')

%take our cortical one
load('cortInOutCCO50')

%divide by out volume to obtain fractional weight estimates
if PERC
    cortOut=bsxfun(@times,Out,1./sum(Out,2));
else
    cortOut=bsxfun(@times,Out,1./sum(In,2));
end
cortIn=In;



%only take into cortical regions
cortregs=(ismember(atlasIds,allIds(ismember(allParents,allIds(isC)))));

% and cocomac
load('cocovals')
%select all connections with at least a number of weighted measurements
noWei=sum(ismember(coco,[101,102,103]),3);
[a,b]=find(noWei>10);
cocoConn=nan(50,length(a));
for i=1:length(a)
    ids=find(ismember(coco(a(i),b(i),:),[101,102,103]));
    if length(ids)>50
        ids=ids(1:50);
    end
    cocoConn(1:length(ids),i)=coco(a(i),b(i),ids);
end
cocoConn=(cocoConn-104)*2;


%% get macaque volumes
load('macaque_volume')

cutoff=.9;

nOuts=3;
Outs=cell(nOuts,1);
Ins=cell(nOuts,1);



%%

%arrange experiments
for l=1:2
    if l==1;
        Inn=cortIn;
        Outt=cortOut;
    elseif l==2;
        Inn=macIn;
        Outt=macOut;
    end
    normIn=bsxfun(@times,Inn,1./sum(Inn,2));
    maxE=ceil(max(sum(normIn>cutoff)));
    regs=find(sum(normIn>cutoff)>0);
    Outs{l}=[];
    for j=1:length(regs)
        exps=normIn(:,regs(j))>cutoff;
        Outs{l}=[Outs{l},[Outt(exps,:);nan(maxE-sum(exps),size(Outt,2))]];
    end
end

%all cocomac values

strCoco=coco;
strCoco(strCoco==101.23)=nan;%cant use 'exists' for strength analysis
strCoco(strCoco==100)=nan;%cant use non-existence for strength analysis
strCoco=strCoco-100;
noStr=squeeze(sum(~isnan(strCoco),3));
maxE=max(noStr(:));
Outs{3}=nan(maxE,sum(noStr(:)>0));
s=size(strCoco,1);
for i=1:s
    for j=1:s
        if noStr(i,j)>0
            exps=~isnan(strCoco(i,j,:));
            Outs{3}=[Outs{3},10.^[squeeze(strCoco(i,j,exps));nan(maxE-sum(exps),1)]];
        end
    end
end

%%

studyNames={['Allen'],['Markov'],'Cocomac'};
levNames={'Log signal fraction','Variance','Weight distribution'};
% plot regions and their std

for i=1:nOuts
    x=log10(Outs{i});
    limit=min(x(x(:)>-Inf));
    x(x==-Inf)=nan;%remove 0's
    %remove connections with 1 or less measurements after saving means
    x=x(:,sum(~isnan(x))>0);
    ms=nanmean(x);
    x=x(:,sum(~isnan(x))>1);
    [~,ord]=sort(nanmean(x));
    x=x(:,ord);
    subplot(2,nOuts,i)
    
    jitter=rand(size(x))*.5-.25;
    if i<3
        jitter=jitter*0;
    end
    plot(repmat(1:size(x,2),size(x,1),1),x+jitter,'.k','markersize',1)
    nC=sum(sum(~isnan(x))>1)
    xlim([.9 nC])
    ylim([-10 0])
    if ~PERC
        ylim([-5 5])
    end
    if i==3
        ylim([0 4])
    end
    %     title(studyNames{i})
    subplot(2,nOuts,i+nOuts)
    varss=nanvar(x);
    varss(sum(~isnan(x))<=1)=nan;%remove those with only one value available
    plot(varss(sum(~isnan(x))>0)/nanvar(ms),'k');
    ylim([0 2.5])
    xlim([.9 nC])
    if any(varss==0)&i<3
        fout
    end
    
    V=(nanmean(varss)/nanvar(ms))
    
end
for i=1:2
    subplot(2,nOuts,(i-1)*nOuts+1)
    ylabel(levNames{i})
end
for i=1:nOuts
    subplot(2,nOuts,(2-1)*nOuts+i)
    xlabel('Region pairs')
end

set(gcf,'PaperUnits','Centimeters','PaperPosition',[0 0 12 9],'PaperSize',[12 9]);
print('fig2.pdf','-dpdf','-r300')
