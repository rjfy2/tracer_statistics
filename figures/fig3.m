%read in csv on noise by Allen (supplied with main text)
data=importdata('Data behind SF7.xlsx');

load('cortInOutCCO50');
%load full 
load('rawCortInOut');%load the full set of experiments, not just those into cortex
Out=OutC;
%take the full In
load('50InOut','In','experimentsIds')

load('atlasData');
%differentiate between hemispheres
dim3=size(atlas);

%also read the 295region brain, which corresponds to the manual analysis
load('exps.mat')
%differentiate between hemispheres
Prcl=[Prcl_ID(:,1);Prcl_ID(:,2)+10000];
vol=[Prcl_Vol(:,1);Prcl_Vol(:,2)];
factor=.05^3;%convert to mm3
vols=vols*factor;
%% plot the atlases
ADD=max(atlas(:))/2;
%%
ids = data.data.manualAnalysis(1,:);
ids(isnan(ids))=[];
N=length(ids);


%match the 20 Allen noise experiments with their experimental data: experimentsIds on ids
experimentsIds=str2num(cell2mat(experimentsIds));
expids=zeros(N,1);
expOut=zeros(N,nreg);
expnormIn=zeros(N,nreg);
expIn=zeros(N,1);
for i=1:N
    id=find(experimentsIds==ids(i));
    if ~isempty(id)
        expids(i)=id;
            % get the Out values for these 20 experiments
         expOut(i,:)=Out(id,:)*factor;
         %total In value
         expIn(i)=sum(In(id,:))*factor;
         %relative In value
         expnormIn(i,:)=InC(id,:)/sum(In(id,:));
    end

end
nE=length(expids);

%% reshape the data points
valsAllen=data.data.manualAnalysis(5:end,:);
valsAllen(isnan(valsAllen))=[];
valsAllen=(reshape(valsAllen,size(data.data.manualAnalysis,1)-4,N));


%% now match the rows to the regions
N=size(atlasData,1);
regids=cell(N,1);
done=false(N,1);
for i=1:N %regids(i) gives the location in the sheet of the first of the two that correspond to region i of the parcellation
    id=find(strcmp(strtrim(atlasData(i,3)),strtrim(data.textdata.manualAnalysis(:,2))),1);
    if ~isempty(id)
        regids{i}=id;
        done(i)=true;
    else
        disp(id)
        if(strcmp(strtrim(atlasData(i,3)),'SSp-un'))%misspelled
            regids{i}=find(strcmp('SSp-uk',strtrim(data.textdata.manualAnalysis(:,2))),1);
            done(i)=true;
        end
    end
end
done=find(done);
%% match each region of interest with the manual assigned ones (match to the new 2015 parcellation)

for i=1:295
    regids=assignRegsUp(done(i),atlasData,regids,allIds,allParents);
    regids=assignRegsDown(done(i),atlasData,regids,allIds,allParents);
end
%% find the regions used in the current atlas
regUsed=false(size(atlasData,1),1);
for i=1:nreg
    regUsed(allIds==atlasIds(i))=true;
end
%% go through the experiments. Get the values for each region, and the manual assignments
PTLpID=find(strcmp(atlasData(:,3),'PTLp'));
PosNeg=nan(nE,nreg);
for ctr=1:nE
    expID=expids(ctr);
    if(expID~=0)
        disp(ctr/nE)
        
        for hemi=0:1
            for i=1:nreg/2
                regionID=find(allIds==atlasIds(i));
                %now check if all its subregions are negative
                idsM=regids{regionID};
                if i>=nreg/2-1 %the two PTLp regions
                    if ~isempty(idsM)
                        noafdsf%error
                    end
                    idsM=PTLpID;
                end
                if isempty(idsM)
                    regions(i)
                end
                hits=NaN(length(idsM),1);
                for jj=1:length(idsM)
                    j=idsM(jj);
                    judgement=(data.textdata.manualAnalysis(j-hemi+1,2*ctr+3));
                   if length(judgement)>1
                       'more than one region assigned' %should do something with this
                   end
                    if length(judgement{1})<8
                        hits(jj)=-1;
                    else
                        if strcmpi(judgement{1}(1:8),'TRUE NEG')
                            hits(jj)=0;%false signal according to Allen
                        elseif strcmpi(judgement{1}(1:8),'TRUE POS')
                            hits(jj)=1;%true signal according to Allen
                        else
                            hits(jj)=-1;
                            judgement{1}
                        end
                    end
                end
                
                %false signal if all are false signal (or unknown)
                PosNeg(ctr,i+hemi*(nreg/2))=1-all(hits<=0);
                if all(hits==-1)
                    PosNeg(ctr,i+hemi*(nreg/2))=-2;
                    if ~isempty(idsM)
                    hits
                    length(idsM)
                    judgement{1}
                    end
                end
                
            end
        end
    end
end
xs=nan(nreg/2,nE*2);%signal
trueSignal=nan(nreg/2,nE*2);%manually assigned: true connection or not?
meas=nan(6,nreg/2,nE*2);%measure different aspects
Dens=bsxfun(@times,expOut,1./vols);
normVol=bsxfun(@times,expOut,1./expIn);
normDens=bsxfun(@times,normVol,1./vols);
Fraction=bsxfun(@times,expOut,1./sum(expOut,2));

%combine the two hemispheres

for i_reg=1:nreg/2
    x=log10([Dens(:,i_reg);Dens(:,i_reg+nreg/2)]);
    vals=[PosNeg(:,i_reg+nreg/2);PosNeg(:,i_reg)]
    
    use=x>-Inf&vals>=0;
    

    % logistic regression
    if(sum(vals(use))==0)%||i_reg==33)%all 0s, just assume one after largest value measured
        %%33, SSs, is strange: only one positive for all values. Assume
        %%almost always segmentation error
        
        %all zeros. Assume highest non-measured is actually positive
        id=x==max(x);
        use(id)=1;
        vals(id)=1;
        i_reg

     elseif(sum(vals(use))==sum(use))%all ones
         i_reg
        
    end
        
% save values to draw from during simulations
xs(i_reg,:)=x;
meas(1,i_reg,:)=x;
meas(2,i_reg,:)=log10([expOut(:,i_reg);expOut(:,i_reg+nreg/2)]);%total volume

meas(3,i_reg,:)=log10([normVol(:,i_reg);normVol(:,i_reg+nreg/2)]);%volume signal/volume inj

meas(4,i_reg,:)=log10([normDens(:,i_reg);normDens(:,i_reg+nreg/2)]);%norm density

 meas(5,i_reg,:)=log10([Fraction(:,i_reg);Fraction(:,i_reg+nreg/2)]);

trueSignal(i_reg,:)=vals;
end
%save('noisePerRegionCort','trueSignal','xs')

%%
noo=5;
for i=1:noo
    MIN=min(meas(i,meas(i,:)>-Inf));
    MIN=-11;%manually assigned
    meas(i,meas(i,:)<MIN-1)=MIN-1;
end
%%
%U testing


ids=[1 2 4];%dens, Vol, norm dens
noo=length(ids);
for ii=1:noo
    i=ids(ii);

%Mann Whitney U statistic: number of wins
	pos=meas(i,trueSignal(:)==1);
	neg=meas(i,trueSignal(:)==0);
	A=repmat(pos,length(neg),1);
	B=repmat(neg',1,length(pos));
	U=mean(A(:)>B(:))
end


%%
    
names={'Log density','Signal volume','Normalised volume','Log normalised density','Fraction'};


ids=[1 4];
noo=length(ids);
for ii=1:noo
    i=ids(ii);
    subplot(2,noo,ii)
	xx=min(meas(i,:)):.1:max(meas(i,meas(i,:)<Inf));
	plot(xx,ksdensity(meas(i,trueSignal(:)==1),xx),'k','linewidth',2)
	hold all

	plot(xx,ksdensity(meas(i,trueSignal(:)<0),xx),'--k')

	plot(xx,ksdensity(meas(i,trueSignal(:)==0),xx),'color',[.5 .5 .5],'linewidth',2)

	xlim([min(xx),max(xx)]);
	ylim([ 0 .5])
	if ii==1
		ylabel('Probability density')
	end

	hold off
	%also show for one particular region
	subplot(2,noo,ii+noo)
	a=15;
	plot(xx,ksdensity(squeeze(meas(i,a,trueSignal(a,:)==1)),xx),'k','linewidth',2)
	hold all
	if any(trueSignal<0)
		plot(xx,ksdensity(squeeze(meas(i,a,trueSignal(a,:)<0)),xx),'--k')
	end
	plot(xx,ksdensity(squeeze(meas(i,a,trueSignal(a,:)==0)),xx),'color',[.5 .5 .5],'linewidth',2)

	scatter(squeeze(meas(i,a,trueSignal(a,:)==1)),ones(sum(trueSignal(a,:)==1),1)*.05,30,[0 0 0],'f')
	scatter(squeeze(meas(i,a,~trueSignal(a,:))),zeros(sum(trueSignal(a,:)==0),1)*.05,30,[0.5 0.5 0.5],'f')
	hold off

	pos=squeeze(meas(i,a,trueSignal(a,:)==1));
	neg=squeeze(meas(i,a,trueSignal(a,:)==0));
	A=repmat(pos',length(neg),1);
	B=repmat(neg,1,length(pos));
	xlabel(names(i));
	xlim([min(xx),max(xx)]);
	ylim([ 0 .4])
	U=mean(A(:)>B(:))
	if i==1
    		ylabel('Probability density')
	end
end
set(gcf,'PaperUnits','Centimeters','PaperPosition',[0 0 9 9]);
print('fig3a_d.pdf','-dpdf','-r300')

%% are hemispheres related?

subplot 111
vals=xs;
 vals(xs<-20)=nan;
vals(trueSignal~=1)=NaN;
connR=nanmedian(vals(:,1:end/2),2);
connL=nanmedian(vals(:,end/2+1:end),2);

[r,p]=corr(connR(~isnan(connR)),connL(~isnan(connR)))
mean(connR(~isnan(connR))>connL(~isnan(connR)))
x=sum(connR(~isnan(connR))>connL(~isnan(connR)));%test
N=sum(~isnan(connR));
pl=binocdf(x,N,.5)

vals=xs;
 vals(xs<-20)=nan;
vals(trueSignal~=0)=NaN;
noiseR=nanmean(vals(:,1:end/2),2);
noiseL=nanmean(vals(:,end/2+1:end),2);

[r,p]=corr(noiseR(~isnan(noiseR)),noiseL(~isnan(noiseR)))
mean(noiseR(~isnan(noiseR))>noiseL(~isnan(noiseR)))%one larger than other?
x=sum(noiseR(~isnan(noiseR))>noiseL(~isnan(noiseR)));%test
N=sum(~isnan(noiseR));
pl=binocdf(x,N,.5)
scatter([connR;noiseR],[connL;noiseL],15,[ones(length(connR),1);repmat(2,length(noiseR),1)],'f')
colormap([0 0 0;.5 .5 .5])
hold on
plot([-7 -1],[-7 -1],'k');
xlim([-7 -1]);
ylim([-7 -1]);
hold off
ylabel({'Median log density','(right hemisphere)'})
xlabel({'Median log density','(left hemisphere)'})
w=6;
h=w-1;
set(gcf,'PaperUnits','Centimeters','PaperPosition',[0 0 w h],'Papersize',[w h]);
print(‘fig3e.pdf’,’-dpdf')
%% are the different regions different? Do anova on all distributions, and pairwise t-tests on contralaterals?
noms=zeros(43*40,1);
no=0;
assigned=zeros(43*40,1);
for i=1:43
    ids=trueSignal(i,1:20)==0;
    noms(no+(1:sum(ids)))=xs(i,ids);
    assigned(no+(1:sum(ids)))=i;
    no=no+sum(ids);
    ids=trueSignal(i,21:40)==0;
    noms(no+(1:sum(ids)))=xs(i,20+find(ids));
    assigned(no+(1:sum(ids)))=i+nreg/2;
    no=no+sum(ids);
end
noms=noms(1:no);
assigned=assigned(1:no);
noms(noms<-10)=-10;

anova1(noms,mod(assigned,nreg/2))
kruskalwallis(noms,mod(assigned,nreg/2))
%plot values in ascending order
medians=zeros(nreg/2,1);
for i=1:nreg/2
    medians(i)=median(noms(assigned==i|assigned==i+nreg/2));
end
[~,ord]=sort(medians);
a=assigned;
a(a>nreg/2)=a(a>nreg/2)-nreg/2;
b=a;
for i=1:nreg/2
    b(a==i)=find(ord==i);
end
scatter(b,noms,3,[0 0 0],'f')
ylim([-10 -2])
set(gca,'XTick',[-9 -6 -3])
xlabel('Region')
ylabel('Signal density')
xlim([1 nreg/2])
w=4;
h=w;
set(gcf,'PaperUnits','Centimeters','PaperPosition',[0 0 w h],'Papersize',[w h]);
print(‘fig3f’,’-dpdf')

%% testing; do we see differences between regions?
for i=1:nreg/2
    [h,p]=ttest2(noms(assigned==i),noms(assigned==i+nreg/2));
    [pu]=ranksum(noms(assigned==i),noms(assigned==i+nreg/2));
    if(min(p,pu)<.1)
        [p pu ]
        regions{i}
    end
end
%% lack of differences may be lack of power; test for randomly selected region pairs
N=1000;
nos=zeros(N,1);
for n=1:N
    n
for i=1:nreg/2
    j=randsample(43,1)+43;
    [pu]=ranksum(noms(assigned==i),noms(assigned==j));
    if(pu<.1)
        nos(n)=nos(n)+1;
    end
end
end
mean(nos<1)
median(nos)
min(nos)
