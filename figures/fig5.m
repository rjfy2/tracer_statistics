%load MCMC results
load('1901_runs');%based on 5 runs
%the connections for which information is available (i.e. there is at least one experiment that gives information on this)
load('knowable')
%there connections should have the prior distribution (uniform in the main text). Set them to nan to not confuse us here
medMu(~knowable)=nan;
MUS(:,~knowable)=nan;


%for entire posterior
mm=10.^MUS(ran,:,:);
%convert to expected total signal
%some tricks to extend vols
v=repmat(vols,[1,1,31]);
v=permute(v,[1,3,2]);
mm=bsxfun(@times,mm,v);

musF=log10(bsxfun(@times,mm,1./nansum(mm,3)));%convert to percentage weights
musF(:,~knowable)=nan;


%% do similar for macaque
%Kennedy
cd('/work/imagingL/autism_graphtheory/Tracers/macaque Kennedy')
[macIn, macOut, FLN, source, unS, target, macNames, data, PNASdata,total] = loadKennedyData();

expSets={1:5,6:8,9:10,[34,37:39]};
 load('volume','macVols')
%% 
%for each source region, find the best selectivity for each target
%(function of injection volume and target region noise)
normIn=bsxfun(@times,In,1./sum(In,2));
injVols=sum(In,2);
outDens=bsxfun(@times,Out,1./vols);
load('noisePerRegionCort')
    
    medNoise=zeros(size(Out));
     %symmetry:
    xs=[xs;xs];
    trueSignal=[trueSignal;trueSignal];
  
    factor = .05^3;%convert to mm3 - already done for other variable
    xs = xs - log10(factor);
    %find the median density of noise distribution for each region
        for j=1:size(Out,2);
            %look at when 50% above noise
            medNoise(:,j)=quantile(xs(j,trueSignal(j,:)==0),.5);
            
        end

%for each experiment, each target region, what mu would be needed to  exceed noise?
%find by dividing the density needed by the injection volume
     noiseThr=log10(bsxfun(@times,10.^medNoise,1./injVols));
     %assuming total signal out (not in target region) fixed; how much of this signal do we need
     %to exceed noise for each region?
     OutnotIn=Out;
     OutnotIn(:,in)=OutnotIn(:,in)-In;
     percNoiseThr=log10(bsxfun(@times,bsxfun(@times,10.^medNoise,vols),1./sum(Out,2)));
 
muNeededN=nan(Ninjected,Nregions);
percNeeded=muNeededN;
for i=1:Ninjected

    %the medMu needed for this region to exceed noise
    needed=log10(bsxfun(@times,10.^medNoise,1./In(:,i)));
    muNeededN(i,:)=min(needed)';
    

    %the experiments this region was injected in
     exps=normIn(:,i)>.5;
        
   %for each region, take the minimum needed
    percNeeded(i,:)=min(percNoiseThr(exps,:),[],1);
   
end



  %%  density with noise
namesDens={'log NCD','% network density'};
upp=0;
bins=-13:.1:upp;
nb=length(bins);
dens=nan(length(bins),1);
densIpsi=nan(length(bins),1);
densContra=densIpsi;

distMu=musF;
Dens=nan(length(bins),size(distMu,1));%probabilistically, look at the full posterior
DensIpsi=Dens;

        
        %this is the masking by other connections 
        percNow=10.^MUNN(ran,:,:).*permute(repmat(vols,[size(distMu,1),1, 31]),[1,3,2]);
        percNow=log10(bsxfun(@times,percNow,1./sum(percNow,3)));
        
        need=max(permute(repmat(percNeeded,[1,1,size(distMu,1)]),[3,1,2]),percNow);

        

musFI=distMu(:,:,end/2+1:end);
DensContra=Dens;
musFC=distMu(:,:,1:end/2);
ni=nan(length(bins),2);
for i=1:length(bins)

    ids=need(:,:,end/2+1:end)<bins(i);
    ni(i,1)=sum(ids(:));%no of connections of sufficient specificity
    if any(ids(:))
        mm=musFI;
        mm(~ids)=nan;
        DensIpsi(i,:)=nansum(nansum(mm>bins(i),3),2)./sum(sum(ids>0,3),2); 
    end
    ids=need(:,:,1:end/2)<bins(i);
    ni(i,2)=sum(ids(:));
    if any(ids(:))
        mm=musFC;
        mm(~ids)=nan;
        
        DensContra(i,:)=nansum(nansum(mm>bins(i),3),2)./sum(sum(ids>0,3),2); 
    end
end

%%
subplot 211
plot(bins,DensContra*100,'color',[.8 .8 1])
hold on
plot(bins,DensIpsi*100,'color',[.4 .4 .8])

bndsMouse=quantile(need(:),[.1 .5 .9]);
fMouse=ksdensity(need(need(:)>-30),bins,'width',.2);

oneNeuronMouse=log10(1./(4.2e6/43)); %from HH: 5e6 total minus PIR and ENT
xlabel('log fractional weight')
xlim([bndsMouse(1)-1 upp])
ylim([0 80])
ylabel(namesDens{2})


% get CI
col=find(abs(bndsMouse(2)-bins)==min(abs(bndsMouse(2)-bins)));
contraCI=quantile(DensContra(col,:),[.025,.5,.975])
ipsiCI=quantile(DensIpsi(col,:),[.025,.5,.975])

plot(bins,fMouse*100,'color',[.4 .4 .6],'linewidth',2)
plot([oneNeuronMouse,oneNeuronMouse],[0 100],'--','color',[.4 .4 .6],'linewidth',3)

hold off
display('size vols to In')
mean(vols)
mean(sum(In*factor,2))
sum(vols(in))
sum(In(:)*factor)
%% look at our estimates from Kennedy data
cd('/work/imagingL/autism_graphtheory/Tracers/2015_data/logliknormflip')
load('2910_50_Mean_Kenn.mat')

%for entire posterior
mm=10.^Mus(end/5:end,:,:);
%convert to expected total signal
%some tricks to extend vols
v=repmat(macVols,[1,1,29]);
v=permute(v,[2,3,1]);
mm=bsxfun(@times,mm,v);

%have to permute due to MATLAB oddness..
musFmac=permute(log10(bsxfun(@times,permute(mm,[3 1 2]),1./nansum(permute(mm,[3,1,2]),1))),[2 3 1]);% is this the same as estimating percentages straight away?

subplot 212
oneNeuron=log10(1./(1.5e9/91)); %from Collins: 1.3billion in right hemisphere (not sure which regions included)
%HH 2007 PNAS: 1.7e9 (whole cortex? hemisphere?)

plot(bins,DensIpsi*100,'color',[.4 .4 .8])
 hold on
macMin=log10(1/(4e8/2/91));%if every neuron would send axon to other cortical region
estMac=nan(29,91);

estTotal=nan(29,1);
mI=macIn(:,sum(macIn)>0);
for i=1:29
    estMac(i,:)=mean(macOut(mI(:,i)==1,:),1);
    estTotal(i)=max(total(mI(:,i)==1));
end
estMac=log10(estMac);

nb=length(bins);
macdensR=nan(size(musFmac,1),nb);
for i=1:nb
    %density at this level
    
     exps=log10(1./estTotal)<bins(i);%only take experiments with enough sensitivity at this level
    if any(exps(:))
    	macdensR(:,i)=nansum(nansum(musFmac(:,exps,:)>bins(i),3),2)./sum(sum(~isnan(musFmac(:,exps,:)),3),2); 
    end

end

hold on
plot(bins,macdensR*100,'color',[.6 .8 .6])

% get CI
bnds=quantile(log10(1./estTotal),[.1 .5 .9]);

f=ksdensity(log10(1./estTotal),bins,'width',.15);
plot(bins,f/3*100,'color',[.4 .6 .4],'linewidth',2)


plot([oneNeuron,oneNeuron],[0 100],'--','color',[.4 .6 .4],'linewidth',3)
    

xlabel('log fractional weight')
xlim([bndsMouse(1)-1 upp])
ylim([0 80])
ylabel(namesDens{2})
hold off
col=find(abs(bnds(2)-bins)==min(abs(bnds(2)-bins)));
macCI=quantile(macdensR(:,col),[.025,.5,.975])
%%
wi=9;
hi=9;
set(gcf,'PaperUnits','Centimeters','PaperPosition',[0 0 wi hi],'PaperSize',[wi hi])
print('fig5.pdf','-dpdf')
%% get upper and lower bounds across density for Kenn based on the above

nb=length(bins);
macdens=nan(nb,1);
macdensNaive=macdens;

for i=1:nb
    %density at this level
    exps=log10(1./estTotal)<bins(i);%only take experiments with enough sensitivity at this level
    macdens(i)=sum(sum(estMac(exps,:)>bins(i)))/sum(exps)/90; %check - zit zelf hier niet in?
    macdensNaive(i)=sum(estMac(:)>bins(i))/29/90;
end
hold all
%show the sensitivity threshold
plot(log10(1./[median(estTotal),median(estTotal)]),[0,0.8],'--','color',[.7 .4 .4])
plot(log10(1./[quantile(estTotal(:),0.025),quantile(estTotal(:),0.025)]),[0,0.8],'--','color',[.7 .4 .4])
plot(log10(1./[quantile(estTotal(:),0.975),quantile(estTotal(:),0.975)]),[0,0.8],'--','color',[.7 .4 .4])
hold off
xlabel('log percentage')
ylabel('density')

