
%find posterior distribution of parameters, with noise as given by Allen manual experiments 
%Assume distrs are lognormal, mixed is still lognormal (errors are
%dependent), std is constant

%In,Out are NinjxNregion matrices, giving the total signal into and out of
%regions for each experiment
%beta=weight for simulated annealing. Set to Inf for MCMC
%cols target regions to look at
%M is number of iterations
%vols gives the volumes of regions
%FPs Nparcelxdx2 gives information on total signal found and whether this was a false
%positive, for each region
%Outids is the id of the Out region, for each of the In regions
%smallest injection fraction we allow. If lower, remove for computational reasons
    SMALLESTIN=.01;

    %% load the data from the new Allen CCO
%%this was made by masking the cortical regions from the downloadable data, and making sure contralateral homologue regions are distinct
    load('/work/imagingL/autism_graphtheory/Tracers/2015_data/cortInOutCCO50')
    
    factor = .05^3;%convert to mm3 before estimation
    Out=Out*factor;
    In=In*factor;
    vols=vols*factor;
        

    outIds=1:size(Out,2);
    %remove 0’s; assume they have below-threshold signal (either from connections or noise). Assign the lowest volume found (corrected for injection and region size
    [a,b]=find(Out==0);
    injVol=sum(In,2);
    %find smallest value measured
    smallest=bsxfun(@times,Out,1./vols);
    MINVALUE=min(nonzeros(smallest));
    
    zeroMeasured=nan(size(Out));
    zeroMeasured(Out==0)=MINVALUE.*vols(b)';%assign low value weighted by injection and region size
    
    nexp=size(In,1);
    volsNonInj=repmat(vols,nexp,1);%volume not injected, per region per experiment. As we assume there is no information for injected regions, this is equal to region volume
    
    
    %compute the pdf+cdf for each Out measurement that it was false
    %positive.
    load('noisePerRegionCort')
    xs=xs+log10(factor);%convert to mm3 (xs are the measured values for the 20 manually curated experiments, on a log-scale)
    pdfNoise=zeros(size(Out)); % initialise
    cdfNoise=pdfNoise; 
    %assume symmetry: contralateral homologue regions have the same noise profile (this is supported by the data)
    xs=[xs;xs];
    trueSignal=[trueSignal;trueSignal]; %trueSignal is whether the signal was deemed to be a real connection
    
    
    DISCARDEXTREMES=1;%reviewer suggested sensitivity analysis on discarding extreme false positives
    if DISCARDEXTREMES
        for i=1:43
            high=max(xs(i,trueSignal(i,:)==0));
            low=min(xs(i,trueSignal(i,:)==1));
            if high>low%swap the two most extreme values (rater error?)
                trueSignal(i,xs(i,:)==high)=1;
                trueSignal(i,xs(i,:)==low)=0;
                trueSignal(i+43,xs(i,:)==high)=1;
                trueSignal(i+43,xs(i,:)==low)=0;
            else
                [i,high,low]
            end
        end
    end
    
    for i=1:size(Out,1);
        for j=1:size(Out,2);
            %take cumulative
            cdfNoise(i,j)=mean(xs(j,trueSignal(j,:)==0)<log10(Out(i,j)./volsNonInj(i,j)));
        end
    end
    
    for i=1:size(Out,2)
        pdfNoise(:,i)=ksdensity(xs(i,trueSignal(i,:)==0),log10(Out(:,i)./volsNonInj(:,i)));
    end
    
    % median of the noise - not used in calculation but in later plots
    medNoiseDensity=zeros(size(Out));
    %find the median density of noise distribution for each region, in
    %units total signal
    for j=1:size(Out,2);
        medNoiseDensity(:,j)=10.^quantile(xs(j,trueSignal(j,:)==0),.5);
    end
   %%
	MACAQUE = 0 % set to one to analyse macaque data
if MACAQUE
    %For the macaque data by Kennedy et al.
    [In, Out, FLN, source, unS, target, names, data, PNASdata,total] = loadKennedyData();
    load('volume');
    SMALLESTIN=.5;%irrelevant
    outIds=1:size(Out,2);
    volsNonInj=repmat(macVols',size(Out,1),1).*(1-In);
    %zeros means smaller than minimum measurable: one neuron. Make it half a
    %neuron
    MINVALUES=repmat(1./total,1,size(Out,2));
    zeroMeasured=nan(size(Out));
    zeroMeasured(Out==0)=MINVALUES(Out==0);
    
    pdfNoise=Out*0;%assume -no- noise. Seems somewhat optimistic…
    cdfNoise=ones(size(Out));
end
    
    %%
    
    beta=Inf%set to lower for simulated annealing
    M=1000000% # iterations



%%
%no of injections
[Ninj] = size(In,1);
%no of regions
Nregions=size(Out,2);


%safety checks
if size(volsNonInj,2)~=Nregions || size(volsNonInj,1)~=Ninj
    disp('not all volumes available!')
end
if size(In,2)~=Nregions
    disp('In and out not same!')
end
%region volumes
vols=max(volsNonInj,[],1);



%for each experiment, record the out regions injected. These are subsequently not
%evaluated in experiments, as we can’t separate signal due to connection or due to direct injection
injOutRegions=zeros(Ninj,Nregions);
for i=1:Ninj
    %remove if In is a significant fraction of Out
    nz=nonzeros(unique(outIds(In(i,:)./Out(i,:)>0.0)));
    injOutRegions(i,nz)=1;
 end
injOutRegions=logical(injOutRegions);

%record which connections can be measured (i.e. are not removed above)
knowable=false(Nregions,Nregions);
for i=1:size(In,1);
    not=(find(In(i,:)>0));
    knowable(not,~ismember(1:Nregions,outIds(not)))=true;
end


%ignore regions that dont have at least x% of an injection somewhere (31 remain)
inEnough=bsxfun(@times,In,1./sum(In,2))>.5;
in=find(sum(inEnough));

Ninjected=length(in);
%remove empty columns from In
In=In(:,in);
outIds=outIds(in);
knowable=knowable(in,:);


%% lower and upper bound for the priors
lowMu=-20;
highMu=5;
lowSigm=0.01;
highSigm=10;
%% initialize
mus=zeros(Ninjected, Nregions);%only have params for injected regions
relIn=bsxfun(@times,In,1./sum(In,2));
ninj=zeros(Ninjected,1);
normOut=bsxfun(@times,Out,1./sum(In,2));

%initialise to sensible values
for i=1:Ninjected
    vals=bsxfun(@times,normOut(relIn(:,i)>0,:),1./vols);
    ninj(i)=size(vals,1);
    for j=1:Nregions
        if(sum(vals(:,j)>0)>0)
            mus(i,j)=log10(mean(vals(vals(:,j)>0,j)));%set the mean to the weighted mean of experiments
        else
            mus(i,j)=lowMu;
        end
        
    end
end

    coef=1;%constant sigma

%hyperparameters for the step size in the MCMC will be tweaked during burnin
    Ccoeff=1/5/10;
C=ones(Ninjected,2,Nregions)/5;



if any(mus<lowMu|mus>highMu)
    priors too small!
end
%% for MCMC
sampFreq=200;
Mus = nan(M/sampFreq,Ninjected,Nregions);%save samples. burnin is discarded later
%save space
Mus=single(Mus);
    Coefs=single(nan(M/sampFreq,1));

lls=zeros(M/sampFreq,1);%for AIC

muNeededNow=single(nan(M/sampFreq,Ninjected,Nregions));


%%
tic
llmax=ones(Ninj,1)*-Inf;%log likelihood for each of the injections
ll=llmax;

%initialise log likelihood
ll=loglikLognormMean(In,Out,...
    mus,coef,pdfNoise,...
    cdfNoise,zeroMeasured,volsNonInj,injOutRegions);

sum(ll)
toc
current=1;
%%
reportTime=5;%how often do you want a printout?
tic
for m=current:M
    current=m;
    time=toc;
    if(mod(m,sampFreq)==0)%sample
        Mus(m/sampFreq,:,:)=mus;
        Coefs(m/sampFreq,:)=coef(:);
        lls(m/sampFreq)=sum(ll);
        
        
        %compute whether we are above threshold
        totMu=permute(10.^repmat(mus,[1,1,Ninj]),[3,1,2]).*...
            repmat(In,[1,1,Nregions]);
        for i=1:Ninjected
            %compute how much the other regions contribute
            muOther=totMu;
            muOther(:,i,:)=0;%set own contribution to zero
            otherContr=squeeze(nansum(muOther,2));
            %1 how much is needed to contribute ate least 50% here?
            %match the total contributino of others
            muNeededOther=bsxfun(@times,otherContr,1./In(:,i));
            %2 how much is needed to be above the median of the noise?
            %match the difference between total signal needed for noise
            %and what is provided by others
            %                 muNeededNoise=bsxfun(@times,(medNoiseDensity-otherContr),1./In(:,i));
            %take the maximum of these two, and then the minimum over
            %all experiments
            %                 muNeededNow(i,:)=log10(min(max(muNeededOther,muNeededNoise)));
            muNeededNow(m/sampFreq,i,:)=log10(min((muNeededOther)));
        end
        

        
        
        
    end
    
    if(time>reportTime)
        %every so often, print a confusing report

        mus([18,19,24],25)
        [log10(sum(In(41,:).*10.^mus(:,25)')),        log10(Out(41,25)/vols(25))]
        [~, noiseComp, connComp, ~] = loglikLognormMean(In,Out,mus,coef,pdfNoise,cdfNoise,zeroMeasured,volsNonInj,injOutRegions);
        tic
        display('mean strength')
        
        [mean(mus(:))]
        display('density (total above noise):')
       nanmean(nanmean(connComp./(noiseComp+connComp)))
        display('lik:')
        sum(ll)
        
        display('acceptance rates:')
        display([acc/(acc+rej),accCoef/(accCoef+rejCoef)])
        display('MCMC parameters:')
        [m/M  mean(mean(C(:,1,:))) mean(mean(C(:,2,:))) mean(mean(Ccoeff))  mean(coef(:))]
        
        if(m>sampFreq)
%plot some information to check on chain
            subplot 321
            plot(lls(1:floor(m/sampFreq)))
            hold all
            %             plot(plls(1:floor(m/sampFreq)))
            hold off
            subplot 323

            plot(squeeze(Mus(1:floor(m/sampFreq),1:min(5,Ninjected),1)))
            subplot 325
            plot(mean(Coefs(1:floor(m/sampFreq),:),2))
            hold all
            plot(max(Coefs(1:floor(m/sampFreq),:),[],2))
            hold off
            
            subplot 322
            imagesc(mus,[-10, highMu]);
            
            subplot 324
            scatter(mus(:),coef,20,knowable(:),'filled')
            subplot 326
            hist(mus(:))

            drawnow
            
        end
    end
    
    acc=0;
    rej=0;
    %update parameters for each injected region
    for i=1:Ninjected
        %         only change a random set of the parameters for this source
        %         region
        toChange=rand(Nregions,1)>.5;
        while sum(toChange)==0%change something
            toChange=rand(Nregions,1)>.5;
        end
        
        nMu=mus;
	%change according to step size (uniform)
        nMu(i,(toChange(:,1)))=mus(i,(toChange(:,1)))+((rand(1,sum(toChange(:,1)))-.5)).*squeeze(C(i,1,(toChange(:,1))))';

        %            lower and upper bounds
        
        %mirror in boundaries. 
        while(sum(nMu(i,:)>highMu|nMu(i,:)<lowMu)>0)
            nMu(i,nMu(i,:)>highMu)=2*highMu-nMu(i,nMu(i,:)>highMu);
            nMu(i,nMu(i,:)<lowMu)=2*lowMu-(nMu(i,nMu(i,:)<lowMu));
        end

        
        
        %have to update the lik for all experiments this region is injected in
        nll=ll;
        k=find(In(:,i)>0);
        nll(k)=loglikLognormMean(In(k,:),Out(k,:),nMu,coef,pdfNoise(k,:),...
            cdfNoise(k,:),zeroMeasured(k,:),volsNonInj(k,:),injOutRegions(k,:));
        if(exp(sum(nll)+priorLoglik(nMu,coef)-sum(ll)-priorLoglik(mus,coef))>1-exp(-m/beta)*rand(1))
            ll=nll;
            mus=nMu;
            acc=acc+1;
            if m<M/10
		%during burnin, update hyper parameters
                C(i,1,toChange(:,1))=C(i,1,toChange(:,1))*1.1;
                C(:,1,:)=min(C(:,1,:),10);
            end
            if(sum(ll)>sum(llmax))
		%record maximum values found - not used in article
                llmax=ll;
                coefmax=coef;
                mumax=mus;
                
            end
        else
            rej=rej+1;
            if m<M/10
                C(i,1,toChange(:,1))=C(i,1,toChange(:,1))/1.1;
                C(:,1,:)=max(C(:,1,:),.001);
            end
        end
        
        
    end
    
    
    accCoef=0;
    rejCoef=0;
    for repe=1:5
%change sigma
        beta_i=1;
        nbeta=coef;
        
        nbeta(beta_i)= (coef(beta_i)+(rand(length(beta_i),1)-.5).*Ccoeff(beta_i));
        nsigm=nbeta;
        %             boundaries
        n=0;
        while(any(nsigm(:)<=lowSigm|nsigm(:)>highSigm))
            nbeta(beta_i)= (coef(beta_i)+(rand(1)-.5)*.1^n);
            nsigm=getSigm(nbeta,mus);
            n=n+1;
        end
        
        
        nll=loglikLognormMean(In,Out,...
            mus,nbeta,pdfNoise,cdfNoise,zeroMeasured,...
            volsNonInj,injOutRegions);
        
        if(exp(sum(nll)+priorLoglik(mus,nbeta)-sum(ll)-priorLoglik(mus,nbeta))>1-exp(-m/beta)*rand(1))
            ll=nll;
            coef=nbeta;
            accCoef=accCoef+1;
            if m<M/10
                Ccoeff(beta_i)=min(Ccoeff(beta_i)*1.1,1);
            end
            if(sum(ll)>sum(llmax))
                llmax=ll;
                coefmax=coef;
                mumax=mus;
                
            end
        else
            rejCoef=rejCoef+1;
            if m<M/10
                Ccoeff(beta_i)=max(Ccoeff(beta_i)/1.1,.001);
            end
        end
        
    end
end




%%provided filename and comment in below to save
%save(file,'atlas','regions','Coefs', 'In', 'Metrics','muNeededNow', 'C','Mus','Ninj','Ninjected','Nregions','m','Out','Qs','in','injOutRegions','volsNonInj','lls','m','pdfNoise','cdfNoise','zeroMeasured','sampFreq','SMALLESTIN','beta','-v7.3')
