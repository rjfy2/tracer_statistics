% adding more and more weak edges, how do graph theoretical measures
% change?

%measures:
%1 betweenness
%2 global efficiency
%3 size max component
noMeas=3;

%methodological:
%3 binary (probably wont make sense)
%2 linear weight
%1 log weight
noMeth=3;


%% get the matrix and thresholds
%this bit uses the Brain Connectivity Toolbox (https://sites.google.com/site/bctnet/)
load('cortInOutCCO50')
load('1901_runs');%based on 5 runs
knowable=false(Ninjected,Nregions);
load('knowable')
medMu(~knowable)=nan;
lowMu(~knowable)=nan;
highMu(~knowable)=nan;

%% convert to percentages

%for entire posterior
mm=10.^MUS;
%convert to expected total signal
%some tricks to extend vols
v=repmat(vols,[1,1,31]);
v=permute(v,[1,3,2]);
mm=bsxfun(@times,mm,v);

musF=log10(bsxfun(@times,mm,1./nansum(mm,3)));% is this the same as estimating percentages straight away?
musF(:,~knowable)=nan;
l=size(musF,1);
ran=floor(l/5):floor(l);

medMu=squeeze(nanmedian(musF(ran,:,:)));

lowMu=squeeze(quantile(musF(ran,:,:),[0.025]));

highMu=squeeze(quantile(musF(ran,:,:),[.975]));

%%
bins=[0:-.2:-4,-5:-1:-8];
noBin=length(bins);
cc=nan(noMeth,noBin);
globEff=nan(noMeth,noBin);
bc=nan(noMeth,noBin,size(medMu,1)*2);
for i_meth=1:noMeth
for i=1:noBin
    binmat = single(medMu>bins(i)); 
    if i_meth==1;
        mat=binmat.*(medMu-min(medMu(:))+.1);
       
    elseif i_meth==2
        mat=binmat.*10.^medMu;
        mat=mat/mean(mat(binmat>0));%make more comparable to the binary case
    elseif i_meth==3
        mat=binmat;
    end
%     mat=mat(:,in);
    mat=mat(:,[in-43,in]);
    mat=[mat(:,end/2+1:end),mat(:,1:end/2);mat];
    mat(isnan(mat))=0;
    mat=mat/max(mat(:));
    
    binmat=binmat(:,[in-43,in]);
    binmat=[binmat(:,end/2+1:end),binmat(:,1:end/2);binmat];
   
     bc(i_meth,i,:) = betweenness_wei(1./mat);
    globEff(i_meth,i) = efficiency_wei(mat);
    if i_meth==3
     globEff(i_meth,i) = efficiency_bin(mat);
    end
       
    [~,siz] = get_components(logical(binmat+binmat'));
    cc(i_meth,i) = max(siz);
end
   
end
   

subplot 131
plot(bins,globEff')
legend({'log','lin','bin'})
subplot 132
plot(bins,cc')
subplot 133
plot(bins,mean(bc,3)')
%% figure: glob eff for both
subplot 111
globEff(globEff==1)=0;
plot(bins, bsxfun(@times,globEff([2,3],:,:),1./max(globEff([2,3],:,:),[],2)),'linewidth',2)
hold all
plot(bins,cc'/max(cc(:)),'--k','linewidth',2)
hold off
xlim([-7.9 0])
xlabel('log percentage weight')
ylabel('global efficiency')
choice=4;

%% only plot one
for choice = 1:2
subplot(2,2,choice)
%choice=1;%2 for linear weights
plot(bins,globEff(choice,:)'/max(globEff(choice,:)),'linewidth',2)
hold all

plot(bins,cc(choice,:)'/max(cc(choice,:)),'linewidth',2)
bb=mean(bc,3);
plot(bins,bb(choice,:)'/max(bb(choice,:)),'linewidth',2)
hold off
% legend('global efficiency','largest connected component','mean betweenness')
xlabel('log percentage weight')
ylabel('normalized value')
xlim([-7.9 -1])
end

%%
wi=9;
hi=9;
set(gcf,'PaperUnits','Centimeters','PaperPosition',[0 0 wi hi],'PaperSize',[wi hi])
print([num2str(choice) '_fig7.pdf'],'-dpdf')

%% split connections up per weight category. different properties?
load('need')
% find centroid
centr=nan(86,3);
for i=1:nreg
    x=find(atlas==atlasIds(i));
    [a,b,c]=ind2sub(size(atlas),x);

    centr(i,:)=mean([a,b,c]);
end
centr=centr([in-43,in],:);
d=squareform(pdist(centr))*.05;

need=squeeze(median(need(:,:,[in-43, in])));
    need=[need(:,end/2+1:end),need(:,1:end/2);need];

    %%


    mat=medMu(:,[in-43,in]);
    mat=[mat(:,end/2+1:end),mat(:,1:end/2);mat];
    regi=regions([in-43,in]);
    siz=4;
hist(medMu(:),50)
bnds=quantile(mat(mat>need),[.05 .95]);
weak=mat<bnds(1)&mat>need;
strong=mat>bnds(2)&mat>need;

%generate random nulls
N=100;
dsrand=[];%distance values found;
distrand=[];%degree values found;
allowed=find(mat>need);
nn=length(allowed);
for n=1:N
    randmat=false(size(mat));
    randmat(allowed(randsample(nn,sum(strong(:)))))=1;%set the same number of edges to 1
    dsrand=[dsrand; d(randmat)];
    distrand=[distrand,sum(randmat+randmat')];%could use bins for performance...
end
%%
subplot('Position',[0 .5 .5 .5])
gplot(weak,centr(:,[2,3]))
hold on
scatter(centr(:,2),centr(:,3),siz+sum(weak+weak')*siz,'f')
hold off
axis equal
axis off

set(gca,'XTick',[])
set(gca,'YTick',[])
% xlabel('Axial view of weak connections')

subplot('Position',[0.5 .5 .5 .5])
gplot(strong,centr(:,[2,3]))
hold on
scatter(centr(:,2),centr(:,3),siz+sum(strong+strong')*siz,'f')
hold off
strongest=find(sum(strong+strong')>6);
%remove the doubles, leave some in each hemi
% strongest([1:2:end/2,end/2+2:2:end])=[];
strongest([1+end/2:end])=[];
text(centr(strongest,[2])+3,centr(strongest,[3])-4,regi(strongest),'FontSize',9,'fontweight','bold')
axis off
set(gca,'XTick',[])
set(gca,'YTick',[])
axis equal
% xlabel('Axial view of strong connections')
subplot 223
histogram(d(weak),'facealpha',.5)
hold on
histogram(d(strong),'facealpha',.5)
h=ylim;
xs=xlim;
dd=.5;
bi=xs(1):dd:xs(2);
f=zeros(length(bi),1);
for i=1:length(bi)
    f(i)=sum(bi(i)-dd/2<dsrand&dsrand<bi(i)+dd/2);
end
plot(bi,f/N/dd,'k','linewidth',2);
%differen distributions?
[~,p] = kstest2(dsrand,d(strong))
[~,p] = kstest2(dsrand,d(weak))
[~,p] = kstest2(d(weak),d(strong))
hold off
ylabel('Frequency')
xlabel('Euclidean distance (mm)')
subplot 224
histogram(sum(weak+weak'),0:1/3:10,'facealpha',.5)
hold on
histogram(sum(strong+strong'),-1/3+.01:1/3:10,'facealpha',.5)
xlim([-.4 10])
h=ylim;
xs=xlim;
dd=1;
bi=xs(1):dd:xs(2);
f=zeros(length(bi),1);
for i=1:length(bi)
    f(i)=sum(bi(i)-dd/2<distrand&distrand<bi(i)+dd/2);
end
plot(bi,f/N,'k','linewidth',2);
hold off
[~,p] = kstest2(distrand,sum(strong+strong'))
[~,p] = kstest2(distrand,sum(weak+weak'))
[~,p] = kstest2(sum(weak+weak'),sum(strong+strong'))
ylabel('Frequency')
xlabel('Degree')
%%
wi=12;
hi=12;
set(gcf,'PaperUnits','Centimeters','PaperPosition',[0 0 wi hi],'PaperSize',[wi hi])
print('fig6.pdf','-dpdf')
