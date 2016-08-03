
clear all
load('estimates.csv')
load('CI_high.csv')
load('CI_low.csv')
load('noise_high.csv')
load('noise_low.csv')
load('coinj_high.csv')
load('coinj_low.csv')
load('coinj_median.csv')
load('noise_median.csv')
as=[11 1];
nreg=86;
for j=1:2
    for hemi=1:2
        subplot(2,2,j+(hemi-1)*2)
        a=as(j);
        
        intra=nreg/2+1:nreg;
        [~,ord]=sort(estimates(a,intra));
        if hemi==1
            ord=ord+nreg/2;
        end
        th=20;
        scatter(1:nreg/2,estimates(a,ord),th,[0 0 0],'f')
        hold on
        for i=1:nreg/2
            %plot the noise threshold
            width=.7;
            lowN=max(squeeze(noise_low(a,ord(i))),-20);
            highN=squeeze(noise_high(a,ord(i)));
            rectangle('Position',[i-width/2 lowN, width highN-lowN],'Curvature',.5,'FaceColor','red','edgecolor','none')
            
            %plot the coinjection threshold
            lowO=max(squeeze(coinj_low(a,ord(i))),-20);
            highO=squeeze(coinj_high(a,ord(i)));
            if ~(highO==-Inf)%if threshold is not absent
                rectangle('Position',[i-width/2 lowO, width highO-lowO],'Curvature',.5,'FaceColor','blue','edgecolor','none')
            end
            
            %if overlap, make purple
            if  ~((lowN>highO)|highN<lowO|isnan(highO))
                rectangle('Position',[i-width/2 max(lowO,lowN), width min(highN,highO)-max(lowO,lowN)],'Curvature',.5,'FaceColor',[.5 0 .5],'edgecolor','none')
            end
            
            %plot CI
            plot(([i i]),[CI_low(a,ord(i)), CI_high(a,ord(i))],'k','linewidth',1)
            yl=[-8 1];
            ylim(yl)
            xlim([1 nreg/2]-.5)
            xlabel('target regions (right hemisphere)')
            ylabel('log NCD (mm^{-3})')
        end
        
        hold off
    end
end