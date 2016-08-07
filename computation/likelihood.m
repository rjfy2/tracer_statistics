function [lle, noiseComp, connComp, connContr] = loglikLognormMean(inExp,outExp,mus,sigm,pdfNoise,cdfNoise,zeroMeasured,volsNonInjExp,outInj)
%assume noise and connection dont both contribute. Then it's one of the two; the other must be smaller
%parameter mu is strength of connection, sigma variability
%note that for multiply injected regions with connection, the new strength is not a convolution
%(which would assume independence)

%number of injections and number of injected regions
[ninj, nin]=size(inExp);
%number of regions read from
nout=size(mus,2);


%weighted mu for each region, for those for which a connection exists (repeated for ease of
%computation later). Note that this is usually nearly the same as the strongest mu.

totMu=log10(nansum(permute(10.^repmat(mus,[1,1,ninj]),[3,1,2]).*...
    repmat(bsxfun(@times,inExp,1./sum(inExp,2)),[1,1,nout]),2));
totSigm=log10(nansum(permute(10.^repmat(sigm,[1,1,ninj]),[3,1,2]).*...
    repmat(bsxfun(@times,inExp,1./sum(inExp,2)),[1,1,nout]),2));

% at every iteration, for every experiment, get probability it is noise. Also probability that certain
%     connection is masked by others? For each connection, we want to have a sense of whether this is a
%     real connection or whether it has been masked
%get the contribution to the overall mu for each connection
connContr=bsxfun(@times,permute(10.^repmat(mus,[1,1,ninj]),[3,1,2]).*...
    repmat(bsxfun(@times,inExp,1./sum(inExp,2)),[1,1,nout]),1./(10.^totMu));



%NCD
outMeas=log10(bsxfun(@times,outExp./volsNonInjExp,1./sum(inExp,2)));
% 0's are below the measuring threshold; threshold given by zeroMeasured
nonz=outMeas>-Inf;

if size(inExp,1)==1 %MATLAB issues with column/row
    totMu=squeeze(totMu)';
    totSigm=squeeze(totSigm)';
end

lik=nan(size(outExp));
noiseComp=lik;
connComp=lik;
%noise larger
noiseComp(nonz)=pdfNoise(nonz).*normcdf(squeeze(outMeas(nonz)),squeeze(totMu(nonz)),squeeze(totSigm(nonz)));

%connection larger
connComp(nonz)=cdfNoise(nonz).*normpdf(outMeas(nonz),totMu(nonz),totSigm(nonz));


lik=noiseComp+connComp;

%if 0 is measured, use cdfs:
%NCD
outMeas=log10(bsxfun(@times,zeroMeasured./volsNonInjExp,1./sum(inExp,2)));
%contribution of cdfNoise here is a constant, so leave it out.
lik(~nonz)=normcdf(outMeas(~nonz),totMu(~nonz),totSigm(~nonz));


lik(outInj)=1;%ignore the out regions that were injected
lle=sum(log(lik),2);
end
