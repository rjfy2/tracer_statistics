%read the cocomac values
fid=fopen('matrix.txt');
i=0;
j=1;
nregs=58;
vals=cell(nregs);
line=fgetl(fid);
while(isempty(line)||any(line~=-1))
    line=fgetl(fid);
    
    if(~isempty(line))
        i=i+1;
        j=0;
    end
    rest=line;
    while(~isempty(rest)&&any(rest~=-1))
        j=j+1;
        [bit,rest]=strtok(rest,',');
        vals{i,j}=bit;
    end
end
fclose(fid);

%% now go through the text fields to get the values out
coco=nan(nregs,nregs,350);
cocoNames=cell(nregs,1);
for i=1:nregs;
    %get the name out
    cocoNames{i}=vals{i,1}(7:end-1);
    for j=1:nregs
        if i~=j
        if(~strcmp(vals{i,j+1},'-Inf'))%not known
            if(strcmp(vals{i,j+1},'1.23'))%X ('exists')
                %single value not interesting for consistency, ignore
                coco(i,j,1)=101.23;
            elseif length(vals{i,j+1})==1
                %single value not interesting for consistency, ignore
            else
                for k=2:length(vals{i,j+1})-1%skip the quotation marks
                    val=(vals{i,j+1}(k));
                    if strcmp(val,'X')
                        val=101.23;
                    else
                        val=str2num(val)+100;
                    end
                    if ~isempty(val)
                     coco(i,j,k-1)=val;
                    end

                end
            end
        end
        end
    end
end
%plot variance
imagesc(squeeze(nanvar(permute(coco,[3,1,2]))))

%no of experiments
noExp=squeeze(sum(~isnan(coco),3));
plot(noExp,squeeze(nanvar(permute(coco,[3,1,2]))),'.')

save('cocovals','coco','vals','noExp','nregs','cocoNames')