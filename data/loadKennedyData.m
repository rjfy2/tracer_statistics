function [In, Out, FLN, source, unS, target, names, data, PNASdata,total] = loadKennedyData()
%load the xls first - obtainable from core-nets.org
data=importdata('data Kennedy.xls');
%remove header
data.data(1,:)=[];
data.textdata(1,:)=[];


%% same for their provided xls
PNASdata=importdata('PNAS_2013.xls');
%remove header
PNASdata.textdata(1,:)=[];

%NOTE MATLAB has difficulties with some area names (INS, PI, that do work
%for the inital file, which has different abbreviations (insula,
%parainsula). So use that one, rest seems to be identical


%% get all the sources and targets
source=data.data(:,3);%ones in numeric format
un=unique(data.textdata(:,3));

for i=2:length(un)%skip first, it's: ''
    same=strcmp(data.textdata(:,3),un(i));
    source(same)=1000+i;%assign id
end
%%
target=data.data(:,4);

%compare to the same uniques to get the same IDs
for i=2:length(un)%skip first, it's: ''
    same=strcmp(data.textdata(:,4),un(i));
    target(same)=1000+i;%assign id
end

%% get names out
namesList=data.textdata(:,3);
for i=1:length(namesList)
    if ~isnan(data.data(i,3))
        namesList{i}=num2str(data.data(i,3));
    end
end
unS=unique(source);
names=cell(91,1);
for i=1:91
    names{i}=namesList{find(source==unS(i),1)};
end
%% make In/Out
    
FLN=data.data(:,5);
Neur=data.data(:,6);
In=zeros(39,91);
Out=zeros(39,91);

%total neurons per experiment
total=zeros(39,1);
unS=unique(source);
for i=1:39
    ids=find(data.data(:,1)==i);
    In(i,unS==target(min(ids)))=1;
    for j=1:length(ids)
        id=source(ids(j));
        Out(i,unS==id)=FLN(ids(j));
        total(i)=total(i)+Neur(ids(j));
    end
end

end