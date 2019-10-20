function [caseType,d,k,lgs,lgbranch,nlgs,nlgbranch,totalprob,totalass,content,a3]=chain(inputfile1,inputfile2)
%function []=chain(inputfile1,inputfile2)
% clc 
% clear all
% 
% inputfile1='single1.txt';
% inputfile2='single2.txt';
fin1 = fopen(inputfile1,'rt');
fin2 = fopen(inputfile2,'rt');


%extract intputfile1

if ~feof(fin1)
    str=fgetl(fin1);
    while isempty(str)
        str=fgetl(fin1);
    end
end
caseType=str;

lgs={}; %loop guards 
if ~feof(fin1)
    str=fgetl(fin1);
    while isempty(str)
        str=fgetl(fin1);
    end
end
lgbranch=0;% loopg guard branch number, alwys be 1
while ~isempty(str)
    templgs=str(2:end-1);
    vector=strsplit(templgs,'and');
    lgbranch=lgbranch+1;

    for i=1:length(vector)
        proc=strsplit(vector{1,i},'>');
        lgs{lgbranch,end+1}=proc{1,1};
    end
   
    str=fgetl(fin1);
end
lgs=str2sym(lgs);



nlgs={}; % the branches with guards which violate the loop
nlgbranch=0; % the number of branches which violate the loop
while isempty(str)
    str=fgetl(fin1);
end


while (~isempty(str)&& isnan(str2double(str)))
    templgs=str(2:end-1);
    vector=strsplit(templgs,'and');
    nlgbranch=nlgbranch+1;
    for i=1:length(vector)
        proc=strsplit(vector{1,i},'>');
        vector{1,i}=proc{1,1};
    end
    tempn=str2sym(vector);
    nlgs{nlgbranch}=tempn;
    str=fgetl(fin1);
   
end

if ~feof(fin1)
    sizePV = str2double(str);
end

% Program variable name  
if ~feof(fin1) 
    pv=fgetl(fin1);
    vecpv=strsplit(pv,' ');
    namePV=sym([]);%sym PV name
    for i=1:sizePV
        namePV(end+1)=sym(vecpv{i});
    end
end

% Degree

if ~feof(fin1)
    d = str2double(fgetl(fin1));
end
% k
if ~feof(fin1)
    k = str2double(fgetl(fin1));
end

% %----process gamma-------%
% Gamma_matrix1={};%store Gamma of the loop guards
% gamma1num=length(lgs); % the number of gammas in loop guards
% 
% Gamma_matrix2={};%store Gamma of  guards violating the loop
% gamma2num=[];% the number of gammas in each violated branch 
% 
% for i=1:length(lgs)
%     [At,bt]=equationsToMatrix(lgs(1,i),namePV);
%     At=double(At);
%     bt=double(bt);
%     Gamma_matrix1{i}=[At -bt];
%     
% end
% 
% for i=1:nlgbranch
%     gamma2num(end+1)=length(nlgs{1,i});
%     tempmatrix2={};
%     for j=1:gamma2num(i)
%         this=nlgs{1,i}(1,j);
%         [At,bt]=equationsToMatrix(this,namePV);
%         At=double(At);
%         bt=double(bt);
%         
%         tempmatrix2{j}=[At -bt];
%         
%     end
%     Gamma_matrix2{i}=tempmatrix2;
% end
% %-----got the gamma----%

%extract inputfile2
global a1 a2 a3;

if ~feof(fin2)
    a1 = str2double(fgetl(fin2)); %max counter
end

if ~feof(fin2)
    a2 = str2double(fgetl(fin2)); % the cfg chain size
end

global A AA;
A={}; % store the whole cfg info
AA=[];%store the cfg chain counter

for i=1:a2
    str = fgetl(fin2);
    vector=strsplit(str,'*');
    for j=1:3
        A{i,j}=vector{1,j};
    end
    for j=1:2
        AA(i,j)=str2double(vector{1,j}); 
    end
end

assnumber=0; % the number of  end 'true' in A
assrownum=[]; %store the rownumbers in A of end 'true'
asscounter=[]; %store all assignments counters in the program

probnumber=0; % the number of  end 'prob' in A
probrownum=[]; %store the rownumber in A of end 'prob'
for i=1:a2
    if strcmp(A{i,3},'true')
        assnumber=assnumber+1;
      
        assrownum(end+1)=i;
        asscounter(end+1)=AA(i,1);
    end
    if strfind(A{i,3},'prob')
        probnumber=probnumber+1;
        probrownum(end+1)=i;
    end
        
end

content={};%store details of all assignments from inputfile1
for i=1:assnumber
    str = fgetl(fin1);
    vector=strsplit(str,'@');
    content{i}=vector;
end



global B;
B=[]; % store the rownumbers in A of end  '1'
for i=1:a2
    if strcmp(A{i,2},'1')
        B(end+1)=i;
    end
end

a3=length(B);
C={}; %store the separate chain (with program counters)
CC={};% store the separate chain rownumbers in A

global D DD;
D=[];
DD=[];

for i=1:a3
    u=B(i);
    D=[];
    D(end+1)=AA(u,2);
    
    DD=[];
    DD(end+1)=u;
    
    recur(AA(u,1));
    D=fliplr(D);
    C{i}=D;
    DD=fliplr(DD);
    CC{i}=DD;
end

%find the assignments and prob for each chain in CC
eachass={};
eachprob={};
for i=1:a3
    eachass{i}=intersect(CC{i},assrownum);
    eachprob{i}=intersect(CC{i},probrownum);
end

%extract probability from A
% probinfo={};
% for i=1:probnumber
%     probinfo{i,1}=probrownum(i);
%     %temprob=regexp(A{probrownum(i),3},'\d*\.\d*','match');
%     lloc1=findstr(A{probrownum(i),3},'(');
%     lloc2=findstr(A{probrownum(i),3},')');
%     temprob=A{probrownum(i),3}(lloc1+1:lloc2-1);
% %     if strfind(temprob,'/')
% %         temploc=strfind(temprob,'/');
% %         frac1=temprob(1:temploc-1);
% %         frac2=temprob(temploc+1:length(temprob));
% %         frac1=str2double(frac1);
% %         frac2=str2double(frac2);
% %         temprob=frac1/frac2;
% %     else
% %         temprob=str2double(temprob);
% %     end
%     probinfo{i,2}=temprob;
% end
probinfo=[];
for i=1:probnumber
    probinfo(i,1)=probrownum(i);
    %temprob=regexp(A{probrownum(i),3},'\d*\.\d*','match');
    lloc1=findstr(A{probrownum(i),3},'(');
    lloc2=findstr(A{probrownum(i),3},')');
    temprob=A{probrownum(i),3}(lloc1+1:lloc2-1);
    if strfind(temprob,'/')
        temploc=strfind(temprob,'/');
        frac1=temprob(1:temploc-1);
        frac2=temprob(temploc+1:length(temprob));
        frac1=str2double(frac1);
        frac2=str2double(frac2);
        temprob=frac1/frac2;
    else
        temprob=str2double(temprob);
    end
    probinfo(i,2)=temprob;
end

%store the total probability of each chain
% totalprob={};
% fzall={};
% fmall={};
% for i=1:a3
%     probsize=length(eachprob{i});
%  
%     fz=[];
%     fm=[];
%     decimal=[];
%     for j=1:probsize
%         for jj=1:probnumber
%             if eachprob{i}(j)==probinfo{jj,1}
%                 if strfind(probinfo{jj,2},'/')
%                     thisfrac=probinfo{jj,2};
%                     temploc=strfind(thisfrac,'/');
%                     frac1=thisfrac(1:temploc-1);
%                     frac2=thisfrac(temploc+1:length(thisfrac));
%                     frac1=str2double(frac1);
%                     frac2=str2double(frac2);
%                     fz(end+1)=frac1;
%                     fm(end+1)=frac2;
%                     
%                 else
%                     thisfloat=str2double(probinfo{jj,2});
%                     decimal(end+1)=thisfloat;
%                 end
%                 
%             end
%         end
%       
%     end
%     
%    if ~isempty(fz)
%      totalfz=1;
%     for ii=1:length(fz)
%         totalfz=totalfz*fz(ii);
%         
%     end
%     totalfm=1;
%     for ii=1:length(fm)
%         totalfm=totalfm*fm(ii);
%     end
%     fzall{i}=fz;
%     fmall{i}=fm;
%     factor=Factor(totalfz,totalfm);
%     simplefrac1=totalfz/factor;
%     simplefrac2=totalfm/factor;
%     thisprob1=simplefrac1/simplefrac2;
%    end
%    if ~isempty(decimal)
%        thisprob2=1;
%        for ii=1:length(decimal)
%            thisprob2=thisprob2*decimal(ii);
%        end
%    end
%    
%     if ~isempty(fz) && ~isempty(decimal)
%         thisprob=thisprob1*thisprob2;
%     else
%        if ~isempty(fz)
%            thisprob=thisprob1;
%        else
%            if ~isempty(decimal)
%              thisprob=thisprob2;
%            else
%                thisprob=1;
%            end
%        end
%     end
%         
%        
% %     simplefrac1=num2str(totalfz/factor);
% %     simplefrac2=num2str(totalfm/factor);
%     %totalprob{end+1}=[simplefrac1,'/',simplefrac2];
%     
%     totalprob{end+1}=thisprob;
% end
%sumprob=sum(totalprob);
%o=regexp(A{6,3},'\d*\.\d*','match');
%h=str2double(o{1});
totalprob=[];
for i=1:a3
    probsize=length(eachprob{i});
    thisprob=1;
    for j=1:probsize
        for jj=1:probnumber
            if eachprob{i}(j)==probinfo(jj,1)
                thisprob=probinfo(jj,2)*thisprob;
            end
        end
    end
    totalprob(end+1)=thisprob;
end
sumprob=sum(totalprob);


%store content rownums for each chain
totalass={};
for i=1:a3
    asssize=length(eachass{i});
    totalass{i}=[];
    for j=1:asssize
        for jj=1:assnumber
            if A{eachass{i}(j),1}==content{1,jj}{1,1}
                totalass{i}(end+1)=jj;
            end
        end
    end
     
end
%return caseType,sizePV,namePV,d,k,Gamma_matrix1,gamma1num,lgbranch,Gamma_matrix2,gamma2num,nlgbranch,totalprob,totalass,content;
end

function r=Factor(m,n)
r=mod(m,n);
while r~=0
    m=n;
    n=r;
    r=mod(m,n);
end
r=n;
end

function []=recur(x)
global A AA;
global B;
global D DD;
global a2;

D(end+1)=x;
if x~=1
    
    for j=1:a2
        if x==AA(j,2)
            DD(end+1)=j;
            recur(AA(j,1));
        end
    end
end
end

