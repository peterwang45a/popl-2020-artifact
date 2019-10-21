function [outputdata] = synthesis(inputfile1,inputfile2,inputfile3)
% clc
% clear all
t1=clock;
%inputfile='american.txt';
%outputformula='output2.txt';

%fin = fopen(inputfile,'rt');
%fout = fopen(outputformula,'wt+');
% Input from a file
% inputfile1='double1.txt';
% inputfile2='double2.txt';
% inputfile3='doubleconfig.txt';
%outputfile='vmini.mat';
%initial_values=[0 5];
[caseType,d,k,lgs,lgbranch,nlgs,nlgbranch,totalprob,totalass,content,a3]=chain(inputfile1,inputfile2);

% per=[1/304 1/152 3/304 1/76 5/304 3/152 6/152 1/152 9/152 1/152 249/304];
% totalprob=per;
%--process the initial valuations---%
fin3 = fopen(inputfile3,'rt');
str = fgetl(fin3);
vector=strsplit(str,' ');
sizePV=length(vector);
namePV=sym([]);
initial_values=[];
for i=1:sizePV
    thisconfig=strsplit(vector{1,i},'=');
    namePV(end+1)=str2sym(thisconfig{1,1});
    initial_values(end+1)=str2double(thisconfig{1,2});
end
%---Done----%
syms K KK epsilon CC
options = optimoptions('linprog','Display','none');

global target;
%target=[5 0];
target =initial_values;% initial values of all program variables,if one PV hasn't initial value,just input zero
lenoftarget=length(target);%it is equal to sizePV
% Number of program variables

% if ~feof(fin)
%       caseType = fgetl(fin);% affine or linear
% end


% if ~feof(fin)
%     sizePV = str2double(fgetl(fin));
% end

% % Program variable name
% if ~feof(fin)
%     pv=fgetl(fin);
%     vecpv=strsplit(pv,' ');
%     namePV=sym([]);%sym PV name
%     for i=1:sizePV
%         namePV(end+1)=sym(vecpv{i});
%     end
% end

% Degree

% if ~feof(fin)
%     d = str2double(fgetl(fin));
% end
% % k
% if ~feof(fin)
%     k = str2double(fgetl(fin));
% end

global kk;
kk=k;

% % Number of chains
% if ~feof(fin)
%     sizeL = str2double(fgetl(fin));
% end
sizeL=a3;% Number of chains

f = [];
lb = [];

global degree;
degree=d;%template degree
global saveii;
saveii=[];
global output;
output=[];%���ɴ�С���ֵ�������f�����������h1��������i.e.X^3,X^2*y�Ĵ�С
global t;
t=lenoftarget;


%only get the template h1's coeff
if d==1
    for i=1:lenoftarget
        
        f(end+1)=target(i);
        
    end
    f(end+1)=1;
    
else
    iter(lenoftarget,1);
    f=output;
end

%f
f1=f;%only h1's coeff
lenoff1=length(f1);%ÿ��label���ı�������


target=namePV;
saveii=[];
output=sym([]);
iter(sizePV,1);
PVs=output;


% Template h's coeff
templates_coeff = sym([]);%��������out��ǰ��template coeff,����һ��lout��û�и�ֵ
%for i=1:sizeL
for j=1:lenoff1
    syms([char(97),'_',num2str(j)]);
    templates_coeff(end+1) = [char(97),'_',num2str(j)];%��ÿһ��label�����в�����ֵ, i.e a_i_j, label i����j������
end
%end
%templates_coeff

templates=templates_coeff*PVs.';%the whole template h,include all the labels
% syms([char(97),'_',num2str(sizeL),'_',num2str(1)]);
% %out����һ��������a_sizeL_1,
% templates(sizeL)=[char(97),'_',num2str(sizeL),'_',num2str(1)];



%��template�Ĳ�����lb��ֵ,����out���ĳ�����
% for i=1:(sizeL-1)
%     for j=1:lenoff1
%
%             if j==1
%                lb(j+(i-1)*lenoff1)=0.001;%ÿһ�е�һ������lb>=0.001, i.e a_1_1,a_2_1,a_3_1...Ϊ�˱�֤��һ��PV����ߴ���������
%             else
%                lb(j+(i-1)*lenoff1)= -inf;
%             end
%
%     end
%     %lb(i*elemNum)= -inf;%ÿһ�е�ĩβ������
% end
% lb((sizeL-1)*lenoff1+1)=-inf;%h(lout)����a_sizeL_1��lb

%for i=1:sizeL
for j=1:lenoff1
    lb(end+1)= -inf;
end
%end
%lb(1)=0.0001;%��һ��label������ߴ��Ϊ0
%lb((sizeL-1)*lenoff1+1)=-inf;%h(lout)����a_sizeL_1��lb

%lb


preExpectation={};%save pre of the whole program,branching label has two pre becasue it has two next labels,
%other labels has one pre.

trueg={};% g=h-pre

gamma={};%the whole gamma of the program,  (sizeL_1)*max(size of next labels)
allmonoid={};%each monoid_i of the corresponding label of the whole program
ucoeff={};%each monoid_i's coeff of the whole program
handelman={};%handelman g'

subtract={};%g-g'
simplesub={};%collect subtract{} according to namePV
equations={};%equations according to simplesub,which is the coeffs of g-g'
terms={};%the corresponding terms of this polynomial with respect to namePV
updatefunction={};
sampleinf={};%ÿ����֧��sample����������
distribution={};%sample �����ķֲ���Ϣ


%solve A1
%  rsm={};
%  rgamma={};
%  rallmonoid={};
%  rucoeff={};
%  rhandelman={};
%  rsubtract={};
%  rsimplesub={};
%  requations={};
%  rterms={};

%%solve A3
%validnum=0;%��skip��������Чassignmentִ�з�֧�ĸ���
eachbound={};
samplegamma={};
thesamples={};
eachv={};% updated valations at each label
diseachv={};%b-b'
newf={};
newlb={};
newub={};
newAneq={};
newbneq={};

%�ȶ�������loop guard��gamma
% lengamma = str2double(fgetl(fin));% number of elements in this Gamma
% thisgamma=sym([]);
% for j=1:lengamma
%     coefinv= fgetl(fin);%the first invariant's coeff in this Gamma
%     vector = eval(['[',coefinv,']']);
%     inv=vector*[namePV,1].';
%
%     thisgamma(j)=inv;
% end
lengamma=length(lgs);
thisgamma=lgs;
gamma=lgs;

strnamePV={};
for i=1:sizePV
    strnamePV{i}=char(namePV(i));
end

finalmax=[];%ÿ��label�������bounded update
collect_initialv={};
collect_exb={};
for i=1:sizeL
    
    %tlen=lenoff1;%template�Ѿ�ռ����ô�������
    samplegamma1=sym([]); %sampling variables' gamma in this chain
    newsample=sym([]); %sampling variables name in this chain
    initialv=namePV; %use initialv to update F(l,b,r)
    exb=initialv;%use the exb to update E_r(F(l,b,r))
  
    num=length(totalass{1,i}); %the number of assignments in this chain
    update_pre=templates;%�������´�������template function
    update_org=templates;%�������²���������template function
    for inum=1:num
        if ~strcmp(content{1,totalass{1,i}(inum)}{1,2},'skip')
            
            str = content{1,totalass{1,i}(inum)}{1,2};
            vector=strsplit(str,'=');
            updatef=str2sym(vector);%ȡ����update function,�����Ż� '='��ߵ�update�����ұߵ�update���ݣ���x=x+r,updatef(1)��x,updatef(2)��x+r
            left_side=updatef(1);%��� x
            right_side=updatef(2);%�ұ� x+r
            
            left=char(left_side);
            
        
            for nn=1:sizePV
                initialv(nn)=subs(initialv(nn),left_side,right_side);
                exb(nn)=subs(exb(nn),left_side,right_side);
            end
         
            
            distributionType=content{1,totalass{1,i}(inum)}{1,3};
            switch distributionType
                case 'normal'
                    %update_pre=update_pre;
                    initialv=initialv;
                    exb=exb;
                case 'abnormal'
                    str =content{1,totalass{1,i}(inum)}{1,4}; %ȡ������sample ����r1��r2...�����ҷ��Ż�
                    vector=strsplit(str,' ');
                    r_sam=str2sym(vector);
                    size_r=length(r_sam);%sample ��������
                    sampleinf{i,inum}=r_sam; %��ʵÿ��updateֻ��һ��sample����
                    %thesamples{i,1}=[thesamples{i,1},r_sam];
                    newsample(end+1)=r_sam; %�ռ���i�µ�����sampling variables
                    %�������Ե����ԣ�����ɢ����ֲ�DU���������ȷֲ�CU
                    str =content{1,totalass{1,i}(inum)}{1,5};
                    vector =strsplit(str,' ');
                    properties= vector;
                    for j=1:size_r
                        
                        type=properties{1,j};
                        switch type
%                             case 'DS'
%                                 [cr,tr]=coeffs(update_pre,r_sam(j));%trΪ��ǰsample ������Ӧ��terms��������r,��tr=[r^4,r^2,...r,1]
%                                 % Prob p1,p2....
%                                 str =content{1,totalass{1,i}(inum)}{1,6};
%                                 vector = eval(['[',str,']']);
%                                 prob = vector;
%                                 % ��ǰsample ��������ɢȡֵ,������ĸ��ʶ�Ӧ
%                                 str =content{1,totalass{1,i}(inum)}{1,7};
%                                 vector = eval(['[',str,']']);
%                                 rv = vector;
%                                 
%                                 for jj=1:(length(tr)-1)  %���ȼ�ȥ���һ���1
%                                     rv_new=double(subs(tr(jj),r_sam(j),rv));
%                                     Er3=prob*rv_new.';%���sample����ָ���ľ�ֵ����E(r^2),E(r^3)...
%                                     update_pre=subs(update_pre,tr(jj),Er3);%����pre�и�sample����Ϊ���ֵ
%                                 end
                            case 'DU'
                                %ȡ���ֲ��Ľ�(a,b)
                                str =content{1,totalass{1,i}(inum)}{1,6};
                                vector = eval(['[',str,']']);
                                bound = vector;
                                samplegamma1=[samplegamma1,r_sam(j)-bound(1),bound(2)-r_sam(j)]; %�ռ���i�µ�����sampling variable�� gamma, ��r~(1,2),����r-1,2-r
                                %eachbound{i,1}=[eachbound{i,1},bound];
                                for z=1:sizePV
                                    try
                                        [cr,tr]=coeffs(exb(z),r_sam(j));%trΪ��ǰsample ������Ӧ��terms��������r,��tr=[r^4,r^2,...r,1]
                                        for jj=1:(length(tr)-1)
                                            s=coeffs(tr(jj),r_sam(j),'All');
                                            size=length(s)-1;%��ǰtr(jj)��ָ����С
                                            %bound_new=subs(tr(jj),r_sam(j),bound);
                                            %Er=(bound(2)^size-bound(1)^size)/size*(bound(2)-bound(1));
                                            Er1=sum(subs(tr(jj),r_sam(j),bound(1):bound(2)));
                                            Er2=Er1/(bound(2)-bound(1)+1);
                                            exb(z)=subs(exb(z),tr(jj),Er2);%����pre�и�sample����Ϊ���ֵ
                                        end
                                    end
                                end
                            case 'CU'
                                %ȡ���ֲ��Ľ�(a,b)
                                str=content{1,totalass{1,i}(inum)}{1,6};
                                vector = eval(['[',str,']']);
                                bound = vector;
                                samplegamma1=[samplegamma1,r_sam(j)-bound(1),bound(2)-r_sam(j)];
                                % eachbound{i,1}=[eachbound{i,1},bound];
                                distribution{i,inum}=bound;
                                %Er=(bound(1)+bound(2))/2;
                                for z=1:sizePV
                                    try
                                        [cr,tr]=coeffs(exb(z),r_sam(j));%trΪ��ǰsample ������Ӧ��terms��������r,��tr=[r^4,r^2,...r,1]
                                        for jj=1:(length(tr)-1)
                                            s=coeffs(tr(jj),r_sam(j),'All');
                                            size=length(s);%��ǰtr(jj)��ָ��+1
                                            %bound_new=subs(tr(jj),r_sam(j),bound);
                                            Er=(bound(2)^size-bound(1)^size)/(size*(bound(2)-bound(1)));
                                            
                                            exb(z)=subs(exb(z),tr(jj),Er);%����pre�и�sample����Ϊ���ֵ
                                        end
                                    end
                                end
                        end
                        
                    end
                    exb=exb;
            end
        end
    end
    
    collect_initialv{i,1}=initialv;
    collect_exb{i,1}=exb;
    
    samplegamma{i,1}=samplegamma1; %sampling variables' gamma in this chain
    thesamples{i,1}=newsample; %sampling variables name in this chain
    
    preExpectation{i,1}=subs(update_pre,namePV,exb);
    updatefunction{i,1}=subs(update_org,namePV,initialv);
    
    eachv{i,1}=initialv;% the F(l,b,r)
    diseachv{i,1}=eachv{i,1}-namePV; %F(l,b,r)-b, ��һ����diseachv��ÿһ��ľ���ֵ�����ֵmax
    
    renewnamePV=[namePV,thesamples{i,1}];
    %all gammas in this chain
    invgamma=-gamma;
    [A,b]=equationsToMatrix(invgamma,renewnamePV);
    newAneq{i,1}=A;
    newbneq{i,1}=b;
    
    invsamplegamma=-samplegamma{i,1};
    [A,b]=equationsToMatrix(invsamplegamma,renewnamePV);
    newAneq{i,1}=[newAneq{i,1};A];
    newbneq{i,1}=[newbneq{i,1};b];
    newAneq{i,1}=double(newAneq{i,1});
    newbneq{i,1}=double(newbneq{i,1});
    
    xx=diseachv{i,1}; %F(l,b,r)-b
    thismax=[]; %ÿ��program variable�����bounded update
    for xnum=1:sizePV
        
        [A,b]=equationsToMatrix(xx(xnum),renewnamePV);
        newf{i,xnum}=A;
        constant=-b;
        newf{i,xnum}=double(newf{i,xnum});
        %min
        [x,fval]=linprog(newf{i,xnum},newAneq{i,1},newbneq{i,1},[],[],[],[],options);
        absmin=abs(fval+constant);
        %max
        [x,fval]=linprog(-newf{i,xnum},newAneq{i,1},newbneq{i,1},[],[],[],[],options);
        absmax=abs(-fval+constant);
        maxabs=max(absmin,absmax);
        thismax(end+1)=maxabs;
    end
    finalmax(end+1)=max(thismax); %��label i������bounded update
    %end
    
end
boundedupdate=max(finalmax);

updateall=[];
for i=1:sizeL
     update1=[];
    for j=1:sizePV
        [al bl]=coeffs(eachv{i,1}(1,j),namePV);
        try
            if double(bl(end))==1
                pal=al(1:length(al)-1);
                pbl=bl(1:length(bl)-1);
                pal=abs(pal);
                sub=find(pal<0.05);
                pal(sub)=[];
                m1=max(pal);
                this=m1*length(m1);
            end
            
        catch
            al=abs(al);
            sub=find(al<0.05);
            al(sub)=[];
            m1=max(al);
            
            this=m1*length(al);
        end
        update1(end+1)=this;
    end
    updateall(end+1)=max(update1);
end
thisL=max(updateall);
if thisL<1.5 && thisL>0.5
    thisL=1;
end



% str = fgetl(fin);
% vector = eval(['[',str,']']);
%prob = vector;%����֧��ִ�и���
prob=totalprob;
expectations=0;%E(h(F()))
for i=1:sizeL
    expectations=expectations+prob(i)*preExpectation{i,1};
end

trueg= templates-expectations-epsilon; %A3�� trueg=h-E(h(F()))-epsilon>=0
trueg=collect(trueg,namePV);


degree=kk;%gamma's max multicand
saveii=[];%clear saveii before this monoid iter
t=lengamma;
output=sym([]);
target=thisgamma;
iter(lengamma,1);
allmonoid=output;
lenofmonoid=length(allmonoid);%length of this monoid
eachu_coeff=sym([]);

%tlen=tlen+lenofmonoid;%�������´�ѭ��������һ��next labelʱ����������������һ������ϵ������������һ��next label��monoid�ķ���ϵ���ص�
for m=1:lenofmonoid
    syms([char(98),'_',num2str(m)]);
    eachu_coeff(end+1)= [char(98),'_',num2str(m)];
    %eachu_coeff(end+1)=coeff;%coeff corresponding to each u in the current monoid
end

ucoeff=eachu_coeff;
%ÿ��u_i���ˣ�ÿ��u_i��ϵ��Ҳ���ˣ����濪ʼ����Handelman��g
handelman=ucoeff*allmonoid.';%g' of Handelman

subtract=trueg-handelman;%g-g'
simplesub=collect(subtract,namePV);%�ϲ�ͬ����
[c,te]=coeffs(simplesub,namePV);
equations=c;
terms=te;


%%end A3




%%solve A1

lengamma=length(lgs);
thisgamma=lgs;
rgamma1=lgs;

degree=kk;%gamma's max multicand
saveii=[];%clear saveii before this monoid iter
t=lengamma;
output=sym([]);
target=thisgamma;
iter(lengamma,1);
rallmonoid1=output;
lenofmonoid=length(rallmonoid1);%length of this monoid
eachu_coeff=sym([]);

%tlen=tlen+lenofmonoid;%�������´�ѭ��������һ��next labelʱ����������������һ������ϵ������������һ��next label��monoid�ķ���ϵ���ص�
for m=1:lenofmonoid
    syms([char(99),'_',num2str(m)]);
    eachu_coeff(end+1)= [char(99),'_',num2str(m)];
    %eachu_coeff(end+1)=coeff;%coeff corresponding to each u in the current monoid
end

rucoeff1=eachu_coeff;
%ÿ��u_i���ˣ�ÿ��u_i��ϵ��Ҳ���ˣ����濪ʼ����Handelman��g
rhandelman1=rucoeff1*rallmonoid1.';%g' of Handelman


rsubtract1=templates-rhandelman1;%g-g'

rsimplesub1=collect(rsubtract1,namePV);%�ϲ�ͬ����
[c,te]=coeffs(rsimplesub1,namePV);
requations1=c;
rterms1=te;
%end A1


%%solve A4
rgamma4=rgamma1;
rallmonoid4={};
rucoeff4={};
rhandelman4={};
rsubtract4={};
rsimplesub4={};
requations4={};
rterms4={};
rgamma4new={};
for i=1:sizeL
    nexteta=updatefunction{i,1};
    sampler= thesamples{i,1};
    newnamePV=[namePV,sampler];
    
    rgamma4new{i,1}=[rgamma4,samplegamma{i,1}];
 
    lengamma=length(rgamma4new{i,1});
    
    degree=kk;%gamma's max multicand
    saveii=[];%clear saveii before this monoid iter
    t=lengamma;
    output=sym([]);
    target=rgamma4new{i,1};
    iter(lengamma,1);
    rallmonoid4{i,1}=output;
    lenofmonoid=length(rallmonoid4{i,1});%length of this monoid
    
    eachu_coeff=sym([]);
    for m=1:lenofmonoid
        syms([char(102),'_',num2str(i),'_',num2str(m)]);
        eachu_coeff(end+1)= [char(102),'_',num2str(i),'_',num2str(m)];
        %eachu_coeff(end+1)=coeff;%coeff corresponding to each u in the current monoid
    end
    
    rucoeff4{i,1}=eachu_coeff;
    %ÿ��u_i���ˣ�ÿ��u_i��ϵ��Ҳ���ˣ����濪ʼ����Handelman��g
    rhandelman4{i,1}=rucoeff4{i,1}*rallmonoid4{i,1}.';%g' of Handelman
    rsubtract4{i,1}=CC-nexteta+templates-rhandelman4{i,1};%g=c-eta(b')+eta(b),  g-g'
    
    rsimplesub4{i,1}=collect(rsubtract4{i,1},newnamePV);%�ϲ�ͬ����
    [c,te]=coeffs(rsimplesub4{i,1},newnamePV);
    requations4{i,1}=c;
    rterms4{i,1}=te;
    
    eachu_coeff=sym([]);
    for m=1:lenofmonoid
        syms([char(103),'_',num2str(i),'_',num2str(m)]);
        eachu_coeff(end+1)= [char(103),'_',num2str(i),'_',num2str(m)];
        %eachu_coeff(end+1)=coeff;%coeff corresponding to each u in the current monoid
    end
    
    rucoeff4{i,2}=eachu_coeff;
    %ÿ��u_i���ˣ�ÿ��u_i��ϵ��Ҳ���ˣ����濪ʼ����Handelman��g
    rhandelman4{i,2}=rucoeff4{i,2}*rallmonoid4{i,1}.';%g' of Handelman
    rsubtract4{i,2}=nexteta-templates+CC-rhandelman4{i,2};%g=eta(b')-eta(b)+c,  g-g'
    
    rsimplesub4{i,2}=collect(rsubtract4{i,2},newnamePV);%�ϲ�ͬ����
    [c,te]=coeffs(rsimplesub4{i,2},newnamePV);
    requations4{i,2}=c;
    rterms4{i,2}=te;
end
%%end A4
%end

%solve new A2,K,K'
rgamma2={};
rallmonoid2={};
rucoeff2={};
rhandelman2={};
rsubtract2={};
rsimplesub2={};
requations2={};
rterms2={};
subsupdate={};%A2�½�һֱPV�����������º�õ��� F(b)
subsup={};
%numofA2=str2double(fgetl(fin)); % b����loop guard, F(b)������loop guard�ķ�֧����
numofA2=nlgbranch;
for i=1:numofA2
  
    lengamma=length(nlgs{1,i});
    thisgamma=nlgs{1,i};
    rgamma2{i,1}=nlgs{1,i};
    
    degree=kk;%gamma's max multicand
    saveii=[];%clear saveii before this monoid iter
    t=lengamma;
    output=sym([]);
    target=thisgamma;
    iter(lengamma,1);
    rallmonoid2{i,1}=output;
    lenofmonoid=length(rallmonoid2{i,1});%length of this monoid
    % update_temp=collect(update_temp,namePV);
    
    eachu_coeff=sym([]);
    %templates>=K
    for m=1:lenofmonoid
        syms([char(100),'_',num2str(i),'_',num2str(m)]);
        eachu_coeff(end+1)= [char(100),'_',num2str(i),'_',num2str(m)];
        %eachu_coeff(end+1)=coeff;%coeff corresponding to each u in the current monoid
    end
    rucoeff2{i,1}=eachu_coeff;
    %ÿ��u_i���ˣ�ÿ��u_i��ϵ��Ҳ���ˣ����濪ʼ����Handelman��g
    rhandelman2{i,1}=eachu_coeff*rallmonoid2{i,1}.';%g' of Handelman
    rsubtract2{i,1}=templates-K-rhandelman2{i,1};%g-g'
    
    rsimplesub2{i,1}=collect(rsubtract2{i,1},namePV);%�ϲ�ͬ����
    [c,te]=coeffs(rsimplesub2{i,1},namePV);
    requations2{i,1}=c;
    rterms2{i,1}=te;
    
    %templates<=K'
    eachu_coeff=sym([]);
    for m=1:lenofmonoid
        syms([char(101),'_',num2str(i),'_',num2str(m)]);
        eachu_coeff(end+1)= [char(101),'_',num2str(i),'_',num2str(m)];
        %eachu_coeff(end+1)=coeff;%coeff corresponding to each u in the current monoid
    end
    
    rucoeff2{i,2}=eachu_coeff;
    %ÿ��u_i���ˣ�ÿ��u_i��ϵ��Ҳ���ˣ����濪ʼ����Handelman��g
    rhandelman2{i,2}=rucoeff2{i,2}*rallmonoid2{i,1}.';%g' of Handelman
    
    rsubtract2{i,2}=KK-templates-rhandelman2{i,2};%g-g'
    
    rsimplesub2{i,2}=collect(rsubtract2{i,2},namePV);%�ϲ�ͬ����
    [c,te]=coeffs(rsimplesub2{i,2},namePV);
    requations2{i,2}=c;
    rterms2{i,2}=te;
    
end

%end

%end
%end A2





% totalnum=sizeL*lenoff1+numh;%total numbers of the whole program coeffs
totalsym=sym([]);%the total coeffs of the whole program,the order is: template h1,h2,..hn ��PV���ֵ�������a1,b1�ȵ�;
%handelman1 ,handelman2,...��next label˳��ÿ��gamma��invariant���ֵ�˳������coeff
totalsym1=templates_coeff;
totalsym2=ucoeff;  % A3 coeffs
totalsym3=rucoeff1; % A1 coeffs
totalsym4=sym([]);  %A2 coeffcients

totalsym5=sym([]); %A4 coeffcients



for i=1:numofA2
    eu=[rucoeff2{i,:}];
    totalsym4=[totalsym4,eu];
end


for i=1:sizeL
    euu=[rucoeff4{i,:}];
    totalsym5=[totalsym5,euu];
end
totalsym=[totalsym1,totalsym2,totalsym3,totalsym4,totalsym5,K,KK,epsilon,CC];

Aneq=[];
bneq=[];
Aeq=[];
beq=[];

%add A3����
[A,b]=equationsToMatrix(equations,totalsym);
Aeq=[Aeq;A];
beq=[beq;b];



%add A1����
[A,b]=equationsToMatrix(requations1,totalsym);
Aeq=[Aeq;A];
beq=[beq;b];



%add A2����
for i=1:numofA2
    equas=[requations2{i,:}];
    [A,b]=equationsToMatrix(equas,totalsym);
    Aeq=[Aeq;A];
    beq=[beq;b];
end




%add A4����
for i=1:sizeL
    equass=[requations4{i,:}];
    [A,b]=equationsToMatrix(equass,totalsym);
    Aeq=[Aeq;A];
    beq=[beq;b];
end


f=[f,zeros(1,length(totalsym)-length(totalsym1))];


lb=[lb,zeros(1,length(totalsym)-length(totalsym1)-4),-inf,0,1,0];


ub=[];
for i=1:length(totalsym1)
    ub(end+1)=inf;
end
for i=1:(length(totalsym)-length(totalsym1)-4)
    ub(end+1)=inf;
end
ub=[ub,0,0,1,inf];


f=f.';
lb=lb.';
ub=ub.';


% Aneq=double(Aneq);
% bneq=double(bneq);
Aeq=double(Aeq);
beq=double(beq);
[x,fval] =linprog(f,[],[],Aeq,beq,lb,ub,options);
%x=linsolve(A,b);
%end
t2=clock;
time=etime(t2,t1);
%fval=fval;

%---output the values in Table 2,3 -----%
outputdata={};
etacoeff=x(1:lenoff1);
for i=1:lenoff1-1
    if abs(etacoeff(i))<1e-4
        etacoeff(i)=0;
    end
end
eta=vpa(PVs*etacoeff);
% etacoeff1=roundn(etacoeff,-4);
% eta=vpa(PVs*etacoeff1);

Boundedupdate=boundedupdate;
thecoeffs=x(length(x)-3:length(x));
theK=thecoeffs(1);
epi=thecoeffs(3);

%thisM=sizePV*max(abs(x(1:lenoff1-1)));
[am bm]=coeffs(eta,namePV);
try
    if double(bm(end))==1
        pam=am(1:length(am)-1);
        pbm=bm(1:length(bm)-1);
        m1=max(abs(pam));
        m2=length(bm)-1;
        thisM=m1*m2;
    end
    
catch
    
    m1=max(abs(am));
    m2=length(bm);
    thisM=m1*m2;
end
thisM=double(thisM);

% pam=am(1:length(am)-1);
% pbm=bm(1:length(bm)-1);
% m1=max(abs(pam));
% m2=length(bm)-1;
% thisM=m1*m2;

outputdata{1}=caseType;
outputdata{2}=time;
outputdata{3}=epi;
outputdata{4}=thisL;
outputdata{5}=eta;
outputdata{6}=theK;
outputdata{7}=Boundedupdate;
outputdata{8}=thisM;
thec=thecoeffs(4);
outputdata{9}=thec;

end




%--------Done---------%


%save american;


%����h1����
function []=iter(n,q)%n is equal to  sizePV or lenoftarget or t, �����γɸ���������һ�����ѭ��Ƕ��
global saveii;%ÿ��ʹ�ö�Ҫ����
global degree;
global output;
global t;
global target;
if n==0
    if sum(saveii)<=degree
        z=1;
        for ii=1:t
            
            z=z*target(ii)^saveii(ii);
        end
        output(end+1)=z;
    end
    
else
    for i=degree:-1:0
        saveii(q)=i;
        iter(n-1,q+1);
    end
    
end

end


%end
