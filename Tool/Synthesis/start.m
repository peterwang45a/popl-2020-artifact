clc
clear all


example=input('Please enter the Example Name:','s')
name=example;


inputfile1=[name,'1.txt'];
inputfile2=[name,'2.txt'];
inputfile3=[name,'config.txt'];
%outputfile1=[name,'_output.mat'];
outputfile2=[name,'_log.txt'];


diary(outputfile2);
diary on;
%try
%     copyfile(['../Inputs/',inputfile1]);
%     copyfile(['../Inputs/',inputfile2]);
%     copyfile(['../Inputs/',inputfile3]);
    
    %     diary(outputfile2);
    %     diary on;
%catch
    
%end

% inputfile1=['../Inputs/',inputfile1];
% inputfile2=['../Inputs/',inputfile2];
% inputfile3=['../Inputs/',inputfile2];
if exist(inputfile1,'file')==0
    fprintf('No such Example Name: %s.\n',name);
    disp('Please try again and enter the correct name.');
    diary off;
    %movefile(outputfile2,'../Outputs');
    %exit
else

[outputdata] = synthesis(inputfile1,inputfile2,inputfile3);
%save(outputfile1);
fprintf('Example:%s\n',name);
disp('Analysis Finished Successfully.');
fprintf('Epsilon = %f\n',outputdata{3});
fprintf('L= %f\n',outputdata{4});
%disp('L = 1');

if outputdata{4}==1 && strcmp(outputdata{1},'affine')
    fprintf('Sensitivity Type: Non-expansive Expected %s-sensitive\n\n',outputdata{1});
end

if outputdata{4}>1 && outputdata{4}<exp(3*outputdata{3}^2/(8*outputdata{9}^2))
    fprintf('Sensitivity Type: Expansive Expected %s-sensitive\n\n','affine');
end

table2={'mini' 'rdwalk' 'vprdwalk' 'prspeed' 'vrace' 'ad2D' 'vad1D' 'american'};
table3={'vmini' 'single' 'double' 'vrdwalk' 'prdwalk' 'vprspeed' 'race' 'simple' 'pollutant' 'vad2D' 'ad1D' 'vamerican'};
flag=1;
for i=1:length(table2)
    if strcmp(table2{i},name)
       disp('Data Generated for Table 2 in the paper:');
       flag=2;
    end
end


for i=1:length(table3)
    if strcmp(table3{i},name)
       fprintf('Sensitivity Type: Non-expansive Expected %s-sensitive\n\n',outputdata{1});
       disp('Data Generated for Table 3 in the paper:');
       %fprintf('Sensitivity Type: Expected %s-sensitivity\n',outputdata{1});
       flag=3;
    end
end


if flag==1
    disp('Data Generated:');
end


fprintf('Runtime = %f seconds\n',outputdata{2});
fprintf('eta(b)=');
disp(outputdata{5});
fprintf('K = %f\n',outputdata{6});
fprintf('d = %f\n',outputdata{7});
fprintf('M = %f\n',outputdata{8});
if ~(outputdata{4}==1 && strcmp(outputdata{1},'affine'))
   fprintf('c = %f\n',outputdata{9});
end

diary off;

% movefile(outputfile1,'../Outputs');
% movefile(outputfile2,'../Outputs');

%delete *.txt
% movefile(inputfile1,'../Intputs');
% movefile(inputfile2,'../Intputs');
% movefile(inputfile3,'../Intputs');
end

%exit
