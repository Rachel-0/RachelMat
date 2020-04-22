function Qresults=Questionnaire_analysis

% load
% average
% make subject matrix

Qheaders={'sub ethnicity tongue UKy UKm age gender fluency date q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12 q13 q14 q15'};

% RC4&5_2_2119demographics.mat

Qfolder='F:\RachelCooperExperiment\Results_Questionnaire\';
% Qfolder='/Volumes/RACHELSTICK/RachelCooperExperiment/Results_Questionnaire/';
exp='RC4&5_';
suffix='demographics.mat';
IndianQs=[1 4 5 8 9 11 12];
EnglishQs=[2 3 6 7 10 13 14];

subs=[2 3 5 92 93 94 97 98 99 911 912 913 914 916]; %  96 915 % update
Nsubs=length(subs);

for subnumber=1:Nsubs; 
   sub=subs(subnumber);

file=[Qfolder exp int2str(sub) '_' suffix];

load(file);

DEMsubs(subnumber,:)=Resultsdemog(:,[1:9 24]);

if strcmp(DEMsubs(subnumber,10), '6')==1 | strcmp(DEMsubs(subnumber,10), '5')==1
    subaccent=1;
elseif strcmp(DEMsubs(subnumber,10), '1')==1 | strcmp(DEMsubs(subnumber,10), '2')==1
    subaccent=2;
else subaccent=0;
end


Indianans(subnumber,:)=Resultsdemog(:,IndianQs+10);
Englishans(subnumber,:)=Resultsdemog(:,EnglishQs);
Indianexp=cellfun(@str2num,Indianans, 'UniformOutput' , false);
Englishexp=cellfun(@str2num,Englishans, 'UniformOutput', false);
Indianexp=cell2mat(Indianexp);
Englishexp=cell2mat(Englishexp);

IndianexpM=mean(Indianexp(subnumber,:));
EnglishexpM=mean(Englishexp(subnumber,:));

Qsubs(subnumber,1:2)=[IndianexpM EnglishexpM];
Qsubs(subnumber,3)=EnglishexpM-IndianexpM;
Qsubs(subnumber,4)=subaccent;
Qsubs(subnumber,5)=sub;
% Qsubs(subnumber,1)=IndianexpM;
% Qsubs(subnumber,2)=EnglishexpM;
end

Qresults.Qsubs=Qsubs;
Qresults.DEMsubs=DEMsubs;
outputfile=[Qfolder 'Qsubs_' int2str(Nsubs)];

save(outputfile, 'Qresults');


fprintf '\n\n\n    * *** Questionnaire_analysis finished *** * \n\n\n';