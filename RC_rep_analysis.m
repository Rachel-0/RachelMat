clear all
clc

% datafolder='/Volumes/RACHELSTICK/RachelCooperExperiment/Results_RC_rep/RC_rep_analysis/';
% datafolder='G:\RachelCooperExperiment\Results_RC_rep\RC_rep_analysis\';
datafolder='C:\Users\rcoopea\Desktop\Results_RC_rep\';
experiment='RC_rep';
ext='.mat';

% Qresults=Questionnaire_analysis;

% vi=find(Qresults.Qsubs(:,4)==1);
% IndianSubs=Qresults.Qsubs(vi,[5]);
% vi=find(Qresults.Qsubs(:,4)==2);
% EnglishSubs=Qresults.Qsubs(vi,[5]);

Ncondvoices=5;
Nvoices=10;
testnum=6;

subs=[2 3 4 6 9 10 11 12 13 14 15 16 18 19 20 21 23 24  27   111 115 116 120 121 122 126 144 145 200 201 888]; % 97 98 % 3 4 8 22 % subs missing blocks: 124 118 113 22 % subs with errors: 97 98
IndianSubs=[3 6 20 126 122 120 121 116 115 111 27 888     144 145 200]; % update these 3 lists % 97 % 3 %%%% 
EnglishSubs=[4 9 10 11 12 13 14 15 16 18 19 21 23 24  201]; % 98 % 4 8
Nsubs=length(subs);

for subnumber=1:Nsubs; 
   sub=subs(subnumber);
     if ~isempty(find(IndianSubs == sub))
        accent=1;
    elseif ~isempty(find(EnglishSubs == sub))
        accent=2;
    else accent=0;
    end
    for stimaccent=1:2
        name=[experiment '_' int2str(sub) '_' int2str(stimaccent) '_' int2str(testnum) ext];
        file=[datafolder name];
        load(file);
        
        totalcorrect=sum(Results.data(:,7));
        pcorrect=totalcorrect/length(Results.data);
        pcorrectboth(:,stimaccent)=pcorrect*100;
        if sub>=9
        vi=find(Results.data(:,7)==1);
        rtcorrect=Results.data(vi,9);
        rt=median(rtcorrect);
        rtboth(:,stimaccent)=rt;
        end
    end
    
    
    pcorrect1=pcorrectboth(:,1); % Indian
    pcorrect2=pcorrectboth(:,2); % English
    if sub>=9
        rt1=rtboth(:,1);
        rt2=rtboth(:,2);

        
        results.rt(sub,:)=[rt1, accent, 1, subnumber]; % Indian
        results.rt2(sub,:)=[rt2, accent, 2, subnumber]; % English
        
        results.header=['pcorrect' 'subaccent' 'stim accent' 'sub'];
        results.data(sub,:)= [ pcorrect1, accent, 1, subnumber]; % Indian 
        results.data2(sub,:)=[ pcorrect2, accent, 2, subnumber]; % English
    else
        results.header=['pcorrect' 'subaccent' 'stim accent' 'sub'];
        results.data(sub,:)= [ pcorrect1, accent, 1, subnumber]; % Indian 
        results.data2(sub,:)=[ pcorrect2, accent, 2, subnumber]; % English 
    end

  
end


adata=sortrows(results.data,2); % Indian
adata2=sortrows(results.data2,2); % English
aRT=sortrows(results.rt,2); % Indian
aRT2=sortrows(results.rt2,2); % English

results.databoth=[adata ; adata2];
results.rtboth=[aRT ; aRT2];

% adata(:,1:2)=[];
vi=find(results.databoth(:,4)==0);
results.databoth(vi,:)=[];
vi=find(results.rtboth(:,4)==0);
results.rtboth(vi,:)=[];
% stimaccent(1:Nsubs,:)=1;  % { ;'Indian' ;'Indian' ;'Indian' ;'Indian' ;'Indian' ;'Indian' ;'Indian' ;'English' ;'English' ;'English' ;'English'; 'English' ;'English' ;'English' ;'English'};
% stimaccent(Nsubs+1:Nsubs*2,:)=2;
% adata(:,3)=stimaccent;

% % % % fprintf ('\n    Percentage Correct Anova    \n');
% % % % anova=mixed_between_within_anova(results.databoth);
% % % % fprintf ('\n\n    Reation Time Anova    \n');
% % % % anova=mixed_between_within_anova(results.rtboth);

% [P,TABLE,STATS] = anovan(adata,accentgroup, 'model','full');
%  [P,TABLE,STATS] = anova2(adata,length(IndianSubs));
%   [H P CI stats]=ttest(results.data(:,2),results.data(:,3));
fprintf ('\n\n  * *** Finished *** *\n\n');