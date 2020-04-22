%Top level m-file that extracts key information from each participants
%report files.
%Each participant's data needs to be first extracted into a fixation
%report, an interest area report, and a trial report.
% initially written by Rachel Cooper
% modified by Steffan Kennett

%TO DO
%make all race and hit/miss allocations in this file as a single two-column
%variable. Pass this on to the two/three functions as needed

clear all
clc

%%%%% Bias switch %%%%%
bias = 1; % 1 = race bias, 2 = uni bias
%%%

% input paths
path='/Users/skennett/Dropbox/Shared/Rachel_Steffan/Eyetracking_Analysis_Files/'; %Steffan's Macbook

% output paths
outputBehav = [path 'AllSubsBehav.mat'];
outputEye = [path 'AllSubsFix.mat'];
outputEye2 = [path 'AllSubsIA.mat'];

% path and file building
resultspath = 'results/Master dissertation/';
fprintf ('Race bias');
macpcfolder='/'; %'\'
BehavExtention = 'RESULTS_FILE_TEST.txt';
fixExtention = 'fix.txt';
trialExtention = 'trial.txt';
iaExtention = 'ia.txt';

%Experiment-level settings
subs=[01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126];
% subs=05;%:126;
Nsubs=length(subs);
nlearnphase=40; %trials in learning phase
ntestphase=80; %trials is test phase

for subnumber=1:Nsubs; % subject loop
    clear C Table trialreport
    sub=subs(subnumber);
    
    if sub<10 %Depends how folders named
        behavname = [path resultspath '0' int2str(sub) macpcfolder BehavExtention];
        fixname = [path resultspath '0' int2str(sub) macpcfolder fixExtention];
        trialname = [path resultspath '0' int2str(sub) macpcfolder trialExtention];
        ianame = [path resultspath '0' int2str(sub) macpcfolder iaExtention];
    else
        behavname = [path resultspath int2str(sub) macpcfolder BehavExtention];
        fixname = [path resultspath int2str(sub) macpcfolder fixExtention];
        trialname = [path resultspath int2str(sub) macpcfolder trialExtention];
        ianame = [path resultspath int2str(sub) macpcfolder iaExtention];
    end
    
    behav = importdata(behavname, '\t', 1); %imports Results_File_Test.txt
    ia = importdata(ianame, '\t',1); %imports ia.txt report file

    %dealing with dots inputs
    h=1;
    del='\t';
    type='text';
    if length(behav.data)~=120
        sub
        crashbehav;
    end
    %trial.txt report data extraction
    fmt='%n%n%s%s%n%n%n%s%n%n%n%n%n%n%n%n%n%n%s%s%n%n%n';
    C=SteffanInput(trialname,h,fmt,type,del);
    Table=zeros(length(C.data{1}),length(C.data));
    for i=[1:2 5:7 9:18 21:23]
        Table(:,i)=cell2mat(C.data(1,i));
    end
    for i=1:length(Table)
        p=str2num(C.data{1,4}{i,1}); %successfully converts the number-strings to numbers.
        q=str2num(C.data{1,3}{i,1});
        r=str2num(C.data{1,8}{i,1});
        p2=str2num(C.data{1,19}{i,1});
        p3=str2num(C.data{1,20}{i,1});
        if isempty(p)
            Table(i,4)=0;
        else
            Table(i,4)=p;
        end
        if isempty(q)
            Table(i,3)=0;
        else
            Table(i,3)=q;
        end
        if isempty(r)
            Table(i,8)=0;
        else
            Table(i,8)=r;
        end
        if isempty(p2)
            Table(i,19)=0;
        else
            Table(i,19)=p2;
        end
        if isempty(p3)
            Table(i,20)=0;
        else
            Table(i,20)=p3;
        end
    end
    trialreport.data=Table; %This makes ia.data correct. What about text?
    if length(trialreport.data)~=120
        sub
        crashtrialreport;
    end
    iasub=zeros(length(subs),1);
    if ia.data(end,1)~=120
        iasub(sub)=1;
        fmt='%n%n%n%s%n%s%s%n%n%n%n%n%n%n%n%n%n%n%n%n%n';
        C=SteffanInput(ianame,h,fmt,type,del);
        Table=zeros(length(C.data{1}),length(C.data)); %only 5 needs special treatment
        for i=[1:3 8:21]
            Table(:,i)=cell2mat(C.data(1,i)); % C = 23 not 24 in sub 12, 14, 17 I think
        end
        for i=1:length(Table)
            p=str2num(C.data{1,4}{i,1}); %successfully converts the number-strings to numbers.
            q=str2num(C.data{1,6}{i,1});
            r=str2num(C.data{1,7}{i,1});
            if isempty(p)
                Table(i,4)=0;
            else
                Table(i,4)=p;
            end
            if isempty(q)
                Table(i,4)=0;
            else
                Table(i,4)=q;
            end
            if isempty(r)
                Table(i,4)=0;
            else
                Table(i,4)=r;
            end
        end
        ia.data=Table; %This makes ia.data correct. What about text?
        if ia.data(end,1)~=120
            sub
            crashia
        end
    end
    
    %%%% HIT vs MISS analysis
    vi=(nlearnphase+1):(nlearnphase+ntestphase);
    Test = behav.data(vi,:); % find all the test trials
    Teststrings=behav.textdata(42:121,[5 7 8]); % [stimcode(test) learntORnot response]
        
    for j=1:length(Teststrings); % for each test trial compare keyboard response with whether the face appeared at learn(1) or was a distractor(2)
        if     strcmp(Teststrings{j,3},'122')==1 & strcmp(Teststrings{j,2},'1')==1; % Hit
            c(j)=1;
        elseif strcmp(Teststrings{j,3},'122')==1 & strcmp(Teststrings{j,2},'2')==1; % False alarm
            c(j)=2;
        elseif strcmp(Teststrings{j,3},'109')==1 & strcmp(Teststrings{j,2},'1')==1; % Miss
            c(j)=3;
        elseif strcmp(Teststrings{j,3},'109')==1 & strcmp(Teststrings{j,2},'2')==1; % Correct rejection
            c(j)=4;
        else
            c(j)=999;
        end
    end
    
    hitvi = c(:,:)==1;
    HitStim = Teststrings(hitvi,1); % stimuli numbers for hit stim
    missvi = c(:,:)==3;
    MissStim = Teststrings(missvi,1); % stimuli numbers for missed stim
    
    blackTest=Test(:,4)==1; % Black
    whiteTest=Test(:,4)==2; % White
    c1=c(blackTest); %Hit Miss etc for Black faces
    c2=c(whiteTest); %Hit Miss etc for White faces

    PresentT = length(Test)/4; % Needed to get hit rate before passing to the dprime function
    AbsentT = length(Test)/4; %
    
%     Consider using this. Currently not in output
    blackh = hist(c1,4); % frequencies for group 1 stims
    whiteh = hist(c2,4); % frequencies for group 2 stims
    
    blackgroup=blackh; % black
    whitegroup=whiteh; % white
    
    WhitepHit = whitegroup(1)/PresentT; % Hit rate
    WhitepFA = whitegroup(2)/PresentT; % FA rate
    [WhitedPrime, WhiteC] = dprime(WhitepHit,WhitepFA,PresentT,AbsentT);
    
    BlackpHit = blackgroup(1)/PresentT;
    BlackpFA = blackgroup(2)/PresentT;
    [BlackdPrime, BlackC] = dprime(BlackpHit,BlackpFA,PresentT,AbsentT); 
    
    ResultsBehav.data(subnumber,:) = [sub BlackpHit BlackpFA BlackdPrime BlackC WhitepHit WhitepFA WhitedPrime WhiteC];
    ResultsBehav.headers = ['sub BlackpHit BlackpFA BlackdPrime BlackC WhitepHit WhitepFA WhitedPrime WhiteC'];
    ResultsBehav.HitStim = HitStim ; %won't work since overwritten by each new participant!
    ResultsBehav.MissStim = MissStim ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Eye tracking analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [Eyesub ] = EyeDataAnalysis(sub, behav, trialreport, bias, HitStim, MissStim);
    iasub = iaDataAnalysis(sub, behav, ia, bias, fixname);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    EyeData(subnumber,:)=[sub Eyesub];
    IAData(subnumber,:) = [sub iasub];
end


save(outputBehav, 'ResultsBehav');
save(outputEye, 'EyeData'); %7 columns: Sub, Face race * [Average Fixtation time, Number of fixations, Total dwell time]
save(outputEye2, 'IAData'); %25 sub, then just ia=1:4 and just counts, dwell and earlydwellT

fprintf ('\n\n  * *** Finished *** *\n\n'); % Prints to the command window when done

