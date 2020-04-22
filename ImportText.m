
path = 'C:\Users\rcoopea\Documents\MATLAB\Eye-tracking_analysis\';
extention = '.txt';
exp = 'p';%'eyegroup';


subs=[1 ]; %
Nsubs=length(subs);

for subnumber=1:Nsubs;
    sub=subs(subnumber);
    
    name = [path exp int2str(sub) extention];
    
    behav = importdata(name, '\t', 1);
    
    vi = behav.data(:,6)>0;
    Test = behav.data(vi,:);
    
    for j=1:length(Test);
        if Test(j,19)==109 & Test(j,10)==1 % Hit
            c(j)=1;
        elseif Test(j,19)==109 & Test(j,10)==2 % False alarm
            c(j)=2;
        elseif Test(j,19)==122 & Test(j,10)==1 % Miss
            c(j)=3;
        elseif Test(j,19)==122 & Test(j,10)==2 % Correct rejection
            c(j)=4;
        else
            c(j)='x';
        end
    end
    
    Test(:,20)=c; % Frequencies of hits, misses, correct rejections and false alarms for all test stimuli.
    
    vi = Test(:,20)==1;
    HitStim = Test(vi,8); % stimuli numbers for hit stim
    vi = Test(:,20)==3;
    MissStim = Test(vi,8); % stimuli numbers for missed stim
    
    
    onei=find(Test(:,16)==1);
    twoi=find(Test(:,16)==2);
    Group1=Test(onei,:);
    Group2=Test(twoi,:);
    
    PresentT = length(Test)/4; % 20
    AbsentT = length(Test)/4; % 20
    
    onefreq = Group1(:,20);
    twofreq = Group2(:,20);
    oneh(subnumber,:) = hist(onefreq,4); % frequencies for group 1 stims
    twoh(subnumber,:) = hist(twofreq,4); % frequencies for group 2 stims
    
    if behav.data(1,4)==1 % match background colour with group (Essex or not)
        InGroup = oneh;
        OutGroup = twoh;
    else behav.data(1,4)==2
        InGroup = twoh;
        OutGroup = oneh;
    end
    
    InpHit = InGroup(subnumber,1)/PresentT; % Hit rate
    InpFA = InGroup(subnumber,2)/PresentT; % FA rate
    [IndPrime, InC] = dprime(InpHit,InpFA,PresentT,AbsentT); % InGroup i.e., uni of Essex
    
    OutpHit = OutGroup(1)/PresentT;
    OutpFA = OutGroup(2)/PresentT;
    [OutdPrime, OutC] = dprime(OutpHit,OutpFA,PresentT,AbsentT); % OutGroup i.e., Anglia
    
    
    
    ResultsBehav.data(subnumber,:) = [sub InpHit InpFA IndPrime InC OutpHit OutpFA OutdPrime OutC];
    ResultsBehav.headers = ['sub InpHit InpFA IndPrime InC OutpHit OutpFA OutdPrime OutC'];
    ResultsBehav.HitStim = HitStim ;
    ResultsBehav.MissStim = MissStim ;
    
    
    
end

 fprintf ('\n\n  * *** Finished *** *\n\n'); % Prints to the command window when done
 
 