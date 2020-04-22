function [DataOut,Table] = FixReport(file, bias, subnumber)
%Takes a fixation report file and outputs some trial-by-trial data in a matrix

% ARGUMENTS:
% file: string specifying a single full path filename
% DataOut: one row per trial. For N regions of interest, N columns of counts & totaltime. N columns indicating first two fixations (1=1st,2=2nd,3=1st&2nd)
% and N indicating time spent on those two fixations. Last three columns show: total fixation time, "missing time" and number of fixations
% Table: a matlab matrix of the input fixrport file..useful if using
% fixreport for condition analysis

%Notes:
% 1. "missing time" is the time not allocated to fixations between time=0 and the end of the last fixation. The actual trial length might be
% known and might be more useful for this (subtract total fixation time from trial length).
% 2. To calculate the average time per fixation from DataOut, divide the 2nd N columns by the 1st N columns. (see code below)
% Nroi=11;
% AveFix=DataOut(:,Nroi+1:Nroi*2)./DataOut(:,1:Nroi)

if nargin==0
    file = '/Users/skennett/Dropbox/shared/rachel_steffan/Eyetracking_Analysis_Files/results/10/fix.txt';
    bias = 1;
    subnumber = 10;
end

%importing information
h=1;
fmt='%n%n%n%n%s%n%n%n%n%n%n%n%n%n%n%n%n%n%n%s';
del='\t';
type='text';

%REGION DEFINITIONS
%each row is a region of interest. The values are the logical definitions
%of regions. so [1 -1 0 0 0 0 0 0 0 0 0] means the region is all parts of A
%that are not in B
% cd('F:\EyeTrack');
cd('/Users/skennett/Dropbox/shared/rachel_steffan/Eyetracking_Analysis_Files');
% cd('C:\Documents and Settings\skennett\My Documents\4Projects\Rachel\EyeTrack')
load roiDefined %loads variable roi

C=SteffanInput(file,h,fmt,type,del);
Table=zeros(length(C.data{1}),length(C.data)); %5 & 20 need special treatment
% MAKE A TABLE OF THE numerical DATA
for i=[1:4 6:19]
    Table(:,i)=cell2mat(C.data(1,i)); % C = 23 not 24 in sub 12, 14, 17 I think
end
for i=1:length(Table)
    p=str2num(C.data{1,20}{i,1}); %successfully converts the number-strings to numbers.
    if isempty(p)
        Table(i,20)=0;
    else
        Table(i,20)=p;
    end
end

iacol=20;
Table(:,5)=Table(:,iacol);
roifix=Table(:,5);
if Table(end,10)~=120
    subnumber
    crashfix;
end
% if bias==1 & subnumber <3
%     iacol = 20;
%     Table(:,5)=Table(:,iacol);
%     roifix = Table(:,5);
% elseif bias==1 & subnumber ==3
%     % CONVERT THE FIXATION ZONES TO A VARIABLE. PAD WITH ZEROS
%     D=C.data{1,5};
%     zones=zeros(length(D),length(roi));
%     for i=1:length(D)
%         E=textscan(D{i},'%s');
%         F=E{1};
%         p=length(F);
%         q=[];
%         for j=2:p
%             G=textscan(F{j},'%n%s');
%             q(j-1)=G{1};
%         end
%         zones(i,1:length(q))=q;
%     end
%     %padded zeros removed
%     a=max(find(max(zones)>0));
%     zones=zones(:,1:a);
%
%     %convert zones into roifix (the roi of each fixation)
%     for i=1:length(roi)
%         a=roi(i,:);
%         b=find(a==1);
%         c=find(a==-1);
%         d=zones==b;
%         d=max(d,[],2);
%         nots=logical(zeros(length(d),1));
%         for j=1:length(c)
%             e=zones==c(j);
%             e=max(e,[],2);
%             nots=or(nots,e);
%         end
%         e=and(d,~nots);
%         roifix(e)=i*1;
%     end
%     roifix=roifix';
%     Table(:,5)=roifix;
% elseif bias==2 & subnumber < 5
%
%     % CONVERT THE FIXATION ZONES TO A VARIABLE. PAD WITH ZEROS
%     D=C.data{1,5};
%     zones=zeros(length(D),length(roi));
%     for i=1:length(D)
%         E=textscan(D{i},'%s');
%         F=E{1};
%         p=length(F);
%         q=[];
%         for j=2:p
%             G=textscan(F{j},'%n%s');
%             q(j-1)=G{1};
%         end
%         zones(i,1:length(q))=q;
%     end
%     %padded zeros removed
%     a=max(find(max(zones)>0));
%     zones=zones(:,1:a);
%
%     %convert zones into roifix (the roi of each fixation)
%     for i=1:length(roi)
%         a=roi(i,:);
%         b=find(a==1);
%         c=find(a==-1);
%         d=zones==b;
%         d=max(d,[],2);
%         nots=logical(zeros(length(d),1));
%         for j=1:length(c)
%             e=zones==c(j);
%             e=max(e,[],2);
%             nots=or(nots,e);
%         end
%         e=and(d,~nots);
%         roifix(e)=i*1;
%     end
%     roifix=roifix';
%     Table(:,5)=roifix;
%
% elseif bias==1 & subnumber>3
%     iacol = 21;
%     Table(:,5)=Table(:,iacol);
%     roifix = Table(:,5);
% elseif bias ==2 & subnumber ==5
%     iacol = 24; %23
%     Table(:,5)=Table(:,iacol);
%     roifix = Table(:,5);
% elseif bias ==2 & subnumber ==6
%     iacol = 24; %23
%     Table(:,5)=Table(:,iacol);
%     roifix = Table(:,5);
% else % bias 2
%     iacol = 24;
%     Table(:,5)=Table(:,iacol);
%     roifix = Table(:,5);
% end

%Time analysis
%average fixation duration, total dwell time, time of first, time of first2
%convert roifix to a timeslice file?
fixstart=Table(:,6);
fixend=Table(:,3);
a=max([max(fixstart) max(fixend)]);

%need one row per trial. Output must be easily combinable. Need 11 columns for each stat [11*duration
%11*count 11*aveFixLength 11*1st2fix[1=1stonly; 2=2nd only; 3=both]
%11*1st2fixLength 1*missingTime 1*totalfixations
for i=1:Table(end,10) %goes to 120!! only want 40!
    b=find(Table(:,10)==i); % find all rows where trial number == i
    counts=zeros(1,length(roi)); % create matrices of the correct size and shape
    dwell=counts;
    earlyN=counts;
    earlyT=counts;
    AvgFixDuration=counts;
    for j=1:length(roi)% for 1:11 areas
        c=find(roifix(b)==j); %WARNING! INDICES ARE OF b **NOT** OF Table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        counts(j)=length(c); %number of fixations in roi
        dwell(j)=sum(Table(b(c),2)); % sum (table(trial(area,2))) i.e., sum of fixation durations for each trial each area
        AvgFixDuration(j) = dwell(j)/counts(j);
        if min(c)< 3 %region fixated in first two
            if find(c==1)
                if find(c==2)
                    earlyN(j)=3;
                    earlyT(j)=sum(Table(b(1:2),2));
                else
                    earlyN(j)=1;
                    earlyT(j)=Table(b(1),2);
                end
            else
                earlyN(j)=2;
                earlyT(j)=Table(b(2),2);
            end
        end
        
    end
    
    %     sum(Table(b,2)) % average fix duration in each trial. Rachel took out
    %     of DataOut and replaced with average fix duration for each area in
    %     each trial
    try
        DataOut(i,:)=[counts dwell earlyN earlyT AvgFixDuration (Table(b(end),3)-sum(Table(b,2))) length(b)]; % add columns for trial identifiers % 57 columns
        trialType(i,:) = Table(b(1),12:end);
        
    catch
        DataOut(i,:)=[counts dwell earlyN earlyT AvgFixDuration 0 length(b)];
        trialType(i,:)= zeros(1,9);
    end
end
stimcol = 19;
s=num2str(Table(:,stimcol));% find which faces were Black or White from the hit trials. Look at 1st digit in stimulus code.
[r,col]=size(Table);

for row = 1:r;
    if strcmp(s(row,1),'1')==1;
        facerace(row)=1; % Black face
    elseif strcmp(s(row,1),'2')==1;
        facerace(row)=2; % White face
    else
        if strcmp(s(row,2),'1')==1; % if decimal point is after 1st rather than 2nd digit
            facerace(row)=1; % Black face
        elseif strcmp(s(row,2),'2')==1;
            facerace(row)=2; % White face
        end
    end
end
Table(:,end) = facerace; % not sure this is column correct

% vi = find(Table(:,4)==1); % fix number column. trial type identifiers. Check correct for race and university biases.
%
% trialType = Table(vi,12:end); % 13 cols
DataOut=[DataOut trialType]; % 70 cols
areas = {'Leye' 'Reye' 'nose' 'mouth' 'chin' 'forehead' 'hair' 'neck' 'Lcheek' 'Rcheek' 'screen'}; % check this order
DataOutHeaders = {'Leyecount', 'Reyecount', 'nosecount', 'mouthcount', 'chincount', 'foreheadcount', 'haircount', 'neckcount', 'Lcheekcount', 'Rcheekcount', 'screencount' ,'Leyedwell', 'Reyedwell', 'nosedwell', 'mouthdwell', 'chindwell', 'foreheaddwell', 'hairdwell', 'neckdwell', 'Lcheekdwell', 'Rcheekdwell', 'screendwell', 'Leye12fix', 'Reye12fix', 'nose12fix', 'mouth12fix', 'chin12fix', 'forehead12fix', 'hair12fix', 'neck12fix', 'Lcheek12fix', 'Rcheek12fix', 'screen12fix',  'Leye12fixL', 'Reye12fixL', 'nose12fixL', 'mouth12fixL', 'chin12fixL', 'forehead12fixL', 'hair12fixL', 'neck12fixL', 'Lcheek12fixL', 'Rcheek12fixL', 'screen12fixL',  'Leyefixdur', 'Reyefixdur', 'nosefixdur', 'mouthfixdur', 'chinfixdur', 'foreheadfixdur', 'hairfixdur', 'neckfixdur', 'Lcheekfixdur', 'Rcheekfixdur', 'screenfixdur' ,'UnaccountedTime' ,'TotalFixationsPerTrial', 'TrialCountLearn', 'TrialCountTest' ,'backcolour_n', 'backcolour_n_t','essexcolourn','facegender','facegender_t','learnt','learnt_t','phase','stimnumber','stimnumber_t','CURRENT_FIX_INTEREST_AREA_INDEX'};
ending = ['FixReport analysed for subject:' num2str(subnumber)]


% Dwell time / counts for each area is the average fixation duration for
% each area
