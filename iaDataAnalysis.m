function iasub = iaDataAnalysis(subnumber, behav, ia, bias, fixname);

% m = trial numbers leading to misses
% h = trial numbers leading to hits

% find rows for in group, out group, hits, misses, learn trials
% average all columns across trial

% iaHeaders = {'TRIAL_INDEX','RECORDING_SESSION_LABEL','IA_AREA','IA_AVERAGE_FIX_PUPIL_SIZE','IA_DWELL_TIME','IA_DWELL_TIME_%','IA_FIXATION_%','IA_FIXATION_COUNT','IA_ID','IA_RUN_COUNT','TRIAL_DWELL_TIME','TRIAL_FIXATION_COUNT','TRIAL_START_TIME','essexcolourn','backcolour_n','backcolour_n_t','TrialCountLearn','TrialCountTest','learnt','learnt_t','stimnumber','stimnumber_t','facegender','facegender_t','phase'};

stimcol = 18;
lastcol = 22;

vi = find(ia.data(:,1)<41);
iaL = ia.data(vi,:); % Learn trials %

% normalise ia size.
screensize = 786432; % L * W in pixels
iaNorm = (ia.data(:,3)/screensize)*100; % proportion of screen each area takes up
% ia.data(:,lastcol+1) = iaNorm;

iaNorm(:,2) = ia.data(:,9); % ia label
iaNorm(:,3) = ia.data(:,1); % trial label

for trial = 1:120
    trialind = find(iaNorm(:,3)==trial);
    sizetrial = iaNorm(trialind,:)'; % get each row as a trial and each ia as a column.
    if length(sizetrial)<11 % if there is no hair area, add an area with size 0 in column 7 and shift others along.
        sizetrial(:,8:11)=sizetrial(:,7:end);
        sizetrial(:,7)=0;
    elseif length(sizetrial)==12 % if there are two hair areas, sum them together and shift other areas back.
        aind=find(sizetrial(2,:)==7);
        a=sizetrial(1,aind);
        b=sum(a);
        sizetrial(:,8:11)=sizetrial(:,9:end);
        sizetrial(1,7)=b;
        sizetrial(:,12)=[];
        %  deal with missing ias and extra ias for some participants
    elseif length(sizetrial) ==13
        aind=find(sizetrial(2,:)==7);
        a=sizetrial(1,aind);
        b=sum(a);
        sizetrial(:,8:11)=sizetrial(:,10:end);
        sizetrial(1,7)=b;
        sizetrial(:,12:end)=[];
    end
    IAsize(trial,:)=sizetrial(1,:);
end


areas = {'Leye' 'Reye' 'nose' 'mouth' 'chin' 'forehead' 'hair' 'neck' 'Lcheek' 'Rcheek' 'screen'}; % check this order


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file = fixname;

DataOut=FixReport(file, bias, subnumber); % fixreport.m for each subject number % each row is a trial

d=1:120;
d=d';
DataOut=[DataOut d];

DataOutL = DataOut(1:40,:);

% % % d = DataOutL(:,71); % find all rows for hit trials (learn phase)
% % % i = ismember(d,h,'rows');
% % % vi = find(i(:,1)==1);
% % % iaHits = DataOutL(vi,:);
% % %
% % % % find all rows for miss trials (learn phase)
% % % i = ismember(d,m,'rows');
% % % vi = find(i(:,1)==1);
% % % iaMisses = DataOutL(vi,:);

% Find all rows for in-group trials (hit and miss)
% Find all rows for out-group trials (hit and miss)

stimcol=64;
sh=num2str(DataOutL(:,stimcol));% find which faces were Black or White. Look at 1st digit in stimulus code.
[r,col]=size(DataOutL);
for row = 1:r;
    if strcmp(sh(row,1),'1')==1;
        facerace(row)=1; % Black face
    elseif strcmp(sh(row,1),'2')==1;
        facerace(row)=2; % White face
    else
        if strcmp(sh(row,2),'1')==1; % if decimal point is after 1st rather than 2nd digit
            facerace(row)=1; % Black face
        elseif strcmp(sh(row,2),'2')==1;
            facerace(row)=2; % White face
        end
    end
end


%%Steffan here (manual dealing with in and out group labels
blackind=find(facerace==1);
whiteind=find(facerace==2);

ia_White=DataOutL(whiteind,:);
ia_Black=DataOutL(blackind,:);
Wh_sizeind = IAsize(whiteind,:);
Bl_sizeind = IAsize(blackind,:);

%STEFFAN TO ADD HITS/MISSES CATEGORIES?
% ia_WhiteHits=DataOutL(whiteindh,:);
% ia_BlackHits=DataOutL(blackindh,:);
% WhHits_sizeind = IAsize(whiteindh,:);
% BlHits_sizeind = IAsize(blackindh,:);

% remove trial identifier info leaving only eye-tracking data. from 71 to 57 columns
ia_Black(:,58:end)=[];
ia_White(:,58:end)=[];

for i = 1:size(ia_White,1) % learn trials
    for area = 1:11
        counts=ia_White(:,1:11);
        %counts(i,area)=counts(i,area)/Wh_sizeind(i,area);%normalises by area....get rid
        dwell = ia_White(:,12:23);
        %dwell(i,area)=dwell(i,area)/Wh_sizeind(i,area);
    end
end
for i = 1:size(ia_Black,1) % learn trials
    for area = 1:11
        counts=ia_Black(:,1:11);
        %counts(i,area)=counts(i,area)/Bl_sizeind(i,area);
        dwell = ia_Black(:,12:23);
        %dwell(i,area)=dwell(i,area)/Bl_sizeind(i,area);
    end
end

if size(ia_Black,1)>1; %averaging
    Bl = mean(ia_Black);
elseif size(ia_Black,1)==0;
    Bl = zeros(1,57);
else
    Bl = ia_Black;
end

if size(ia_White,1)>1; %averaging
    Wh = mean(ia_White);
elseif size(ia_White,1)==0;
    Wh = zeros(1,57);
else
    Wh = ia_White;
end

column_selection=[1:4, 12:15 34:37]; %only really want ia=1,2,3,4 only want counts, dwell and early stuff.

Bl = Bl(1,column_selection);
Wh = Wh(1,column_selection);

%iasub.headers = {'Leyecount', 'Reyecount', 'nosecount', 'mouthcount', 'chincount', 'foreheadcount', 'haircount', 'neckcount', 'Lcheekcount', 'Rcheekcount', 'screencount' ,'Leyedwell', 'Reyedwell', 'nosedwell', 'mouthdwell', 'chindwell', 'foreheaddwell', 'hairdwell', 'neckdwell', 'Lcheekdwell', 'Rcheekdwell', 'screendwell', 'Leye12fix', 'Reye12fix', 'nose12fix', 'mouth12fix', 'chin12fix', 'forehead12fix', 'hair12fix', 'neck12fix', 'Lcheek12fix', 'Rcheek12fix', 'screen12fix',  'Leye12fixL', 'Reye12fixL', 'nose12fixL', 'mouth12fixL', 'chin12fixL', 'forehead12fixL', 'hair12fixL', 'neck12fixL', 'Lcheek12fixL', 'Rcheek12fixL', 'screen12fixL',  'Leyefixdur', 'Reyefixdur', 'nosefixdur', 'mouthfixdur', 'chinfixdur', 'foreheadfixdur', 'hairfixdur', 'neckfixdur', 'Lcheekfixdur', 'Rcheekfixdur', 'screenfixdur' ,'UnaccountedTime' ,'TotalFixationsPerTrial','Leyecount', 'Reyecount', 'nosecount', 'mouthcount', 'chincount', 'foreheadcount', 'haircount', 'neckcount', 'Lcheekcount', 'Rcheekcount', 'screencount' ,'Leyedwell', 'Reyedwell', 'nosedwell', 'mouthdwell', 'chindwell', 'foreheaddwell', 'hairdwell', 'neckdwell', 'Lcheekdwell', 'Rcheekdwell', 'screendwell', 'Leye12fix', 'Reye12fix', 'nose12fix', 'mouth12fix', 'chin12fix', 'forehead12fix', 'hair12fix', 'neck12fix', 'Lcheek12fix', 'Rcheek12fix', 'screen12fix',  'Leye12fixL', 'Reye12fixL', 'nose12fixL', 'mouth12fixL', 'chin12fixL', 'forehead12fixL', 'hair12fixL', 'neck12fixL', 'Lcheek12fixL', 'Rcheek12fixL', 'screen12fixL',  'Leyefixdur', 'Reyefixdur', 'nosefixdur', 'mouthfixdur', 'chinfixdur', 'foreheadfixdur', 'hairfixdur', 'neckfixdur', 'Lcheekfixdur', 'Rcheekfixdur', 'screenfixdur' ,'UnaccountedTime' ,'TotalFixationsPerTrial','Leyecount', 'Reyecount', 'nosecount', 'mouthcount', 'chincount', 'foreheadcount', 'haircount', 'neckcount', 'Lcheekcount', 'Rcheekcount', 'screencount' ,'Leyedwell', 'Reyedwell', 'nosedwell', 'mouthdwell', 'chindwell', 'foreheaddwell', 'hairdwell', 'neckdwell', 'Lcheekdwell', 'Rcheekdwell', 'screendwell', 'Leye12fix', 'Reye12fix', 'nose12fix', 'mouth12fix', 'chin12fix', 'forehead12fix', 'hair12fix', 'neck12fix', 'Lcheek12fix', 'Rcheek12fix', 'screen12fix',  'Leye12fixL', 'Reye12fixL', 'nose12fixL', 'mouth12fixL', 'chin12fixL', 'forehead12fixL', 'hair12fixL', 'neck12fixL', 'Lcheek12fixL', 'Rcheek12fixL', 'screen12fixL',  'Leyefixdur', 'Reyefixdur', 'nosefixdur', 'mouthfixdur', 'chinfixdur', 'foreheadfixdur', 'hairfixdur', 'neckfixdur', 'Lcheekfixdur', 'Rcheekfixdur', 'screenfixdur' ,'UnaccountedTime' ,'TotalFixationsPerTrial','Leyecount', 'Reyecount', 'nosecount', 'mouthcount', 'chincount', 'foreheadcount', 'haircount', 'neckcount', 'Lcheekcount', 'Rcheekcount', 'screencount' ,'Leyedwell', 'Reyedwell', 'nosedwell', 'mouthdwell', 'chindwell', 'foreheaddwell', 'hairdwell', 'neckdwell', 'Lcheekdwell', 'Rcheekdwell', 'screendwell', 'Leye12fix', 'Reye12fix', 'nose12fix', 'mouth12fix', 'chin12fix', 'forehead12fix', 'hair12fix', 'neck12fix', 'Lcheek12fix', 'Rcheek12fix', 'screen12fix',  'Leye12fixL', 'Reye12fixL', 'nose12fixL', 'mouth12fixL', 'chin12fixL', 'forehead12fixL', 'hair12fixL', 'neck12fixL', 'Lcheek12fixL', 'Rcheek12fixL', 'screen12fixL',  'Leyefixdur', 'Reyefixdur', 'nosefixdur', 'mouthfixdur', 'chinfixdur', 'foreheadfixdur', 'hairfixdur', 'neckfixdur', 'Lcheekfixdur', 'Rcheekfixdur', 'screenfixdur' ,'UnaccountedTime' ,'TotalFixationsPerTrial'};
iasub = [Bl Wh];







