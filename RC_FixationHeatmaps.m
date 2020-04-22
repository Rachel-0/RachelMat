function [x_in y_in x_out y_out] = RC_FixationHeatmaps(subnumber, behav, fix, bias, HitStim, MissStim);

% ----------------------- Learn and hit/miss trials ---------------------%

Linds = find(fix.data(:,5)<41);
fixreportL = fix.data(Linds,:); % learn trials

if bias==1;
    if subnumber==64
        stimcol = 14;
    else
        stimcol=14; %18 for sub 64
    end
elseif bias==2;
    stimcol=16;
end

% hit_size = length(HitStim);
% miss_size = length(MissStim);
stim = [HitStim ; MissStim];

% convert HitStims and MissStims into numbers
for i= 1:length(stim); % find all hit stims
    hits(i)=str2num(stim{i});
end
% for i= 1:length(MissStim); % find all missed stims
%     misses(i)=str2num(MissStim{i});
% end

% find the indices for hit and missed stim in the learn trials
% c=ismember(fixreportL(:,stimcol),misses); % c = a boolean matrix. true(row) if trialreportL(row) is included in misses(row). c is the same size as trialreportL.
% m=find(c(:,1)==1);
c=ismember(fixreportL(:,stimcol),hits);
h=find(c(:,1)==1);

Trials =fixreportL(h,:); % data for trials leading to hits
% MissedTrials=fixreportL(m,:); % data for trials leading to misses

% -------------------------- Group membership -------------------------- %

% % Split trials in to in and out groups depending on experiment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bias ==1 % RACE BIAS
    
    sh=num2str(fixreportL (:,stimcol)); % find which faces were Black or White from the hit trials. Look at 1st digit in stimulus code.
    [r,col]=size(fixreportL );
    for row = 1:r;
        if strcmp(sh(row,1),'1')==1;
            faceraceh(row)=1; % Black face
        elseif strcmp(sh(row,1),'2')==1;
            faceraceh(row)=2; % White face
        else
            if strcmp(sh(row,2),'1')==1; % if decimal point is after 1st rather than 2nd digit
                faceraceh(row)=1; % Black face
            elseif strcmp(sh(row,2),'2')==1;
                faceraceh(row)=2; % White face
            end
        end
    end
%     sm=num2str(MissedTrials(:,stimcol)); % find which faces were Black or White from the miss trials
%     [r,col]=size(MissedTrials);
%     for row = 1:r;
%         if strcmp(sm(row,1),'1')==1;
%             faceracem(row)=1; % Black face
%         elseif strcmp(sm(row,1),'2')==1;
%             faceracem(row)=2; % White face
%         else
%             if strcmp(sm(row,2),'1')==1; % if decimal point is after 1st rather than 2nd digit
%                 faceracem(row)=1; % Black face
%             elseif strcmp(sm(row,2),'2')==1;
%                 faceracem(row)=2; % White face
%             end
%         end
%     end
    blackindh=find(faceraceh==1);
    whiteindh=find(faceraceh==2);
%     blackindm=find(faceracem==1);
%     whiteindm=find(faceracem==2);
% %     
    trial_ingroupHits=fixreportL (whiteindh,:);
%     trial_ingroupMisses=MissedTrials(whiteindm,:);
    trial_outgroupHits=fixreportL (blackindh,:);
%     trial_outgroupMisses=MissedTrials(blackindm,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif bias ==2 % GROUP BIAS
    if fixreportL (1,11)==1; % match background colour with group (Essex or not)
        ini = find(fixreportL (:,9)==1);
        trial_ingroupHits = fixreportL (ini,:);
        outi = find(fixreportL (:,9)==2);
        trial_outgroupHits = fixreportL (outi,:);
    elseif fixreportL (1,11)==2; % match background colour with group (Essex or not)
        ini = find(fixreportL (:,9)==2);
        trial_ingroupHits = fixreportL (ini,:);
        outi = find(fixreportL (:,9)==1);
        trial_outgroupHits = fixreportL (outi,:);
    end
%     if MissedTrials(1,9)==1; % match background colour with group (Essex or not)
%         ini = find(MissedTrials(:,12)==1); % Essex
%         trial_ingroupMisses = MissedTrials(ini,:);
%         outi = find(MissedTrials(:,12)==2); % Anglia
%         trial_outgroupMisses = MissedTrials(outi,:);
%     elseif MissedTrials(1,9)==2; % match background colour with group (Essex or not)
%         ini = find(MissedTrials(:,12)==2);
%         trial_ingroupMisses = MissedTrials(ini,:);
%         outi = find(MissedTrials(:,12)==1);
%         trial_outgroupMisses = MissedTrials(outi,:);
%     end
end

% proportion_indices = [1];

% a = mean(trial_ingroupHits,1); % 26 columns each
% % a(proportion_indices)=a(proportion_indices)/hit_size; % for loop needed?
% b = mean(trial_outgroupHits,1);
% % b(proportion_indices)=b(proportion_indices)/hit_size;
% c = mean(trial_ingroupMisses,1);
% % c(proportion_indices)=c(proportion_indices)/miss_size;
% d = mean(trial_outgroupMisses,1);
% % d(proportion_indices)=d(proportion_indices)/miss_size;

% ---------- Build x and y vectors for all fixations (all subs per group, each cond) ---------- %

% AllSubsXY_in = zeros(2,2400); % prebuild vectors
% AllSubsXY_out = zeros(2,2400);
% AllSubsX_in = [];
% AllSubsX_out = [];
% AllSubsY_in = [];
% AllSubsY_out = [];

x_in = trial_ingroupHits(:,2);
x_out = trial_outgroupHits(:,2);
y_in = trial_ingroupHits(:,3);
y_out = trial_outgroupHits(:,3);




% if subnumber == 2 | subnumber ==3
%     Eyesub = [a 999 b 999 c 999 d]; % inHits outHits inMisses outMisses









%     Eyesub = [a b c d]; % inHits outHits inMisses outMisses
