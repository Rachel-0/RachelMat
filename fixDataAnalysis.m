function fixsub = fixDataAnalysis(subnumber, behav, fix, bias, HitStim, MissStim);

% find rows for in group, out group, hits, misses, learn trials
% average all columns across trial

fixsub = 1;

fixHeaders = {'RECORDING_SESSION_LABEL','CURRENT_FIX_DURATION','CURRENT_FIX_END','CURRENT_FIX_INDEX','CURRENT_FIX_INTEREST_AREAS','CURRENT_FIX_START','CURRENT_FIX_X','CURRENT_FIX_Y','TRIAL_FIXATION_TOTAL','TRIAL_INDEX','TRIAL_START_TIME','TrialCountLearn','TrialCountTest','backcolour_n','backcolour_n_t','essexcolourn','facegender','facegender_t','learnt','learnt_t','phase','stimnumber','stimnumber_t','CURRENT_FIX_INTEREST_AREA_INDEX';};

a = fix.textdata(2:end,4);
a(1,:)=[];
a=cell2mat(a);
fix.data(:,20) = a; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bias==1;
%     stimcol=18;
%     phasecol = 22;
elseif bias==2;
%     stimcol=21;
    phasecol = 16;
end

vi = find(fix.data(:,phasecol)==1);
fixL = fix.data(vi,:); % learn trials

a = fixL(:,5);
h=h'; % find all rows for hit trials (learn phase)
i = ismember(a,h,'rows');
vi = find(i(:,1)==1);
fixHits = fixL(vi,:);

% m=m'; % find all rows for miss trials (learn phase)
i = ismember(a,m,'rows');
vi = find(i(:,1)==1);
iaMisses = iaL(vi,:);

% Find all rows for in-group trials (hit and miss)
% Find all rows for out-group trials (hit and miss)

if bias ==1 % RACE BIAS
    
    sh=num2str(iaHits(:,stimcol));% find which faces were Black or White from the hit trials. Look at 1st digit in stimulus code.
    [r,col]=size(iaHits);
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
    sm=num2str(iaMisses(:,stimcol)); % find which faces were Black or White from the miss trials
    [r,col]=size(iaMisses);
    for row = 1:r;
        if strcmp(sm(row,1),'1')==1;
            faceracem(row)=1; % Black face
        elseif strcmp(sm(row,1),'2')==1;
            faceracem(row)=2; % White face
        else
            if strcmp(sm(row,2),'1')==1; % if decimal point is after 1st rather than 2nd digit
                faceracem(row)=1; % Black face
            elseif strcmp(sm(row,2),'2')==1;
                faceracem(row)=2; % White face
            end
        end
    end
    blackindh=find(faceraceh==1);
    whiteindh=find(faceraceh==2);
    blackindm=find(faceracem==1);
    whiteindm=find(faceracem==2);
    
    trial_ingroupHits=iaHits(whiteindh,:);
    trial_ingroupMisses=iaMisses(whiteindm,:);
    trial_outgroupHits=iaHits(blackindh,:);
    trial_outgroupMisses=iaMisses(blackindm,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif bias ==2; % UNIVERSITY BIAS
    
    if iaHits(1,14)==1; % match background colour with group (Essex or not)
        ini = find(iaHits(:,15)==1);
        trial_ingroupHits = iaHits(ini,:);
        outi = find(iaHits(:,15)==2);
        trial_outgroupHits = iaHits(outi,:);
    elseif iaHits(1,14)==2; % match background colour with group (Essex or not)
        ini = find(iaHits(:,15)==2);
        trial_ingroupHits = iaHits(ini,:);
        outi = find(iaHits(:,15)==1);
        trial_outgroupHits = iaHits(outi,:);
    end
    if iaMisses(1,14)==1; % match background colour with group (Essex or not)
        ini = find(iaMisses(:,15)==1); % Essex
        trial_ingroupMisses = iaMisses(ini,:);
        outi = find(iaMisses(:,15)==2); % Anglia
        trial_outgroupMisses = iaMisses(outi,:);
    elseif iaMisses(1,14)==2; % match background colour with group (Essex or not)
        ini = find(iaMisses(:,15)==2);
        trial_ingroupMisses = iaMisses(ini,:);
        outi = find(iaMisses(:,15)==1);
        trial_outgroupMisses = iaMisses(outi,:);
    end
end








