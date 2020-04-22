
% Run 1st

%%% TO DO:
%%% ONLY AVERAGE WITHIN ACCENT TYPE (across blocks)     *
%%% CORRECT BEHAVIOURAL RESPONSES? *DONT HAVE THEM*
%%% INTERPOLATE BAD ELECTRODES                          *
%%% AVERAGE ACROSS PARTICIPANTS (GRAND AVERAGE)         *
%%% FIND WAY TO DEAL WITH EYE ARTEFACTS FOR SUBS 1:10
% get attend and accent info from RC3_analysis.m        *


% Make create basic event list work
% make use of perm theshold


clear EEG;
clear ERP;
clear all;
clc;


directory = 'C:\';
cntfolder='RachelC\Study3\';
ev2folder='RachelC\';
eventlistfolder='erp_analysis\';
fileprefix='p';
cntsuffix='.cnt';
ev2suffix='.ev2';
eventlistsuffix='eventlist.txt';
Nblocks = 16;
textsuffix = '.txt';
EEGSetDir = ['C:\Users\rcoopea\Desktop\EEGAnalysis_classicERPS\EEGsets\'];
ICAWeightsDir = ['C:\Users\rcoopea\Desktop\EEGAnalysis_classicERPS\ICAweights\'];
Nsamples = 555;
Nchannels = 36; % after mastoid averaging

% % % %         Channel indices
% % % %         1 = FP1
% % % %         2 = FP2
% % % %         3 = F7
% % % %         4 = F3
% % % %         5 = FZ
% % % %         6 = F4
% % % %         7 = F8
% % % %         8 = FT7
% % % %         9 = FC3
% % % %         10 = FCZ
% % % %         11 = FC4
% % % %         12 = FT8
% % % %         13 = T7
% % % %         14 = C3
% % % %         15 = CZ
% % % %         16 = C4
% % % %         17 = T8
% % % %         18 = TP7
% % % %         19 = CP3
% % % %         20 = CPZ
% % % %         21 = CP4
% % % %         22 = TP8
% % % %         23 = P7
% % % %         24 = P3
% % % %         25 = PZ
% % % %         26 = P4
% % % %         27 = P8
% % % %         28 = O1
% % % %         29 = OZ
% % % %         30 = O2
% % % %         31 = M1
% % % %         32 = M2
% % % %         33 = VEOG
% % % %         34 = E1
% % % %         35 = E3
% % % %         36 = averageMastoids

subs = [1:40];
% subs = [1:6 8:14 16:40];
Nsubs = length(subs);
% subs = [  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40]; % % 1 2 3 4 5 6 7 8
% subs = [27 28   30 31 32 33 34 35 36 37 38 39 40 ]; % 15 8 % for debugging
% subs = [29];
% prebuild AllSubs matrices
% ALLSubs_ac1at1_DiffWaves = zeros(Nchannels,Nsamples,Nsubs);
% ALLSubs_ac1at2_DiffWaves = zeros(Nchannels,Nsamples,Nsubs);
% ALLSubs_ac2at1_DiffWaves = zeros(Nchannels,Nsamples,Nsubs);
% ALLSubs_ac2at2_DiffWaves = zeros(Nchannels,Nsamples,Nsubs);

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % initiate eeglab



subcount=0;
for sub=subs;
    fprintf(['\n\n\n****   ' int2str(sub)  '   ****\n\n\n']);
    subcount=subcount+1;
    fprintf(['\n\n\n' subcount '\n\n\n']);
    
    
    
    %                 %%%% Load 40 participants
    %
    EEG = pop_loadset('filename',['p' int2str(sub) 'AllBlocks_art_Interp_ICA_BC_dipfit_2_finefit_artICA_mast.set'],'filepath',EEGSetDir);
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    eeglab redraw;
    % -----------------------------------------------
    LatencyRange = 211:251; % latency1
%     LatencyRange = 251:326; % latency2
%     LatencyRange = 211:246; % latency3
%     LatencyRange = 246:311; % latency4
%     LatencyRange = 311:351; % latency5
    % -----------------------------------------------
    minL = min(LatencyRange);
    maxL = max(LatencyRange);
    minL_Time = EEG(1).times(minL);
    maxL_Time = EEG(1).times(maxL);
    LatencyRangeTimes = [minL_Time maxL_Time];
    MeasurementWindowMS = maxL_Time - minL_Time;
    MeasurementWindowSecs = MeasurementWindowMS / 1000;
    
    %     % Average mastoid channels creating new channel 36
    %     EEG = pop_eegchanoperator( EEG, {  'ch36 =ch31+ch32/2'} , 'ErrorMsg', 'popup', 'Warning', 'on' ); % GUI: 07-Aug-2014 20:48:18
    %     [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    
    % update bin columns used by ERPlab
    if sub == 7
        EEG.event(476) = [];
    elseif sub == 15;
        EEG.event(120) = [];
    end
    
    b=[EEG.event.AccentCond; EEG.event.AttendCond; EEG.event.PersonCond; EEG.event.bini];
    %     b1=EEG.event.AccentCond;
    %     b2=EEG.event.AttendCond;
    %     b3=EEG.event.PersonCond;
    %     b4=EEG.event.bini;
    %     b = [b1 b2 b3 b4];
    b=b';
    Aa=find(b(:,1)==1 & b(:,2)==1 & b(:,4)==1); % accent 1, attend 1, standard
    Ab=find(b(:,1)==1 & b(:,2)==1 & b(:,4)==2); % accent 1, attend 1, odd
    
    Ba=find(b(:,1)==1 & b(:,2)==2 & b(:,4)==1); % accent 1, attend 2, standard
    Bb=find(b(:,1)==1 & b(:,2)==2 & b(:,4)==2); % accent 1, attend 2, odd
    
    Ca=find(b(:,1)==2 & b(:,2)==1 & b(:,4)==1); % accent 2, attend 1, standard
    Cb=find(b(:,1)==2 & b(:,2)==1 & b(:,4)==2); % accent 2, attend 1, odd
    
    Da=find(b(:,1)==2 & b(:,2)==2 & b(:,4)==1); % accent 2, attend 2, standard
    Db=find(b(:,1)==2 & b(:,2)==2 & b(:,4)==2); % accent 2, attend 2, odd
    
    newbini = zeros(length(b),1);
    newbini(Aa)=1; % accent 1, attend 1, standard
    newbini(Ab)=2; % accent 1, attend 1, odd
    newbini(Ba)=3; % accent 1, attend 2, standard
    newbini(Bb)=4; % accent 1, attend 2, odd
    newbini(Ca)=5; % accent 2, attend 1, standard
    newbini(Cb)=6; % accent 2, attend 1, odd
    newbini(Da)=7; % accent 2, attend 2, standard
    newbini(Db)=8; % accent 2, attend 2, odd
    
    for row=1:length(newbini);
        EEG.event(row).bini = newbini(row);
        EEG.event(row).binlabel = newbini(row);
        ALLEEG(1).event(row).bini = newbini(row);
        ALLEEG(1).event(row).binlabel = newbini(row);
        ALLEEG(1).epoch(row).eventbini = newbini(row);
        ALLEEG(1).epoch(row).eventbinlabel = newbini(row);
    end
    
    EEG.EVENTLIST.nbin=8;
    
    eeglab redraw;
    erplab redraw;
    
    % Create ERPs
    % --------------------------------------------------------------------
    % average across trials for each sample, each channel within each condition
    
    % create matrices for each condition
    Ac1At1stan = EEG.data(:,:,Aa);
    Ac1At1odd = EEG.data(:,:,Ab);
    Ac1At2stan = EEG.data(:,:,Ba);
    Ac1At2odd = EEG.data(:,:,Bb);
    Ac2At1stan = EEG.data(:,:,Ca);
    Ac2At1odd = EEG.data(:,:,Cb);
    Ac2At2stan = EEG.data(:,:,Da);
    Ac2At2odd = EEG.data(:,:,Db);
    
    % Get indices for each attention task (not differences waves)
    At1 = EEG.data(:,:,[Aa;Ab;Ca;Cb]); % attention 1 indices
    At2 = EEG.data(:,:,[Ba;Bb;Da;Db]); % attention 2 indices
    
    % create matrices for each task, oddballs and standards seperately,
    % collapsed across accent
    
    % Average across trials per condition
    Ac1At1stan = mean(Ac1At1stan,3);
    Ac1At1odd = mean(Ac1At1odd,3);
    Ac1At2stan = mean(Ac1At2stan,3);
    Ac1At2odd = mean(Ac1At2odd,3);
    Ac2At1stan = mean(Ac2At1stan,3);
    Ac2At1odd = mean(Ac2At1odd,3);
    Ac2At2stan = mean(Ac2At2stan,3);
    Ac2At2odd = mean(Ac2At2odd,3);
    
    % Create oddball minus standard difference wave
    Ac1At1diff = Ac1At1odd - Ac1At1stan;
    Ac1At2diff = Ac1At2odd - Ac1At2stan;
    Ac2At1diff = Ac2At1odd - Ac2At1stan;
    Ac2At2diff = Ac2At2odd - Ac2At2stan;
    
    % build AllSubs matrices for each condition
    ALLSubs_ac1at1_DiffWaves(:,:,sub) = Ac1At1diff; % Indian voice, attend to video
    ALLSubs_ac1at2_DiffWaves(:,:,sub) = Ac1At2diff; % Indian voice, attend to voice
    ALLSubs_ac2at1_DiffWaves(:,:,sub) = Ac2At1diff; % English voice, attend to video
    ALLSubs_ac2at2_DiffWaves(:,:,sub) = Ac2At2diff; % English voice, attend to voice
    
    % ---- For plots ----
    plotting_ac1at1_DiffWaves2(:,:,sub)= ALLSubs_ac1at1_DiffWaves;
    plotting_ac1at2_DiffWaves2(:,:,sub)= ALLSubs_ac1at2_DiffWaves;
    plotting_ac2at1_DiffWaves2(:,:,sub)= ALLSubs_ac2at1_DiffWaves;
    ploting_ac2at2_DiffWaves2(:,:,sub)= ALLSubs_ac2at2_DiffWaves;
    
    % 
    plotting_Ac1At1stan(:,:,sub) = Ac1At1stan;
    plotting_Ac1At1odd(:,:,sub) = Ac1At1odd;
    plotting_Ac1At2stan(:,:,sub) = Ac1At2stan;
    plotting_Ac1At2odd(:,:,sub) = Ac1At2odd;
    plotting_Ac2At1stan(:,:,sub) = Ac2At1stan;
    plotting_Ac2At1odd(:,:,sub) = Ac2At1odd;
    plotting_Ac2At2stan(:,:,sub) = Ac2At2stan;
    plotting_Ac2At2odd(:,:,sub) = Ac2At2odd;
    
%     at1Fig = ALLSubs_ac1at1_DiffWaves;
%     at1Fig(:,:,41:80) = ALLSubs_ac2at1_DiffWaves;
%     at1Fig = mean(at1Fig,3);

% % standard vs oddball per attention condition
% at1_stan = Ac1At1stan;
% at1_stan(:,:,41:80) = Ac2At1stan;
% plotting_at1_stan = mean(at1_stan,3);
% 
% at1_odd = Ac1At1odd;
% at1_odd(:,:,41:80) = Ac2At1odd;
% plotting_at1_odd = mean(at1_odd,3);
% 
% at2_stan = Ac1At2stan;
% at2_stan(:,:,41:80) = Ac2At2stan;
% plotting_at2_stan = mean(at2_stan,3);
% 
% at2_odd = Ac1At2odd;
% at2_odd(:,:,41:80) = Ac2At2odd;
% plotting_at2_odd = mean(at2_odd,3);
    
    
    % rectify numbers so that they are absolute values (no sign)
    ALLSubs_ac1at1_DiffWaves= abs(ALLSubs_ac1at1_DiffWaves);
    ALLSubs_ac1at2_DiffWaves= abs(ALLSubs_ac1at2_DiffWaves);
    ALLSubs_ac2at1_DiffWaves= abs(ALLSubs_ac2at1_DiffWaves);
    ALLSubs_ac2at2_DiffWaves= abs(ALLSubs_ac2at2_DiffWaves);
    
    % get average amplitude between two latencies(samples).  (latencies partly
    % determined by grand average plots for each attention task condition)
    
    % sum
    Amp_Ac1At1diff = ALLSubs_ac1at1_DiffWaves(:,LatencyRange,sub); % create matrix containing data from the latency range of interest only
    SumAmp_Ac1At1diff = sum(Amp_Ac1At1diff,2); % average amplitude across these chosen samples
    
    Amp_Ac1At2diff = ALLSubs_ac1at2_DiffWaves(:,LatencyRange,sub); % create matrix containing data from the latency range of interest only
    SumAmp_Ac1At2diff = sum(Amp_Ac1At2diff,2); % average amplitude across these chosen samples
    
    Amp_Ac2At1diff = ALLSubs_ac2at1_DiffWaves(:,LatencyRange,sub); % create matrix containing data from the latency range of interest only
    SumAmp_Ac2At1diff = sum(Amp_Ac2At1diff,2); % average amplitude across these chosen samples
    
    Amp_Ac2At2diff = ALLSubs_ac2at2_DiffWaves(:,LatencyRange,sub); % create matrix containing data from the latency range of interest only
    SumAmp_Ac2At2diff = sum(Amp_Ac2At2diff,2); % average amplitude across these chosen samples
    
    
    
    % multiply by length of measurement window (seconds)
    AreaAmp_Ac1At1diff = SumAmp_Ac1At1diff.*MeasurementWindowSecs;
    AreaAmp_Ac1At2diff = SumAmp_Ac1At2diff.*MeasurementWindowSecs;
    AreaAmp_Ac2At1diff = SumAmp_Ac2At1diff.*MeasurementWindowSecs;
    AreaAmp_Ac2At2diff = SumAmp_Ac2At2diff.*MeasurementWindowSecs;
    
    % build AllSubs matrices for the means for each channel
    ALLSubs_AreaAmp_Ac1At1diff(sub,:) = AreaAmp_Ac1At1diff;
    ALLSubs_AreaAmp_Ac1At2diff(sub,:) = AreaAmp_Ac1At2diff;
    ALLSubs_AreaAmp_Ac2At1diff(sub,:) = AreaAmp_Ac2At1diff;
    ALLSubs_AreaAmp_Ac2At2diff(sub,:) = AreaAmp_Ac2At2diff;
    
    % % % ALL SUBS % % %
    % GRAND AVERAGE: average across all participants, all trials for each sample, each channel
    
    % TASK AVERAGE: average across all participants for each sample, each channel, each attention task condition
    
    
    
    
end

% end
AllSubs_AllConds_AreaAmplitude.Ac1At1diff = ALLSubs_AreaAmp_Ac1At1diff;
AllSubs_AllConds_AreaAmplitude.Ac1At2diff = ALLSubs_AreaAmp_Ac1At2diff;
AllSubs_AllConds_AreaAmplitude.Ac2At1diff = ALLSubs_AreaAmp_Ac2At1diff;
AllSubs_AllConds_AreaAmplitude.Ac2At2diff = ALLSubs_AreaAmp_Ac2At2diff;


% [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
% eeglab redraw;



fprintf('\n\n\n**** FINISHED ****\n\n\n');
fprintf('\n**** Enjoy your data! ****\n');

