% Figure making

% Use to inform LatencyRange in AverageMastoids.m
% Uses data from AverageMastoids.m (don't clear variables)
% Made to suit a sample rate of 500 samples per second

close all;

% ------------------ Switches ---------------------------
topography = 0; % Topography (1) or waveforms (0)
topoplotRange = 1; % Plot topography at the mean of LatencyRange (1) or one time point (0)

% which waveforms would you like to plot?
WaveOption = 6;
% 1 = Attention 1 - Oddball, standard, difference *
% 2 = Attention 2 - Oddball, standard, difference *
% 3 = Attention 1 - accent1 difference, accent2 difference
% 4 = Attention 2 - accent1 difference, accent2 difference
% 5 = Attention 1 - accent1  difference, accent2 difference - IndianSubs and EnglishSubs
% 6 = Attention 1 - accent1  difference, accent2 difference - IndianSubs and EnglishSubs


% ---------------------- Define constants -------------------------
ChannelOfInterest = 15; % number of channel

LatencyRange = 211:251; % latency1
%     LatencyRange = 251:326; % latency2
%     LatencyRange = 211:246; % latency3
%     LatencyRange = 246:311; % latency4
%     LatencyRange = 311:351; % latency5

minL = min(LatencyRange);
maxL = max(LatencyRange);
minL_Time = EEG(1).times(minL);
maxL_Time = EEG(1).times(maxL);
LatencyRangeTimes = [minL_Time maxL_Time];
MeasurementWindowMS = maxL_Time - minL_Time;
MeasurementWindowSecs = MeasurementWindowMS / 1000;
xfig = EEG.times;
electrodes = {'FP1' 'FP2' 'F7' 'F3' 'FZ' 'F4' 'F8' 'FT7' 'FC3' 'FCZ' 'FC4' 'FT8' 'T7' 'C3' 'CZ' 'C4' 'T8' 'TP7' 'CP3' 'CPZ' 'CP4' 'TP8' 'P7' 'P3' 'PZ' 'P4' 'P8' 'O1' 'OZ' 'O2' 'M1' 'M2' 'VEOG' 'E1' 'E3' 'Average Mastoids'};
SubAccent = [1 1 2 1 2 2 2 1 2 1 1 2 2 2 1 1 1 1 2 1 1 2 2 2 1 2 2 2 2 2 2 1 2 1 2 1 1 1 1 1 ]';
IndianSubInds = find(SubAccent(:,1)==1);
EnglishSubInds = find(SubAccent(:,1)==2);
% ----------------- Variable preparation -------------------------

% Get Indian and English sub indices to use below

% get oddball and standard waves for each attention condition i.e.,
% collapse accross voice accent type
at1_stan = plotting_Ac1At1stan;
at1_stan(:,:,41:80) = plotting_Ac2At1stan;
plotting_at1_stan = mean(at1_stan,3);
plotting_at1_stan = plotting_at1_stan(ChannelOfInterest,:);% Attention 1 standard
plotting_at1_Indian_stan = mean(plotting_Ac1At1stan,3);
plotting_at1_Indian_stan = plotting_at1_Indian_stan(ChannelOfInterest,:);%% Attention 1 standard - Indian voices
plotting_at1_English_stan = mean(plotting_Ac2At1stan,3);
plotting_at1_English_stan = plotting_at1_English_stan(ChannelOfInterest,:);%% Attention 1 standard - English voices

at1_odd = plotting_Ac1At1odd;
at1_odd(:,:,41:80) = plotting_Ac2At1odd;
plotting_at1_odd = mean(at1_odd,3);
plotting_at1_odd = plotting_at1_odd(ChannelOfInterest,:);% Attention 1 odd
plotting_at1_Indian_odd = mean(plotting_Ac1At1odd,3);
plotting_at1_Indian_odd = plotting_at1_Indian_odd(ChannelOfInterest,:);%% Attention 1 odd - Indian voices
plotting_at1_English_odd = mean(plotting_Ac2At1odd,3);
plotting_at1_English_odd = plotting_at1_English_odd(ChannelOfInterest,:);%% Attention 1 odd - English voices

at2_stan = plotting_Ac1At2stan;
at2_stan(:,:,41:80) = plotting_Ac2At2stan;
plotting_at2_stan = mean(at2_stan,3);
plotting_at2_stan = plotting_at2_stan(ChannelOfInterest,:);% Attention 2 standard
plotting_at2_Indian_stan = mean(plotting_Ac1At2stan,3);
plotting_at2_Indian_stan = plotting_at2_Indian_stan(ChannelOfInterest,:);%% Attention 1 standard - Indian voices
plotting_at2_English_stan = mean(plotting_Ac2At2stan,3);
plotting_at2_English_stan = plotting_at2_English_stan(ChannelOfInterest,:);%% Attention 1 standard - English voices

at2_odd = plotting_Ac1At2odd;
at2_odd(:,:,41:80) = plotting_Ac2At2odd;
plotting_at2_odd = mean(at2_odd,3);
plotting_at2_odd = plotting_at2_odd(ChannelOfInterest,:);% Attention 2 odd
plotting_at2_Indian_odd = mean(plotting_Ac1At2odd,3);
plotting_at2_Indian_odd = plotting_at2_Indian_odd(ChannelOfInterest,:);%% Attention 1 odd - Indian voices
plotting_at2_English_odd = mean(plotting_Ac2At2odd,3);
plotting_at2_English_odd = plotting_at2_English_odd(ChannelOfInterest,:);%% Attention 1 odd - English voices

% Difference wave for attention 1
at1Fig = plotting_ac1at1_DiffWaves2;
at1Fig(:,:,41:80) = plotting_ac2at1_DiffWaves2;
at1Fig = mean(at1Fig,3);
plotting_at1_diff = at1Fig(ChannelOfInterest,:); %
plotting_at1_Indian_diff = mean(plotting_ac1at1_DiffWaves2,3);
plotting_at1_Indian_diff = plotting_at1_Indian_diff(ChannelOfInterest,:);%% Attention 1 diff - Indian voices
plotting_at1_English_diff = mean(plotting_ac2at1_DiffWaves2,3);
plotting_at1_English_diff = plotting_at1_English_diff(ChannelOfInterest,:);%% Attention 1 diff - English voices

% Difference wave for attention 2
at2Fig = plotting_ac1at2_DiffWaves2;
at2Fig(:,:,41:80) = plotting_ac2at2_DiffWaves2;
at2Fig = mean(at2Fig,3);
plotting_at2_diff = at2Fig(ChannelOfInterest,:); %
plotting_at2_Indian_diff = mean(plotting_ac1at2_DiffWaves2,3);
plotting_at2_Indian_diff = plotting_at2_Indian_diff(ChannelOfInterest,:);%% Attention 1 diff - Indian voices
plotting_at2_English_diff = mean(plotting_ac2at2_DiffWaves2,3);
plotting_at2_English_diff = plotting_at2_English_diff(ChannelOfInterest,:);%% Attention 1 diff - English voices





% Indian subs difference waves - attention 1
% Indian voices
IndianSubs_plotting_ac1at1_DiffWaves2 = plotting_ac1at1_DiffWaves2(:,:,IndianSubInds);
IndianSubs_plotting_ac1at1_DiffWaves2 = mean(IndianSubs_plotting_ac1at1_DiffWaves2,3);
IndianSubs_plotting_ac1at1_DiffWaves2 = IndianSubs_plotting_ac1at1_DiffWaves2(ChannelOfInterest,:);
% English voices
IndianSubs_plotting_ac2at1_DiffWaves2 = plotting_ac2at1_DiffWaves2(:,:,IndianSubInds);
IndianSubs_plotting_ac2at1_DiffWaves2 = mean(IndianSubs_plotting_ac2at1_DiffWaves2,3);
IndianSubs_plotting_ac2at1_DiffWaves2 = IndianSubs_plotting_ac2at1_DiffWaves2(ChannelOfInterest,:);


% English subs difference waves - attention 1
% Indian voices
EnglishSubs_plotting_ac1at1_DiffWaves2 = plotting_ac1at1_DiffWaves2(:,:,EnglishSubInds);
EnglishSubs_plotting_ac1at1_DiffWaves2 = mean(EnglishSubs_plotting_ac1at1_DiffWaves2,3);
EnglishSubs_plotting_ac1at1_DiffWaves2 = EnglishSubs_plotting_ac1at1_DiffWaves2(ChannelOfInterest,:);
% English voices
EnglishSubs_plotting_ac2at1_DiffWaves2 = plotting_ac2at1_DiffWaves2(:,:,EnglishSubInds);
EnglishSubs_plotting_ac2at1_DiffWaves2 = mean(EnglishSubs_plotting_ac2at1_DiffWaves2,3);
EnglishSubs_plotting_ac2at1_DiffWaves2 = EnglishSubs_plotting_ac2at1_DiffWaves2(ChannelOfInterest,:);

% Indian subs difference waves - attention 2
% Indian voices
IndianSubs_plotting_ac1at2_DiffWaves2 = plotting_ac1at2_DiffWaves2(:,:,IndianSubInds);
IndianSubs_plotting_ac1at2_DiffWaves2 = mean(IndianSubs_plotting_ac1at2_DiffWaves2,3);
IndianSubs_plotting_ac1at2_DiffWaves2 = IndianSubs_plotting_ac1at2_DiffWaves2(ChannelOfInterest,:);
% English voices
IndianSubs_plotting_ac2at2_DiffWaves2 = plotting_ac2at2_DiffWaves2(:,:,IndianSubInds);
IndianSubs_plotting_ac2at2_DiffWaves2 = mean(IndianSubs_plotting_ac2at2_DiffWaves2,3);
IndianSubs_plotting_ac2at2_DiffWaves2 = IndianSubs_plotting_ac2at2_DiffWaves2(ChannelOfInterest,:);

% English subs difference waves - attention 2
% Indian voices
EnglishSubs_plotting_ac1at2_DiffWaves2 = plotting_ac1at2_DiffWaves2(:,:,EnglishSubInds);
EnglishSubs_plotting_ac1at2_DiffWaves2 = mean(EnglishSubs_plotting_ac1at2_DiffWaves2,3);
EnglishSubs_plotting_ac1at2_DiffWaves2 = EnglishSubs_plotting_ac1at2_DiffWaves2(ChannelOfInterest,:);
% English voices
EnglishSubs_plotting_ac2at2_DiffWaves2 = plotting_ac2at2_DiffWaves2(:,:,EnglishSubInds);
EnglishSubs_plotting_ac2at2_DiffWaves2 = mean(EnglishSubs_plotting_ac2at2_DiffWaves2,3);
EnglishSubs_plotting_ac2at2_DiffWaves2 = EnglishSubs_plotting_ac2at2_DiffWaves2(ChannelOfInterest,:);



% ----- Variable preparation for topography plots
SampleOfInterest = median(LatencyRange);
% SampleOfInterest = min(LatencyRange)+max(LatencyRange)/2; % find midpoint of measurement window
topat2Fig = at2Fig(:,SampleOfInterest); % One time point, scalp topography

topat2FigR = at2Fig(:,LatencyRange);
topat2FigR = mean(at2Fig,2); % Mean over a time window, scalp topography

% ------------------------------- Plotting ---------------------------------

% % % % Topography
if topography ==1;
    if topoplotRange == 1
        topoplot(topat2FigR, EEG.chanlocs); % Plot scalp topography! yeah!
        colormap(spring)
        set(gcf,'Units','normalized');
        annotation('textbox', [0.38,0,0.1,0.1], 'String', ['Time range ' int2str(LatencyRangeTimes) 'ms']);%
    elseif topoplotRange == 0
        topoplot(topat2Fig, EEG.chanlocs); % Plot scalp topography! yeah!
        colormap(spring)
        set(gcf,'Units','normalized');
        annotation('textbox', [0.45,0,0.1,0.1], 'String', ['Time ' int2str(EEG.times(SampleOfInterest)) 'ms']);%
         
    end
    colorbar;
    
else
    % % % % Waveforms
    
    newplot;
    hold all;
    
    if WaveOption ==1
        % Plot data
        plot(xfig,plotting_at1_stan, '-r');
        plot(xfig,plotting_at1_odd, '-b');
        plot(xfig,plotting_at1_diff, '-k');
        legend('Standard voices', 'Oddball voices', 'Oddball - standard difference', 'Location','southwest')
    elseif WaveOption ==2
        % Plot data
        plot(xfig,plotting_at2_stan, '-r');
        plot(xfig,plotting_at2_odd, '-b');
        plot(xfig,plotting_at2_diff, '-k');
        legend('Standard voices', 'Oddball voices', 'Oddball - standard difference', 'Location','southwest')
    elseif WaveOption ==3
        % Plot data
        plot(xfig,plotting_at1_Indian_diff, '-r');
        plot(xfig,plotting_at1_English_diff, '-b');
        legend('Indian voices (odd - standard difference)', 'English voices (odd - standard difference)', 'Location','southwest')
    elseif WaveOption ==4
        % Plot data
        plot(xfig,plotting_at2_Indian_diff, '-r');
        plot(xfig,plotting_at2_English_diff, '-b');
        legend('Indian voices (odd - standard difference)', 'English voices (odd - standard difference)', 'Location','southwest')
        
    elseif WaveOption ==5
        % Plot data
        plot(xfig,IndianSubs_plotting_ac1at1_DiffWaves2, '-r');
        plot(xfig,IndianSubs_plotting_ac2at1_DiffWaves2, '--r');
        plot(xfig,EnglishSubs_plotting_ac2at1_DiffWaves2, '--b');
        plot(xfig,EnglishSubs_plotting_ac1at1_DiffWaves2, '-b'); 
        legend('Indian participants own-group voices', 'Indian participants other-group voices','English participants other-group voices', 'English participants own-group voices' , 'Location','southwest')
        
    elseif WaveOption ==6
        % Plot data
        plot(xfig,IndianSubs_plotting_ac1at2_DiffWaves2, '-r');
        plot(xfig,IndianSubs_plotting_ac2at2_DiffWaves2, '--r');
        plot(xfig,EnglishSubs_plotting_ac2at2_DiffWaves2, '--b');
        plot(xfig,EnglishSubs_plotting_ac1at2_DiffWaves2, '-b');
        legend('Indian participants own-group voices', 'Indian participants other-group voices','English participants other-group voices', 'English participants own-group voices' , 'Location','southwest')
        
    elseif WaveOption ==7
        % Plot data
        
    elseif WaveOption ==8
        % Plot data
        
    end
    
    
    set(gcf,'Units','normalized');
    
    % Annotate plot
    
    plot(zeros(length(plotting_at1_stan)),plotting_at1_stan*.4, '--k');
    plot(xfig,zeros(555), '-k');
    xlabel('Time (ms)', 'FontName', 'Courier', 'FontSize', 14);
    ylabel('Mean amplitude ($u$v)', 'Interpreter', 'LaTex', 'FontName', 'Courier', 'FontSize', 14);
    % title({['Area amplitude within latencies ' int2str(min(LatencyRangeTimes)) ' and ' int2str(max(LatencyRangeTimes)) ' ms']; ['Channel ' int2str(ChannelOfInterest)]});
    annotation('textbox', [0.17,0.6,0.1,0.1], 'String', ['Electrode ' electrodes(ChannelOfInterest)]);%
%     annotation('rectangle',[0.2 0.1 .5 .85]) % rectangle highlighting measurement window
    hold off;
    
end
