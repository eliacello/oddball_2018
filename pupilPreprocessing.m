%% Elia Benhamou (01/04/2019) elia.benhamou.16@ucl.ac.uk

% Extract pupil timeseries from Eyelink - preprocessing - store data 

%% general set-up
clear all 
dbstop if error;
addpath /Users/Eliacello/Documents/MATLAB/edf-converter-master
addpath /Users/Eliacello/Dropbox/study1_oddball/analysis_pupil_122019/
addpath /Users/Eliacello/Documents/MATLAB/ICM_Eyelink_Oct2018/pupil_preprocessing_MBB/

filename='subjects_list';
sheet='V3';

[~,~,name_subj] = xlsread([filename '.xlsx'], sheet);

dirPupil='/Users/Eliacello/Documents/EBfiles_STUDY1/results_pupil_TC_Sijia/';
fs              = filesep;


dir_base = '/Users/Eliacello/Documents/EBfiles_STUDY1/';

% which sheet? 
V1 = 1; 
V2 = 0; 
V3 = 0;
% steps 
preprocessing = 1;
oddball = 0; %timeseries each condition each position
familiarity = 0; %timeseries fam unf (for 12 trials only)
dprime = 0; 
overall = 1; %overall stats of baseline and overall pupil reactivty to sounds
IC = 1; % extract discrete values from each trial (max/AUC/mean) 
entropy = 0; %timeseries for fam/unf (from start to odd) and AUC 
saveFile = 1;

timecourseLong = 1; % zscored based on pupil stats over whole timeseries
timecourseZscore =0; % zscored based on pupil stats over pre-onset

for s0= [1:38, 40:45, 48:51, 54:60, 63:67, 69:91, 93:106, 109:114, 117:length(name_subj)]
    %V1: 
    %       [1:38, 40:45, 48:51, 54:60, 63:67, 69:91, 93:106, 109:114, 117:length(name_subj)]
    %       remove 38-39 >45% rejected data ?
    %       so=39 -> change playsounds (end-1)
    %       s0 = [1:38, 40:45, 48:51, 54:60, 63:67, 69:71, 73:91, 93:106, 109:114, 117:length(name_subj)]
    %       s0= [46,47, 52:53,  61:62, 68, 92,107:108, 115,116] for subjects with delete_peripheral_fixations bug
    
    %V2: 
    %       s0=6: only 35 trials (if IC ran)
    %       s0= [1:5, 7:19, 22:27, 30:36]
    %       s0 = [20,21, 28, 29,37:38]
    %       for subjects with delete_peripheral_fixations bug
    
    %V3:
    %       s0 = [1,2,9:12] normal 
    %       s0 = 3, 13: missing trials ->adapt script and run separately
    %       s0 = 7: end - 1 (IC) 
    %       s0= [4:6, 8] for subjects with delete_peripheral_fixations bug
    
    func_dir = [dir_base name_subj{s0}];
    cd(func_dir);
    
    
    disp(name_subj{s0})
    
    %read edf file and convert to mat
    edf=[name_subj{s0} '.edf'];
    edfTim=Edf2Mat(edf);


    
    % read results file
    fid = fopen('RESULTS_FILE.txt','r');
    txt = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s');
    fclose(fid);
 
    if V1 
   
    oddball_trials=str2double(txt{1,3}(2:end));
    duration_trials=str2double(txt{1,4}(2:end));
    response=str2double(txt{1,5}(2:end));
    position=str2double(txt{1,6}(2:end));
    % songs
    songs=txt{1,1}(2:end);
    %conditions
    type=txt{1,2}(2:end);
    
    %events from edf file:
    events=edfTim.Events.Messages.info;
    
    playString='PLAY_SOUND';
    oddString='ODDBALL';
    
    end 
    
    if V2 
    
    
    index_trials=str2double(txt{1,2}(2:end));
    oddball_trials=str2double(txt{1,6}(2:end));
    duration_trials=str2double(txt{1,7}(2:end));
    response=str2double(txt{1,8}(2:end));
    position=str2double(txt{1,9}(2:end));
    % songs
    songs=txt{1,4}(2:end);
    %conditions
    type=txt{1,5}(2:end);
    
    %events from edf file:
    events=edfTim.Events.Messages.info;
    
    playString='PLAY_SOUND';
    oddString='0 ODDBALL';
    
    end 
    
    if V3
    oddball_trials=str2double(txt{1,4}(2:end));
    duration_trials=str2double(txt{1,5}(2:end));
    response=str2double(txt{1,6}(2:end));
    position=str2double(txt{1,7}(2:end));
    % songs
    songs=txt{1,1}(2:end);
    %conditions
    type=txt{1,2}(2:end);
    
    %events from edf file:
    events=edfTim.Events.Messages.info;
    
    playString='PLAY_SOUND';
    oddString='ODDBALL';
    
    end
    
    
    index_playsound = ~cellfun('isempty',strfind(events,playString)); %find any cells containing string 'TRIAL_END'
    index_oddsound = ~cellfun('isempty',strfind(events, oddString)); %find any cells containing string 'ODDBALL'
    
    time_events=edfTim.Events.Messages.time;
    playsound_onsets=time_events(index_playsound);
    oddball_onsets=time_events(index_oddsound);
    
    
    %if any oddball index missing, replace with NaN:
    sc=[];
    if length(playsound_onsets)~=length(oddball_onsets)
        for i=1:length(oddball_onsets)
            sc=[sc find(playsound_onsets(i+1)>oddball_onsets(i))];
        end
        indexMiss=length(sc==1);
        oddball_onsets=[oddball_onsets(1:indexMiss) 0 oddball_onsets(indexMiss+1:end)];
    end
    
    %fixations
    gx=edfTim.Samples.gx;
    gy=edfTim.Samples.gy;
    
    %pupil
    pupil=edfTim.Samples.pupilSize;
    
    %timestamps
    timestamps=edfTim.Samples.time;
    
    
    index_time_oddsound=[];
    
    for i=1:length(oddball_onsets)
        if mod(oddball_onsets(i),2)==0
            index_time_oddsound= [index_time_oddsound find(timestamps==oddball_onsets(i))];
        elseif mod(oddball_onsets(i),2)==1
            index_time_oddsound= [index_time_oddsound find(timestamps==oddball_onsets(i)+1)];
        end
    end
    % case timestamp odd
    for i=1:length(oddball_onsets)
        if mod(oddball_onsets(i),2)==1
            index_time_oddsound= [index_time_oddsound find(timestamps==oddball_onsets(i))];
        elseif mod(oddball_onsets(i),2)==0
            index_time_oddsound= [index_time_oddsound find(timestamps==oddball_onsets(i)+1)];
        end
    end
    
    
    index_time_playsound=[];
    
    for i=1:length(playsound_onsets)
        if mod(playsound_onsets(i),2)==0
            index_time_playsound= [index_time_playsound find(timestamps==playsound_onsets(i))];
        elseif mod(playsound_onsets(i),2)==1
            index_time_playsound= [index_time_playsound find(timestamps==playsound_onsets(i)+1)];
        end
    end
    % case timestamp odd
    for i=1:length(playsound_onsets)
        if mod(playsound_onsets(i),2)==1
            index_time_playsound= [index_time_playsound find(timestamps==playsound_onsets(i))];
        elseif mod(playsound_onsets(i),2)==0
            index_time_playsound= [index_time_playsound find(timestamps==playsound_onsets(i)+1)];
        end
    end
    
    
    % index_time_playsound=sort(index_time_playsound);
    if length(playsound_onsets)~=length(oddball_onsets)
        index_time_oddsound=[index_time_oddsound(1:indexMiss) 0 index_time_oddsound(indexMiss+1:end)];
    end

    
    %% Preprocessing
    if preprocessing
        % initialize 
        pupil_size_noz = [];
        %  1) Calculate the recorded sampling rate
        %  2) 'do_mad_outlier':     Remove outliers outside 3 MAD
        %  3) 'do_bandpass_filter': Bandpass filter
        %  4) 'do_interpolation':   Interpolate bad states
        %  5) 'do_resampling':      Resampling to 60Hz if needed
        %  6) 'do_zscore':          Z-score
        
        do_mad_outlier=1;
        do_fixationCorrection=1;
        do_bandpass_filter=1;
        do_interpolation=1;
        do_resampling=1;
        
        do_zscore=1; %no z-score for max and mean but only for timecourse analysis
        
        [pupil_size, timestampsnew, sampling_factor, pSamplesRejected] = pupil_preprocessing_canonical_reward (timestamps, pupil, gx, gy, do_mad_outlier, do_fixationCorrection, do_bandpass_filter, do_interpolation, do_resampling, do_zscore);
        
        do_zscore=0;
        [pupil_size_noz, timestampsnew, sampling_factor, pSamplesRejected] = pupil_preprocessing_canonical_reward (timestamps, pupil, gx, gy, do_mad_outlier, do_fixationCorrection, do_bandpass_filter, do_interpolation, do_resampling, do_zscore);
        
        index_time_oddsound=floor(index_time_oddsound./10); % new timestamps, only keep 50
        index_time_playsound=floor(index_time_playsound./10);
        
        TC.badSamples=sum(pSamplesRejected) * 100;
    end
    
    %% behavioral
    
    schema=find(strcmp(type,'schematic')==1);
    sem=find(strcmp(type,'semantic')==1);
    fs=find(strcmp(type,'foils_schema')==1);
    famW = find(strcmp(type,'famous_w')==1);
    fam=find(strcmp(type,'famous')==1);
    unf=find(strcmp(type,'foils')==1);
    %Hits
    semH=intersect(sem, find(response==1));
    schemaH=intersect(schema, find(response==1));
    fsH=intersect(fs, find(response==1));
    famWH=intersect(famW, find(response==1));
    %correct rejections
    famCR=intersect(fam, find(response==0));
    unfCR=intersect(unf, find(response==0));
    %false alarms
    fa_fam=intersect(fam, find(response==1));
    fa_foils=intersect(unf, find(response==1));
    %miss
    semM=intersect(sem, find(response==0));
    schemaM=intersect(schema, find(response==0));
    fsM=intersect(fs, find(response==0));
    famWM = intersect(famW, find(response==0));
    
    %familiarity
    if V3
        if s0 == 7  
            %   for fudge (n=7)
            schema (6:7) = NaN; fam (4:7)=NaN; unf(6)=NaN; 
        elseif s0 == 13
            %    for tudjoel2 (n=13)
            sem (12) = NaN; fs (12) = NaN; unf(11:12)=NaN;
            %    for boyle (n=3):
        elseif s0 == 3
            unf(10:15)=NaN;
        end
    end

    familiar_all=[schema sem fam famW];
    unfamiliar_all = [fs unf];
    
    %% compute dprimes
    
    if dprime
        h_sem=length(intersect(sem, find(response==1)));
        h_schema=length(intersect(schema, find(response==1)));
        h_fs=length(fsH);
        h_famw=length(famWH);
        fa_fam=length(intersect(fam, find(response==1)));
        fa_foils=length(intersect(unf, find(response==1)));
        h_sem= h_sem/length(sem);
        h_fs= h_fs/length(fs);
        h_schema=h_schema/length(schema);
        h_famw=h_famw/length(famW);
        fa_fam=fa_fam/length(fam);
        fa_foils=fa_foils/length(unf);
        
        perc_beh=[h_famw h_fs h_schema h_sem fa_fam fa_foils];
        
        if h_famw==1 h_famw=length(famWH)-0.5;
            h_famw=h_famw/length(famW);
        end
        if fa_fam==0 fa_fam=0.5/length(fam);end
        if fa_fam==1 fa_fam=(length(fam)-0.5)/length(fam);end
        if fa_foils==0 fa_foils=0.5/length(unf);end
        if fa_foils==1 fa_foils=(length(unf)-0.5)/length(unf);end
        if h_schema==1 h_schema=(length(schema)-0.5)/length(schema); end
        if h_schema==0 h_schema=0.5/length(schema); end
        if h_sem==1 h_sem=(length(sem)-0.5)/length(sem); end
        if h_sem==0 h_sem=0.5/length(sem); end
        if h_fs==1 h_fs=(length(fs)-0.5)/length(fs);end
        if h_fs==0 h_fs=0.5/length(fs);end
        
        dp_famw = norminv(h_famw)-norminv(fa_fam);
        dp_sem = norminv(h_sem)-norminv(fa_fam);
        dp_schema = norminv(h_schema)-norminv(fa_fam);
        dp_fs = norminv(h_fs)-norminv(fa_foils);
        dp_famw = norminv(h_famw)-norminv(fa_fam);
        dp=[dp_sem dp_schema dp_fs dp_famw];
        
    end
    %% position oddballs
    odd2=find(position==2);
    odd3=find(position==3);
    odd4=find(position==4);
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%OVERALL PUPIL METRICS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if overall
        pupilAvg = mean (pupil_size_noz); %change to pupil_zise_noz if zscore =1 during preprocessing
        pupilStd = std (pupil_size_noz) ;
        for i = 1 : length (songs) %s0=39 : change to length(songs)-1 
            baseline (i) = mean(pupil_size_noz(index_time_playsound(i)-0.5*50 ...
                :index_time_playsound(i)));
            baseline_odd (i) = mean(pupil_size_noz(index_time_oddsound(i)-0.5*50 ...
                :index_time_oddsound(i)));
        end
        pupilBasAvg = mean (baseline);
        pupilBasStd = std (baseline);
        pupilBasOddAvg = mean (baseline_odd);
        pupilBasOddStd = std (baseline_odd);
        
        
%         overall = [pupilAvg pupilStd pupilBasAvg pupilBasStd pupilBasOddAvg pupilBasOddStd];
        overall = [pupilAvg pupilStd pupilBasOddAvg pupilBasOddStd];
        TC.overall=overall;
    end
    

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ODDBALL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if oddball
        %% timecourse oddball baseline 0.5sec plotting -> epoching
        if timecourseLong %with pupil zscore with mean (X) and std (X) from whole timecourse
            
            if isempty(semH)==0 semH_timecourse_long=timecourse_long(semH, {odd2, odd3, odd4}, index_time_oddsound, pupil_size, 50); else semH_timecourse_long=[]; end
            if isempty(schemaH)==0 schemaH_timecourse_long=timecourse_long(schemaH, {odd2, odd3, odd4},index_time_oddsound, pupil_size, 50); else schemaH_timecourse_long=[]; end
            if isempty(fsH)==0 fsH_timecourse_long=timecourse_long(fsH,{odd2, odd3, odd4}, index_time_oddsound, pupil_size, 50);else fsh_timecourse_long=[];end
            if isempty(famWH)==0 famWH_timecourse_long=timecourse_long(famWH,{odd2, odd3, odd4}, index_time_oddsound, pupil_size, 50);else famWH_timecourse_long=[];end
            
            if isempty(famCR)==0 famCR_timecourse_long=timecourse_long(famCR,{odd2, odd3, odd4}, index_time_oddsound, pupil_size, 50);else famCR_timecourse_long=[]; end
            if isempty(unfCR)==0 unfCR_timecourse_long=timecourse_long(unfCR, {odd2, odd3, odd4},index_time_oddsound, pupil_size, 50);else unfCR_timecourse_long=[];end %s0=39: 1:end-1 playsound
            
            if isempty(fa_fam)==0 fa_fam_timecourse_long=timecourse_long(fa_fam,{odd2, odd3, odd4}, index_time_oddsound, pupil_size, 50);else fa_fam_timecourse_long=[];end
            if isempty(fa_foils)==0 fa_foils_timecourse_long=timecourse_long(fa_foils,{odd2, odd3, odd4}, index_time_oddsound, pupil_size, 50); else fa_foils_timecourse_long=[];end
            
            if isempty(semM)==0 semM_timecourse_long=timecourse_long(semM, {odd2, odd3, odd4},index_time_oddsound, pupil_size, 50);else semM_timecourse_long=[]; end
            if isempty(schemaM)==0 schemaM_timecourse_long=timecourse_long(schemaM, {odd2, odd3, odd4},index_time_oddsound, pupil_size, 50);else schemaM_timecourse_long=[]; end
            if isempty(fsM)==0 fsM_timecourse_long=timecourse_long(fsM,{odd2, odd3, odd4}, index_time_oddsound, pupil_size, 50);else fsM_timecourse_long=[]; end
            if isempty(famWM)==0 famWM_timecourse_long=timecourse_long(famWM,{odd2, odd3, odd4}, index_time_oddsound, pupil_size, 50);else famWM_timecourse_long=[]; end
            
        end
        
        if timecourseZscore %with pupil zscore with mean (X) and std (X) from baseline
            if isempty(semH)==0
                semH_timecourse_zscore=timecourse_zscore(semH, {odd2, odd3, odd4}, index_time_oddsound, pupil_size, pupilBasOddAvg, pupilBasOddStd, 50);
            else semH_timecourse_zscore=[];
            end
            if isempty(schemaH)==0 schemaH_timecourse_zscore=timecourse_zscore(schemaH, {odd2, odd3, odd4},index_time_oddsound, pupil_size,pupilBasOddAvg, pupilBasOddStd, 50); else schemaH_timecourse_zscore=[]; end
            if isempty(fsH)==0 fsH_timecourse_zscore=timecourse_zscore(fsH,{odd2, odd3, odd4}, index_time_oddsound, pupil_size,pupilBasOddAvg, pupilBasOddStd, 50);else fsh_timecourse_zscore=[];end
            if isempty(famWH)==0 famWH_timecourse_zscore=timecourse_zscore(famWH,{odd2, odd3, odd4}, index_time_oddsound, pupil_size,pupilBasOddAvg, pupilBasOddStd, 50);else famWH_timecourse_zscore=[];end
            
            if isempty(famCR)==0 famCR_timecourse_zscore=timecourse_zscore(famCR,{odd2, odd3, odd4}, index_time_oddsound, pupil_size,pupilBasOddAvg, pupilBasOddStd, 50);else famCR_timecourse_zscore=[]; end
            if isempty(unfCR)==0 unfCR_timecourse_zscore=timecourse_zscore(unfCR(1:end-1) , {odd2, odd3, odd4},index_time_oddsound, pupil_size,pupilBasOddAvg, pupilBasOddStd, 50);else unfCR_timecourse_zscore=[];end %s0=39: 1:end-1 playsound
            
            if isempty(fa_fam)==0 fa_fam_timecourse_zscore=timecourse_zscore(fa_fam,{odd2, odd3, odd4}, index_time_oddsound, pupil_size,pupilBasOddAvg, pupilBasOddStd, 50);else fa_fam_timecourse_zscore=[];end
            if isempty(fa_foils)==0 fa_foils_timecourse_zscore=timecourse_zscore(fa_foils,{odd2, odd3, odd4}, index_time_oddsound, pupil_size,pupilBasOddAvg, pupilBasOddStd, 50); else fa_foils_timecourse_zscore=[];end
            
            if isempty(semM)==0 semM_timecourse_zscore=timecourse_zscore(semM, {odd2, odd3, odd4},index_time_oddsound, pupil_size,pupilBasOddAvg, pupilBasOddStd, 50);else semM_timecourse_zscore=[]; end
            if isempty(schemaM)==0 schemaM_timecourse_zscore=timecourse_zscore(schemaM, {odd2, odd3, odd4},index_time_oddsound, pupil_size,pupilBasOddAvg, pupilBasOddStd, 50);else schemaM_timecourse_zscore=[]; end
            if isempty(fsM)==0 fsM_timecourse_zscore=timecourse_zscore(fsM,{odd2, odd3, odd4}, index_time_oddsound, pupil_size,pupilBasOddAvg, pupilBasOddStd, 50);else fsM_timecourse_zscore=[]; end
            if isempty(famWM)==0 famWM_timecourse_zscore=timecourse_zscore(famWM,{odd2, odd3, odd4}, index_time_oddsound, pupil_size,pupilBasOddAvg, pupilBasOddStd, 50);else famWM_timecourse_zscore=[]; end
            
        end
        
        %% sample conditions schematic-semantic-famous-foils-fs w/o position effect
        
        % schematic oddballs in foils (fs)
        if timecourseLong
            fs_sample=fsH_timecourse_long.position;
        end
        if timecourseZscore
            fs_sample=fsH_timecourse_zscore.position;
        end
        S=size(fs_sample,1);
        empties=cellfun('isempty', fs_sample);
        fs_sample(empties)={NaN([226,1])};
        fs_sample=cell2mat(fs_sample);
        if size(fs_sample,2)==3
            fs_pos1=reshape(fs_sample(:,1), [226,S]); fs_pos1M=nanmean(fs_pos1,2);
            fs_pos2=reshape(fs_sample(:,2), [226,S]); fs_pos2M=nanmean(fs_pos2,2);
            fs_pos3=reshape(fs_sample(:,3), [226,S]); fs_pos3M=nanmean(fs_pos3,2);
            
        elseif size(fs_sample,2)==2
            fs_pos1=reshape(fs_sample(:,1), [226,S]); fs_pos1M=nanmean(fs_pos1,2);
            fs_pos2=reshape(fs_sample(:,2), [226,S]); fs_pos2M=nanmean(fs_pos2,2);
            fs_pos3M=[];
        elseif size(fs_sample,2)==1
            fs_pos1=reshape(fs_sample(:,1), [226,S]); fs_pos1M=nanmean(fs_pos1,2);
            fs_pos2M=[];
            fs_pos3M=[];
        end
        fsA=[fs_pos1M fs_pos2M fs_pos3M];
        
        %
        TC.fs=fsA;
        
        %semantic oddballs in famous songs (sem)
        if timecourseLong
            if isempty (semH_timecourse_long) ==0
                sem_sample=semH_timecourse_long.position;
            else
                TC.sem = [];
                
            end
        end
        if timecourseZscore
            if isempty (semH_timecourse_zscore) ==0
                sem_sample=semH_timecourse_zscore.position;
                
                S=size(sem_sample,1);
                empties=cellfun('isempty', sem_sample);
                sem_sample(empties)={NaN([226,1])};
                sem_sample=cell2mat(sem_sample);
                if size(sem_sample,2)==3
                    sem_pos1=reshape(sem_sample(:,1), [226,S]); sem_pos1M=nanmean(sem_pos1,2);
                    sem_pos2=reshape(sem_sample(:,2), [226,S]); sem_pos2M=nanmean(sem_pos2,2);
                    sem_pos3=reshape(sem_sample(:,3), [226,S]); sem_pos3M=nanmean(sem_pos3,2);
                elseif size(sem_sample,2)==2
                    sem_pos1=reshape(sem_sample(:,1), [226,S]); sem_pos1M=nanmean(sem_pos1,2);
                    sem_pos2=reshape(sem_sample(:,2), [226,S]); sem_pos2M=nanmean(sem_pos2,2);
                    sem_pos3M=[];
                elseif size(sem_sample,2)==1
                    sem_pos1=reshape(sem_sample(:,1), [226,S]); sem_pos1M=nanmean(sem_pos1,2);
                    sem_pos2M=[];
                    sem_pos3M=[];
                end
                semA=[sem_pos1M sem_pos2M sem_pos3M];

                TC.sem=semA;
            else
                TC.sem = [];
            end
        end
        
        % schematic oddballs in famous songs (schema)
        if timecourseLong
            if isempty (schemaH) ==0
                sem_sample=semH_timecourse_long.position;
            else
                TC.schema = [];
    
            end
        end
        if timecourseZscore
            if isempty (schemaH) ==0
                schema_sample=schemaH_timecourse_zscore.position;
                
                S=size(schema_sample,1);
                empties=cellfun('isempty', schema_sample);
                schema_sample(empties)={NaN([226,1])};
                schema_sample=cell2mat(schema_sample);
                if size(schema_sample,2)==3
                    schema_pos1=reshape(schema_sample(:,1), [226,S]); schema_pos1M=nanmean(schema_pos1,2);
                    schema_pos2=reshape(schema_sample(:,2), [226,S]); schema_pos2M=nanmean(schema_pos2,2);
                    schema_pos3=reshape(schema_sample(:,3), [226,S]); schema_pos3M=nanmean(schema_pos3,2);
                elseif size(schema_sample,2)==2
                    schema_pos1=reshape(schema_sample(:,1), [226,S]); schema_pos1M=nanmean(schema_pos1,2);
                    schema_pos2=reshape(schema_sample(:,2), [226,S]); schema_pos2M=nanmean(schema_pos2,2);
                    schema_pos3M=[];
                elseif size(schema_sample,2)==1
                    schema_pos1=reshape(schema_sample(:,1), [226,S]); schema_pos1M=nanmean(schema_pos1,2);
                    schema_pos2M=[];
                    schema_pos3M=[];
                end
                
                schemaA=[schema_pos1M schema_pos2M schema_pos3M];
                
                TC.schema=schemaA;
            else
                TC.sem = [];
                
                
            end
        end

        % famous songs with white-noise (famW)
        if timecourseLong
            if isempty (famWH) ==0
                famW_sample=famWH_timecourse_long.position;
            else
                TC.famW = [];
                
            end
        end
        if timecourseZscore
            if isempty (famWH) ==0
                famW_sample=famWH_timecourse_zscore.position;
                S=size(famW_sample,1);
                empties=cellfun('isempty', famW_sample);
                famW_sample(empties)={NaN([226,1])};
                famW_sample=cell2mat(famW_sample);
                if size(famW_sample,2)==3
                    famW_pos1=reshape(famW_sample(:,1), [226,S]); famW_pos1M=nanmean(famW_pos1,2);
                    famW_pos2=reshape(famW_sample(:,2), [226,S]); famW_pos2M=nanmean(famW_pos2,2);
                    famW_pos3=reshape(famW_sample(:,3), [226,S]); famW_pos3M=nanmean(famW_pos3,2);
                elseif size(famW_sample,2)==2
                    famW_pos1=reshape(famW_sample(:,1), [226,S]); famW_pos1M=nanmean(famW_pos1,2);
                    famW_pos2=reshape(famW_sample(:,2), [226,S]); famW_pos2M=nanmean(famW_pos2,2);
                    famW_pos3M=[];
                elseif size(famW_sample,2)==1
                    famW_pos1=reshape(famW_sample(:,1), [226,S]); famW_pos1M=nanmean(famW_pos1,2);
                    famW_pos2M=[];
                    famW_pos3M=[];
                end
                
                famWA=[famW_pos1M famW_pos2M famW_pos3M];
                
                TC.famW=famWA;
            else
                TC.famW = [];
                
                
            end
        end
        
        
        %famous songs with no oddball (famous)
        if timecourseLong
            if isempty (famCR) ==0
                famous_sample=famCR_timecourse_long.position;
            else
                TC.famous = [];
                
            end
        end
        if timecourseZscore
            if isempty (famCR) ==0
                famous_sample=famCR_timecourse_zscore.position;
                
                
                S=size(famous_sample,1);
                empties=cellfun('isempty', famous_sample);
                famous_sample(empties)={NaN([226,1])};
                famous_sample=cell2mat(famous_sample);
                if size(famous_sample,2)==3
                    famous_pos1=reshape(famous_sample(:,1), [226,S]); famous_pos1M=nanmean(famous_pos1,2);
                    famous_pos2=reshape(famous_sample(:,2), [226,S]); famous_pos2M=nanmean(famous_pos2,2);
                    famous_pos3=reshape(famous_sample(:,3), [226,S]); famous_pos3M=nanmean(famous_pos3,2);
                elseif size(famous_sample,2)==2
                    famous_pos1=reshape(famous_sample(:,1), [226,S]); famous_pos1M=nanmean(famous_pos1,2);
                    famous_pos2=reshape(famous_sample(:,2), [226,S]); famous_pos2M=nanmean(famous_pos2,2);
                    famous_pos3M=[];
                elseif size(famous_sample,2)==1
                    famous_pos1=reshape(famous_sample(:,1), [226,S]); famous_pos1M=nanmean(famous_pos1,2);
                    famous_pos2M=[];
                    famous_pos3M=[];
                end
                famousA=[famous_pos1M famous_pos2M famous_pos3M];
                
                
                TC.famous=famousA;
            else
                TC.famous = [];
                
            end
        end
        
        % foils with no oddballs 
        
        if timecourseLong
            foils_sample = unfCR_timecourse_long.position;
        end
        if timecourseZscore
            foils_sample = unfCR_timecourse_zscore.position;
        end
        
        S=size(foils_sample,1);
        empties=cellfun('isempty', foils_sample);
        foils_sample(empties)={NaN([226,1])};
        foils_sample=cell2mat(foils_sample);
        if size(foils_sample,2)==3
            foils_pos1=reshape(foils_sample(:,1), [226,S]); foils_pos1M=nanmean(foils_pos1,2);
            foils_pos2=reshape(foils_sample(:,2), [226,S]); foils_pos2M=nanmean(foils_pos2,2);
            foils_pos3=reshape(foils_sample(:,3), [226,S]); foils_pos3M=nanmean(foils_pos3,2);
        elseif size(foils_sample,2)==2
            foils_pos1=reshape(foils_sample(:,1), [226,S]); foils_pos1M=nanmean(foils_pos1,2);
            foils_pos2=reshape(foils_sample(:,2), [226,S]); foils_pos2M=nanmean(foils_pos2,2);
            foils_pos3M=[];
        elseif size(foils_sample,2)==1
            foils_pos1=reshape(foils_sample(:,1), [226,S]); foils_pos1M=nanmean(foils_pos1,2);
            foils_pos2M=[];
            foils_pos3M=[];
        end
        
        foilsA=[foils_pos1M foils_pos2M foils_pos3M];
        
        TC.foils=foilsA;
        
        % false alarms in famous songs (fa_fam) 
        if timecourseLong
            if isempty (fa_fam_timecourse_long) ==0
                fa_fam_sample=fa_fam_timecourse_long.position;
            else
                TC.fa_fam = [];
                
            end
        end
        if timecourseZscore
            if isempty (fa_fam_timecourse_zscore) ==0
                fa_fam_sample=fa_fam_timecourse_zscore.position;
                
                
                S=size(fa_fam_sample,1);
                empties=cellfun('isempty', fa_fam_sample);
                fa_fam_sample(empties)={NaN([226,1])};
                fa_fam_sample=cell2mat(fa_fam_sample);
                if size(fa_fam_sample,2)==3
                    fa_fam_pos1=reshape(fa_fam_sample(:,1), [226,S]); fa_fam_pos1M=nanmean(fa_fam_pos1,2);
                    fa_fam_pos2=reshape(fa_fam_sample(:,2), [226,S]); fa_fam_pos2M=nanmean(fa_fam_pos2,2);
                    fa_fam_pos3=reshape(fa_fam_sample(:,3), [226,S]); fa_fam_pos3M=nanmean(fa_fam_pos3,2);
                elseif size(fa_fam_sample,2)==3
                    fa_fam_pos1=reshape(fa_fam_sample(:,1), [226,S]); fa_fam_pos1M=nanmean(fa_fam_pos1,2);
                    fa_fam_pos2=reshape(fa_fam_sample(:,2), [226,S]); fa_fam_pos2M=nanmean(fa_fam_pos2,2);
                    fa_fam_pos3M=[];
                elseif size(fa_fam_sample,2)==1
                    fa_fam_pos1=reshape(fa_fam_sample(:,1), [226,S]); fa_fam_pos1M=nanmean(fa_fam_pos1,2);
                    fa_fam_pos2M=[];
                    fa_fam_pos3M=[];
                    
                end
                fa_famA=[fa_fam_pos1M fa_fam_pos2M fa_fam_pos3M];
                
                TC.fa_fam=fa_famA;
            else
                TC.fa_fam = [];
                
            end
        end
        
        
        % false-alarms in foils (fa_foils)
        if timecourseLong
            if isempty (fa_foils_timecourse_long) ==0
                fa_foils_sample=fa_foils_timecourse_long.position;
            else
                TC.fa_foils = [];
                
            end
        end
        if timecourseZscore
            if isempty (fa_foils_timecourse_zscore) ==0
                fa_foils_sample=fa_foils_timecourse_zscore.position;
                
                
                S=size(fa_foils_sample,1);
                empties=cellfun('isempty', fa_foils_sample);
                fa_foils_sample(empties)={NaN([226,1])};
                fa_foils_sample=cell2mat(fa_foils_sample);
                if size(fa_foils_sample,2)==3
                    fa_foils_pos1=reshape(fa_foils_sample(:,1), [226,S]); fa_foils_pos1M=nanmean(fa_foils_pos1,2);
                    fa_foils_pos2=reshape(fa_foils_sample(:,2), [226,S]); fa_foils_pos2M=nanmean(fa_foils_pos2,2);
                    fa_foils_pos3=reshape(fa_foils_sample(:,3), [226,S]); fa_foils_pos3M=nanmean(fa_foils_pos3,2);
                elseif size(fa_foils_sample,2)==2
                    fa_foils_pos1=reshape(fa_foils_sample(:,1), [226,S]); fa_foils_pos1M=nanmean(fa_foils_pos1,2);
                    fa_foils_pos2=reshape(fa_foils_sample(:,2), [226,S]); fa_foils_pos2M=nanmean(fa_foils_pos2,2);
                    fa_foils_pos3M=[];
                elseif size(fa_foils_sample,2)==1
                    fa_foils_pos1=reshape(fa_foils_sample(:,1), [226,S]); fa_foils_pos1M=nanmean(fa_foils_pos1,2);
                    fa_foils_pos2M=[];
                    fa_foils_pos3M=[];
                    
                    
                end
                fa_foilsA=[fa_foils_pos1M fa_foils_pos2M fa_foils_pos3M];
                
                
                TC.fa_foils=fa_foilsA;
            else
                TC.fa_foils = [];
                
            end
        end
        
        
        
        % miss semantic oddballs
        if timecourseLong
            if isempty (semM_timecourse_long) ==0
                semM_sample=semM_timecourse_long.position;
            else
                TC.semM = [];
                
            end
        end
        if timecourseZscore
            if isempty (semM_timecourse_zscore) ==0
                semM_sample=semM_timecourse_zscore.position;
                
                
                S=size(semM_sample,1);
                empties=cellfun('isempty', semM_sample);
                semM_sample(empties)={NaN([226,1])};
                semM_sample=cell2mat(semM_sample);
                if size(semM_sample,2)==3
                    semM_pos1=reshape(semM_sample(:,1), [226,S]); semM_pos1M=nanmean(semM_pos1,2);
                    semM_pos2=reshape(semM_sample(:,2), [226,S]); semM_pos2M=nanmean(semM_pos2,2);
                    semM_pos3=reshape(semM_sample(:,3), [226,S]); semM_pos3M=nanmean(semM_pos3,2);
                elseif size(semM_sample,2)==2
                    semM_pos1=reshape(semM_sample(:,1), [226,S]); semM_pos1M=nanmean(semM_pos1,2);
                    semM_pos2=reshape(semM_sample(:,2), [226,S]); semM_pos2M=nanmean(semM_pos2,2);
                    semM_pos3M=[];
                elseif size(semM_sample,2)==1
                    semM_pos1=reshape(semM_sample(:,1), [226,S]); semM_pos1M=nanmean(semM_pos1,2);
                    semM_pos2M=[];
                    semM_pos3M=[];
                end
                semMA=[semM_pos1M semM_pos2M semM_pos3M];
                
                TC.semM=semMA;
            else
                TC.semM = [];
                
            end
        end
        
        
        % miss schematic oddballs
        if timecourseLong
            if isempty (schemaM_timecourse_long) ==0
                schemaM_sample=schemaM_timecourse_long.position;
            else
                TC.schemaM = [];
                
            end
        end
        if timecourseZscore
            if isempty (schemaM_timecourse_zscore) ==0
                schemaM_sample=schemaM_timecourse_zscore.position;
                
                
                S=size(schemaM_sample,1);
                empties=cellfun('isempty', schemaM_sample);
                schemaM_sample(empties)={NaN([226,1])};
                schemaM_sample=cell2mat(schemaM_sample);
                if size(schemaM_sample,2)==3
                    schemaM_pos1=reshape(schemaM_sample(:,1), [226,S]); schemaM_pos1M=nanmean(schemaM_pos1,2);
                    schemaM_pos2=reshape(schemaM_sample(:,2), [226,S]); schemaM_pos2M=nanmean(schemaM_pos2,2);
                    schemaM_pos3=reshape(schemaM_sample(:,3), [226,S]); schemaM_pos3M=nanmean(schemaM_pos3,2);
                elseif size(schemaM_sample,2)==2
                    schemaM_pos1=reshape(schemaM_sample(:,1), [226,S]); schemaM_pos1M=nanmean(schemaM_pos1,2);
                    schemaM_pos2=reshape(schemaM_sample(:,2), [226,S]); schemaM_pos2M=nanmean(schemaM_pos2,2);
                    schemaM_pos3M=[];
                elseif size(schemaM_sample,2)==1
                    schemaM_pos1=reshape(schemaM_sample(:,1), [226,S]); schemaM_pos1M=nanmean(schemaM_pos1,2);
                    schemaM_pos2M=[];
                    schemaM_pos3M=[];
                end
                schemaMA=[schemaM_pos1M schemaM_pos2M schemaM_pos3M];
                
                TC.schemaM=schemaMA;
            else
                TC.schemaM = [];
                
            end
        end
        
        
        %miss schematic oddballs in foils
        if timecourseLong
            if isempty (fsM_timecourse_long) ==0
                fsM_sample=fsM_timecourse_long.position;
            else
                TC.fsM = [];
                
            end
        end
        if timecourseZscore
            if isempty (fsM_timecourse_zscore) ==0
                fsM_sample=fsM_timecourse_zscore.position;
                
                
                S=size(fsM_sample,1);
                empties=cellfun('isempty', fsM_sample);
                fsM_sample(empties)={NaN([226,1])};
                fsM_sample=cell2mat(fsM_sample);
                if size(fsM_sample,2)==3
                    fsM_pos1=reshape(fsM_sample(:,1), [226,S]); fsM_pos1M=nanmean(fsM_pos1,2);
                    fsM_pos2=reshape(fsM_sample(:,2), [226,S]); fsM_pos2M=nanmean(fsM_pos2,2);
                    fsM_pos3=reshape(fsM_sample(:,3), [226,S]); fsM_pos3M=nanmean(fsM_pos3,2);
                elseif size(fsM_sample,2)==2
                    fsM_pos1=reshape(fsM_sample(:,1), [226,S]); fsM_pos1M=nanmean(fsM_pos1,2);
                    fsM_pos2=reshape(fsM_sample(:,2), [226,S]); fsM_pos2M=nanmean(fsM_pos2,2);
                    fsM_pos3M=[];
                elseif size(fsM_sample,2)==1
                    fsM_pos1=reshape(fsM_sample(:,1), [226,S]); fsM_pos1M=nanmean(fsM_pos1,2);
                    fsM_pos2M=[];
                    fsM_pos3M=[];
                end
                fsMA=[fsM_pos1M fsM_pos2M fsM_pos3M];
                
                TC.fsM=fsMA;
            else
                TC.fsM = [];
                
            end
        end
        
        
        %miss white noise in famous songs
        
        if timecourseLong
            if isempty (famWM_timecourse_long) ==0
                famWM_sample=famWM_timecourse_long.position;
            else
                TC.famWM = [];
                
            end
        end
        if timecourseZscore
            if isempty (famWM_timecourse_zscore) ==0
                famWM_sample=famWM_timecourse_zscore.position;
                
                
                
                S=size(famWM_sample,1);
                empties=cellfun('isempty', famWM_sample);
                famWM_sample(empties)={NaN([226,1])};
                famWM_sample=cell2mat(famWM_sample);
                if size(famWM_sample,2)==3
                    famWM_pos1=reshape(famWM_sample(:,1), [226,S]); famWM_pos1M=nanmean(famWM_pos1,2);
                    famWM_pos2=reshape(famWM_sample(:,2), [226,S]); famWM_pos2M=nanmean(famWM_pos2,2);
                    famWM_pos3=reshape(famWM_sample(:,3), [226,S]); famWM_pos3M=nanmean(famWM_pos3,2);
                elseif size(famWM_sample,2)==2
                    famWM_pos1=reshape(famWM_sample(:,1), [226,S]); famWM_pos1M=nanmean(famWM_pos1,2);
                    famWM_pos2=reshape(famWM_sample(:,2), [226,S]); famWM_pos2M=nanmean(famWM_pos2,2);
                    famWM_pos3M=[];
                elseif size(famWM_sample,2)==1
                    famWM_pos1=reshape(famWM_sample(:,1), [226,S]); famWM_pos1M=nanmean(famWM_pos1,2);
                    famWM_pos2M=[];
                    famWM_pos3M=[];
                end
                famWMA=[famWM_pos1M famWM_pos2M famWM_pos3M];
                
                TC.famWM=famWMA;
            else
                TC.famWM = [];
                
            end
        end
        
    end
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%INFORMATION_CONTENT ODDBALL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if IC
        
        % use timecourseODDIC and timecourseAUC if you want to baseline-correct pupil already
        % zscored (at subject-level)
        % use timecourseODDIC_zscore and timecourseAUC_zscore if you want to z-score using mu and
        % sigma from baseline_odd and then baseline-correct 
        % attention!!!! initialize IC_ODD
        IC_ODD = []; 
        
         IC_timecourse = timecourseODDIC_zscore( index_time_oddsound, pupil_size, pupilBasOddAvg,pupilBasOddStd, 50); %add pupilBasOddAvg,pupilBasOddStd if using timecourseODDIC_zscore
         IC_timecourse_noz = timecourseODDIC_zscore ( index_time_oddsound, pupil_size_noz, pupilBasOddAvg,pupilBasOddStd, 50);
        
        AUC_timecourse = timecourseAUC_zscore (index_time_playsound, index_time_oddsound, pupil_size, pupilBasAvg,pupilBasStd,  50); %s0=39(V1) and s0=7(V3): 1:end-1 playsound
         AUC_timecourse_noz = timecourseAUC_zscore (index_time_playsound, index_time_oddsound, pupil_size_noz,pupilBasAvg,pupilBasStd, 50); %s0=39(V1) and s0=7 (V3): 1:end-1 playsound
        
        IC_sample = cell2mat(IC_timecourse);
        IC_sample_noz = cell2mat(IC_timecourse_noz);
        
        for i=1:length(songs) %s0=39 (V1), s0=6 (V2) :-1, s0=7 (V3)
            if strcmp(type(i),'schematic')==1 && response (i)==1
                accuracy{i,1} = 1;
            elseif strcmp(type(i),'semantic')==1 && response (i)==1
                accuracy{i,1} = 1;
            elseif strcmp(type(i),'foils_schema')==1 && response (i)==1
                accuracy{i,1} = 1;
            elseif strcmp(type(i),'famous_w')==1 && response (i)==1
                accuracy{i,1} = 1;
            elseif strcmp(type(i),'famous')==1 && response (i)==0
                accuracy{i,1} = 1;
            elseif strcmp(type(i),'foils')==1 && response (i)==0
                
                accuracy{i,1} = 1;
            else
                accuracy{i,1} = 0;
            end
            
            [a0,atime]=max(IC_sample(25:125,i));
            
            IC_ODD{i,1} = name_subj{s0};
            IC_ODD{i,2} = songs{i,1};
            IC_ODD{i,3} = max (IC_sample(25:125,i)); %only 2 first sec
            IC_ODD{i,4} = max (IC_sample_noz(25:125,i)); %only 2 first sec
            IC_ODD{i,5} = accuracy{i,1};
            IC_ODD{i,6} = type{i,1};
            IC_ODD{i,7} = trapz(AUC_timecourse{1,i});
            IC_ODD{i,8} = trapz(AUC_timecourse_noz{1,i});
            IC_ODD{i,9} = mean(AUC_timecourse{1,i});
            IC_ODD{i,10} = mean(AUC_timecourse_noz{1,i});
            IC_ODD{i,11} = atime/50; 
            IC_ODD{i,12} = IC_sample_noz(:, i); 
            IC_ODD{i,13} = IC_sample(:, i); 
            IC_ODD{i,14} = trapz(IC_sample(25:125,i)); 
            IC_ODD{i,15} = trapz(IC_sample_noz(25:125,i)); 
            IC_ODD{i,16} = AUC_timecourse_noz;
            IC_ODD{i,17} = AUC_timecourse;

        end
        
        TC.ICODD=IC_ODD;
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%FAMILIARITY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % timecourse familiarity (play_sound + maximum duration possible) baseline 0.5sec plotting
    % only for fam and unf trials (12 each per subject)
    if familiarity
        
        
        dur = min (min(duration_trials (fam)), min(duration_trials (unf))); 
        
        fam_timecourse = timecourse_fam_zscore(fam,  index_time_playsound(1:end-1), pupil_size, dur, pupilBasAvg, pupilBasStd, 50);%s0=39 (V1) and s0=7(V3): 1:end-1 playsound
        unf_timecourse = timecourse_fam_zscore(unf(1:end-1),  index_time_playsound(1:end-1), pupil_size, dur, pupilBasAvg, pupilBasStd, 50);%s0=39 and s0=7(V3): 1:end-1 playsound and unf
        
        
        % sample conditions familiarity melodies (+3sec)
        
        
        fam_sample = cell2mat(fam_timecourse);
        unf_sample = cell2mat(unf_timecourse);
        
        TC.fam = fam_sample;
        TC.unf = unf_sample;
        
    end
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%ENTROPY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if entropy
        % get entropy scores from excel file
        %     filename='IC_SONGS_ODD';
        %     sheet1='IC_SONG';
        %     sheet2='IC_ODD';
        
        %
        %sheet 1: familiarity
        %extract timecourse from start song to oddball and apply AUC
        % for all
        
        % attention!!!!! initialize AUC 
        AUC = []; 
        
        AUC_timecourse = timecourseAUC_zscore (index_time_playsound(1:end-1), index_time_oddsound, pupil_size,pupilBasAvg,pupilBasStd, 50); 
        
        
        for i=1:length(songs)-1
            
            AUC{i,1} = trapz(AUC_timecourse{1,i});
            AUC{i,2} = songs{i,1};
            
        end
        
        
        famSO_timecourse = timecourseSO_zscore (familiar_all, index_time_playsound(1:end-1), index_time_oddsound, pupil_size,pupilBasAvg,pupilBasStd, 50); 
        unfSO_timecourse = timecourseSO_zscore (unfamiliar_all(1:end-1), index_time_playsound(1:end-1), index_time_oddsound, pupil_size,pupilBasAvg,pupilBasStd, 50); 
       
        
        TC.famSO = famSO_timecourse;
        TC.unfSO = unfSO_timecourse;
        
        
        TC.AUC=AUC;
    end
    
    
    %% save variables
    
    if saveFile
        
        matdatafile = fullfile(dirPupil,name_subj{s0});
        save(matdatafile,'TC')
    end
    
    
end


