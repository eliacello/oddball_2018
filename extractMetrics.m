%% extract behavioural and pupil metrics from trials_allsubjects
% Jan/Feb2020

addpath /Users/Eliacello/Dropbox/study1_oddball/chapter3_PW/
filename='trials_allSubjects_famAnalysis'; %or 'trials_allSubjects_famAnalysis'; 'trials_allSubjects'
sheet='Sheet1';

[~,~,file] = xlsread([filename '.xlsx'], sheet);
name = file (2:end, 1);
type = file (2:end, 9);
[ID, who] = findgroups (name);
group = cell2mat (file (2:end, 3));
entropysong = file (2:end, 14); 
entropysong = cell2mat(entropysong);
ICODD = file(2:end, 13); 
ICODD = cell2mat(ICODD); 


songs = file (2:end, 4); 
[song, nameSong] = findgroups (songs); 
[condition, t] = findgroups (type);


acc = file (2:end, 7);
accuracy = cell2mat(acc);

max1 = cell2mat(file(2:end, 5)); 
max2 = cell2mat(file(2:end, 6));

trap1 = cell2mat(file(2:end, 9)); 
trap2 = cell2mat(file(2:end, 10));

meanSO1 = cell2mat(file(2:end, 11)); 
meanSO2 = cell2mat(file(2:end, 12));

UNFRATING = cell2mat (file(2:end, 26));
UNFSURVEY = cell2mat (file(2:end, 24));

FAMRATING = cell2mat (file(2:end, 25));
FAMSURVEY = cell2mat (file(2:end, 23));

fam = find (condition == 1);
fam_w = find (condition ==2);
foils = find (condition ==3);
fs = find (condition ==4);
schema = find (condition ==7);
sem = find (condition ==8);

fam_all = [fam; fam_w; schema; sem];
fam_schema = [fam; schema]; % for chapter3
unf_all = [foils; fs];

%only keep trials rated as familiar: 
fam_rating = 1; 

% generate score for fam: 
score = zeros(length(ID), 1); 
for i = 1:length (ID)
    if FAMRATING(i) == 1 && FAMSURVEY(i) == 1
    score(i) = 1; 
    elseif FAMRATING(i) == 0 && FAMSURVEY(i) == 0
        score(i) = 1; 
    else 
        score (i) = 0; 
    end 
end 

%% only keep maxzscore <2*SD 

sigma = nanstd (max1); 
noisyTrials = find(max1>2*sigma); 

%% average each variable for each subject WITHOUT NOISY TRIALS 

% 1) max (IC_sample(25:125,i)) for each condition (+4sec after deviant)
% 2) max (IC_sample_noz(25:125,i)) for each condition (+4sec after deviant)
% 3) trapz(AUC_timecourse{1,i}) for fam/unf (start to deviant)
% 4) trapz(AUC_timecourse_noz{1,i}) for fam/unf (start to deviant)
% 5) mean(AUC_timecourse{1,i}) for fam/unf (start to deviant)
% 6) mean(AUC_timecourse_noz{1,i}) for fam/unf (start to deviant)
var1.fam = []; var2.fam = []; var3.fam = []; var4.fam = []; var5.fam = []; var6.fam = [];
 var1.famw = []; var1.foils = []; var1.fs = []; var1.sem = []; var1.schema=[];
 var2.famw = []; var2.foils = []; var2.fs = []; var2.sem = []; var2.schema=[];
var3.unf = []; var4.unf = []; var5.unf = []; var6.unf = [];
groupID=[]; 

% without noisy trials: 
max1 (noisyTrials) = [];  
group (noisyTrials) = []; 
ID (noisyTrials) = []; 
condition  (noisyTrials) = []; 

fam = find (condition == 1);
fam_w = find (condition ==2);
foils = find (condition ==3);
fs = find (condition ==4);
schema = find (condition ==7);
sem = find (condition ==8);

fam_all = [fam; fam_w; schema; sem];
fam_schema = [fam; schema]; % for chapter3
unf_all = [foils; fs];
for i=1:length (unique(ID))
        j = find(ID==i) ;
    
    groupID = [groupID; group(j(1))];
    
    var1.fam = [var1.fam; nanmean(max1 (intersect(fam, find(ID==i))))];
    var1.famw = [var1.famw; nanmean(max1 (intersect(fam_w, find(ID==i))))];
    var1.foils = [var1.foils; nanmean(max1 (intersect(foils, find(ID==i))))];
    var1.fs = [var1.fs; nanmean(max1 (intersect(fs, find(ID==i))))];
    var1.sem = [var1.sem; nanmean(max1 (intersect(sem, find(ID==i))))];
    var1.schema = [var1.schema; nanmean(max1 (intersect(schema, find(ID==i))))];
    
    
    var2.fam = [var2.fam; nanmean(max2 (intersect(fam, find(ID==i))))];
    var2.famw = [var2.famw; nanmean(max2 (intersect(fam_w, find(ID==i))))];
    var2.foils = [var2.foils; nanmean(max2 (intersect(foils, find(ID==i))))];
    var2.fs = [var2.fs; nanmean(max2 (intersect(fs, find(ID==i))))];
    var2.sem = [var2.sem; nanmean(max2 (intersect(sem, find(ID==i))))];
    var2.schema = [var2.schema; nanmean(max2 (intersect(schema, find(ID==i))))];
    
    var3.fam = [var3.fam; nanmean(trap1 (intersect(fam_all, find(ID==i))))];
    var3.unf = [var3.unf; nanmean(trap1 (intersect(unf_all, find(ID==i))))];
    
    var4.fam = [var4.fam; nanmean(trap2 (intersect(fam_all, find(ID==i))))];
    var4.unf = [var4.unf; nanmean(trap2 (intersect(unf_all, find(ID==i))))];
    
    var5.fam = [var5.fam; nanmean(meanSO1 (intersect(fam_all, find(ID==i))))];
    var5.unf = [var5.unf; nanmean(meanSO1 (intersect(unf_all, find(ID==i))))];
    
    var6.fam = [var6.fam; nanmean(meanSO2 (intersect(fam_all, find(ID==i))))];
    var6.unf = [var6.unf; nanmean(meanSO2 (intersect(unf_all, find(ID==i))))];
    
    

end
%% accuracy fam 
fam_acc = [];
unf_acc = [];
FAM_RATING=find(FAMRATING==1);
UNF_RATING=find(UNFRATING==1);
FAM_SURVEY=find(FAMSURVEY==1);
UNF_SURVEY=find(UNFSURVEY==1);
UNF_RATING_W=find(UNFRATING==0);
FAM_RATING_W=find(FAMRATING==0);

unfincorrect=intersect(UNF_RATING_W, UNF_SURVEY);
unfcorrect=intersect(UNF_RATING, UNF_SURVEY);

famincorrect=intersect(FAM_RATING_W, FAM_SURVEY);
famcorrect=intersect(FAM_RATING, FAM_SURVEY);

fam_total=[];
fam_correct=[];
groupID = [];

for i = 1: length (unique(ID))
    j = find(ID==i) ;
    
    groupID = [groupID; group(j(1))];
    
    s = intersect (find (ID ==i), FAM_SURVEY);
    r = intersect (find (ID ==i), UNF_SURVEY);
    
    fam_total = [fam_total; length(s)];
    fam_correct= [fam_correct; length(intersect(famcorrect, s))];
    
    fam_acc = [fam_acc; length(intersect(famcorrect, s))/length(s)];
    unf_acc = [unf_acc; length(intersect(unfcorrect, r))/length(r)];
    
    
end


%% accuracy oddball 
sem_acc = [];
schema_acc = [];
wn_acc = []; 
fs_acc = []; 
fa_odd = []; %only for familiar tunes
fa_foils = []; 
groupID = [];

if fam_rating
    FAM_RATING=find(FAMRATING==1);
    UNF_RATING=find(UNFRATING==1);
    FAM_SURVEY=find(FAMSURVEY==1);
    UNF_SURVEY=find(UNFSURVEY==1);
    UNF_RATING_W=find(UNFRATING==0);
    FAM_RATING_W=find(FAMRATING==0);
    
    famincorrect=intersect(FAM_RATING_W, FAM_SURVEY);
    famcorrect=intersect(FAM_RATING, FAM_SURVEY);
    
    famC = intersect (famcorrect, fam);
    fam_wC = intersect (famcorrect, fam_w);
    semC = intersect (famcorrect, sem);
    schemaC = intersect (famcorrect, schema);
    fam_allC = intersect (famcorrect, fam_all);
    
    sem = semC;
    schema = schemaC;
    fam_w = fam_wC;
    fam = famC;
end

for i = 1: length (unique(ID))
    j = find(ID==i) ;
    
    groupID = [groupID; group(j(1))];
    
    odd1 = intersect (find (ID ==i), sem);
    odd2 = intersect (find (ID ==i), schema);
    odd3 = intersect (find (ID ==i), fam_w);
    odd4 = intersect (find (ID ==i), fs);
    noodd = intersect (find (ID ==i), fam);
    noodd_foils = intersect (find (ID ==i), foils);
    
    sem_acc = [sem_acc; length(find(accuracy (odd1)==1))/length(odd1)];
    schema_acc = [schema_acc; length(find(accuracy (odd2)==1))/length(odd2)];
    wn_acc = [wn_acc; length(find(accuracy (odd3)==1))/length(odd3)];
    fs_acc = [fs_acc; length(find(accuracy (odd4)==1))/length(odd4)];
    

    fa_odd = [fa_odd; length(find(accuracy (noodd)==0))/length(noodd)];
    fa_foils = [fa_foils; length(find(accuracy (noodd_foils)==0))/length(noodd_foils)];
    
end


%% get mean AUCzscore and mean maxzscore for each song for each group 
HC_song_trap1 = []; SD_song_trap1 = []; PNFA_song_trap1=[];  BV_song_trap1 =[]; AD_song_trap1=[]; 
HC_song_max1 = []; SD_song_max1 = []; PNFA_song_max1=[];  BV_song_max1 =[]; AD_song_max1=[]; 

    for j= 1:length(nameSong)
         HC_song_trap1 = [HC_song_trap1; nanmean(trap1 (intersect(find(group==1), find(song==j))))]; 
         SD_song_trap1 = [SD_song_trap1; nanmean(trap1 (intersect(find(group==2), find(song==j))))];
         PNFA_song_trap1 = [PNFA_song_trap1; nanmean(trap1 (intersect(find(group==3), find(song==j))))]; 
         BV_song_trap1 = [BV_song_trap1; nanmean(trap1 (intersect(find(group==4), find(song==j))))]; 
         AD_song_trap1 = [AD_song_trap1; nanmean(trap1 (intersect(find(group==5), find(song==j))))]; 
        
         HC_song_max1 = [HC_song_max1; nanmean(max1 (intersect(find(group==1), find(song==j))))]; 
         SD_song_max1 = [SD_song_max1; nanmean(max1 (intersect(find(group==2), find(song==j))))];
         PNFA_song_max1 = [PNFA_song_max1; nanmean(max1 (intersect(find(group==3), find(song==j))))]; 
         BV_song_max1 = [BV_song_max1; nanmean(max1 (intersect(find(group==4), find(song==j))))]; 
         AD_song_max1 = [AD_song_max1; nanmean(max1 (intersect(find(group==5), find(song==j))))]; 
        

    end
%% average each variable for each subject only if detection fam = 1 

FAM_RATING=find(FAMRATING==1);
UNF_RATING=find(UNFRATING==1);
FAM_SURVEY=find(FAMSURVEY==1);
UNF_SURVEY=find(UNFSURVEY==1);
UNF_RATING_W=find(UNFRATING==0);
FAM_RATING_W=find(FAMRATING==0);

famincorrect=intersect(FAM_RATING_W, FAM_SURVEY);
famcorrect=intersect(FAM_RATING, FAM_SURVEY);

var1.fam = []; var2.fam = []; var3.fam = []; var4.fam = []; var5.fam = []; var6.fam = [];
 var1.famw = []; var1.foils = []; var1.fs = []; var1.sem = []; var1.schema=[];
 var2.famw = []; var2.foils = []; var2.fs = []; var2.sem = []; var2.schema=[];
var3.unf = []; var4.unf = []; var5.unf = []; var6.unf = [];
famC = intersect (famcorrect, fam); 
fam_wC = intersect (famcorrect, fam_w); 
semC = intersect (famcorrect, sem);
schemaC = intersect (famcorrect, schema); 
fam_allC = intersect (famcorrect, fam_all); 
for i=1:length (unique(ID))
    
    var1.fam = [var1.fam; nanmean(max1 (intersect(famC, find(ID==i))))];
    var1.famw = [var1.famw; nanmean(max1 (intersect(fam_wC, find(ID==i))))];
    var1.sem = [var1.sem; nanmean(max1 (intersect(semC, find(ID==i))))];
    var1.schema = [var1.schema; nanmean(max1 (intersect(schemaC, find(ID==i))))];
    

    var3.fam = [var3.fam; nanmean(trap1 (intersect(fam_allC, find(ID==i))))];

    var5.fam = [var5.fam; nanmean(meanSO1 (intersect(fam_allC, find(ID==i))))];
   
    
    

end



%% get max pupil for correct detection of oddballs vs miss trials 
incorrect=find(accuracy==0);
correct= find(accuracy==1);

var1.fam = []; var2.fam = []; 
 var1.famw = []; var1.foils = []; var1.fs = []; var1.sem = []; var1.schema=[];
 var2.famw = []; var2.foils = []; var2.fs = []; var2.sem = []; var2.schema=[];


fam_wC = intersect (correct, fam_w); 
semC = intersect (correct, sem);
schemaC = intersect (correct, schema); 
fsC = intersect (correct, fs); 

fam_wI = intersect (incorrect, fam_w); 
semI = intersect (incorrect, sem);
schemaI = intersect (incorrect, schema); 
fsI = intersect (incorrect, fs); 

groupID = [];

for i=1:length (unique(ID))
            j = find(ID==i) ;
    
    groupID = [groupID; group(j(1))];
    
    var1.fam = [var1.fam; nanmean(max1 (intersect(fam, find(ID==i))))]; % not for fam (no deviant) 
    var1.famw = [var1.famw; nanmean(max1 (intersect(fam_wC, find(ID==i))))];
    var1.sem = [var1.sem; nanmean(max1 (intersect(semC, find(ID==i))))];
    var1.schema = [var1.schema; nanmean(max1 (intersect(schemaC, find(ID==i))))];
%     var1.foils = [var1.foils; nanmean(max1 (intersect(foils, find(ID==i))))];
    var1.fs = [var1.fs; nanmean(max1 (intersect(fsC, find(ID==i))))];
        
%     var2.fam = [var2.fam; nanmean(max2 (intersect(fam, find(ID==i))))];
    var2.famw = [var2.famw; nanmean(max1 (intersect(fam_wI, find(ID==i))))];
    var2.sem = [var2.sem; nanmean(max1 (intersect(semI, find(ID==i))))];
    var2.schema = [var2.schema; nanmean(max1 (intersect(schemaI, find(ID==i))))];
%     var2.foils = [var2.foils; nanmean(max2 (intersect(foilsI, find(ID==i))))];
    var2.fs = [var2.fs; nanmean(max1 (intersect(fsI, find(ID==i))))];
end 

%% get max pupil for correct detection of oddballs vs miss trials if detection fam = 1 only 
incorrect=find(accuracy==0);
correct= find(accuracy==1);

FAM_RATING=find(FAMRATING==1);
UNF_RATING=find(UNFRATING==1);
FAM_SURVEY=find(FAMSURVEY==1);
UNF_SURVEY=find(UNFSURVEY==1);
UNF_RATING_W=find(UNFRATING==0);
FAM_RATING_W=find(FAMRATING==0);

famincorrect=intersect(FAM_RATING_W, FAM_SURVEY);
famcorrect=intersect(FAM_RATING, FAM_SURVEY);

var1.fam = []; var2.fam = []; 
 var1.famw = []; var1.foils = []; var1.fs = []; var1.sem = []; var1.schema=[];
 var2.famw = []; var2.foils = []; var2.fs = []; var2.sem = []; var2.schema=[];

famC_correct = intersect (correct,famcorrect); 

famC_incorrect = intersect (incorrect,famcorrect); 

fam_wC = intersect (fam_w, famC_correct); 
semC = intersect (sem, famC_correct);
schemaC = intersect (schema, famC_correct); 

fam_wMiss= intersect (fam_w, famC_incorrect); 
semMiss = intersect (sem,famC_incorrect);
schemaMiss = intersect (schema,famC_incorrect); 

fam_CR = intersect (fam, famC_correct);
fam_FA = intersect (fam, famC_incorrect);

for i=1:length (unique(ID))
    
    var1.fam = [var1.fam; nanmean(max1 (intersect(fam_CR, find(ID==i))))]; % not for fam (no deviant) 
    var1.famw = [var1.famw; nanmean(max1 (intersect(fam_wC, find(ID==i))))];
    var1.sem = [var1.sem; nanmean(max1 (intersect(semC, find(ID==i))))];
    var1.schema = [var1.schema; nanmean(max1 (intersect(schemaC, find(ID==i))))];
%     var1.foils = [var1.foils; nanmean(max1 (intersect(foils, find(ID==i))))];
%     var1.fs = [var1.fs; nanmean(max1 (intersect(fs, find(ID==i))))];
        
    var2.fam = [var2.fam; nanmean(max1 (intersect(fam_FA, find(ID==i))))];
    var2.famw = [var2.famw; nanmean(max1 (intersect(fam_wMiss, find(ID==i))))];
    var2.sem = [var2.sem; nanmean(max1 (intersect(semMiss, find(ID==i))))];
    var2.schema = [var2.schema; nanmean(max1 (intersect(schemaMiss, find(ID==i))))];
%     var2.foils = [var2.foils; nanmean(max2 (intersect(foilsI, find(ID==i))))];
%     var2.fs = [var2.fs; nanmean(max2 (intersect(fsI, find(ID==i))))];
end 

%% behavioural scores (accuracy for oddball and accuracy for fam) for each song 

[G, nameSongGroup]=findgroups(song(group==3));

mean_accuracy{1,:}=nameSong(nameSongGroup);
mean_accuracy{2,:}=splitapply(@mean,accuracy(group==5),G);
mean_accuracy{3,:}=splitapply(@mean,score(group==5),G);
