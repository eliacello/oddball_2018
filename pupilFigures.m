%% Elia Benhamou (8 avril 2020) - elia.benhamou.16@ucl.ac.uk
% create a cell array structure for each diagnosis group 
% and each type of oddball
% plot time series and scatter plot of maximum values for hits and miss trials

addpath /Users/Eliacello/Documents/MATLAB/edf-converter-master
addpath /Users/Eliacello/Dropbox/study1_oddball/analysis_pupil_122019/
addpath /Users/Eliacello/Documents/MATLAB/ICM_Eyelink_Oct2018/pupil_preprocessing_MBB/

filename='subjects_list';
sheet='groups';

[~,~,name_subj] = xlsread([filename '.xlsx'], sheet);

dirPupil='/Users/Eliacello/Documents/EBfiles_STUDY1/results_pupil_TC/final/';
fs              = filesep;

dirPupil2Blocks = '/Users/Eliacello/Documents/EBfiles_STUDY1/results_pupil_TC/final/';

dir_base = '/Users/Eliacello/Documents/EBfiles_STUDY1/';

groups = cell2mat(name_subj(:,2));

dirFinal= [dirPupil2Blocks 'allSubjects'];

%% create MAT: one column per subject

i = 1; % choose diagnostic group (1:controls - 2: SD - 3: PNFA - 4: bvFTD - 5: AD)

index = find (groups ==i);

mat1 = []; MAT1 = []; %famous
mat2 = []; MAT2 = []; %foils
mat3 = []; MAT3 = []; %semantic
mat4 = []; MAT4 = []; %schematic
mat5 = []; MAT5 = []; %white-noise
mat6 = []; MAT6 = []; %foils-schematic

ID = [];
name = [];
lengthID_schema = []; 
lengthID_fs = [];

for i = 1:length(index)

    fam = 0; foils = 0; sem = 0; schema = 0; wn = 0; fs = 0;
    filename = [dirPupil2Blocks name_subj{index(i), 1} '.mat'];
    load (filename);
    
    songToRemove = find(strcmp(TC.ICODD(:, 2), 'Fly.wav')); 
    TC.ICODD (songToRemove,:) = [];
    for j = 1:size(TC.ICODD, 1)
  
        ID = [ID; i];
        name = [name; name_subj(index(i), 1)];
        if strcmp(TC.ICODD{j, 6}, 'famous')==1
            mat1 = [mat1 TC.ICODD{j, 12}];
            fam = fam + 1;
        elseif strcmp(TC.ICODD{j, 6}, 'foils')==1
            mat2 = [mat2 TC.ICODD{j, 12}];
            foils = foils + 1;
        elseif strcmp(TC.ICODD{j, 6}, 'semantic')==1
            mat3 = [mat3 TC.ICODD{j, 12}];
            sem = sem + 1;
        elseif strcmp(TC.ICODD{j, 6}, 'schematic')==1
            mat4 = [mat4 TC.ICODD{j, 12}];
            schema = schema + 1;
        elseif strcmp(TC.ICODD{j, 6}, 'famous_w')==1
            mat5 = [mat5 TC.ICODD{j, 12}];
            wn = wn + 1;
        elseif strcmp(TC.ICODD{j, 6}, 'foils_schema')==1
            mat6 = [mat6 TC.ICODD{j, 12}];
            fs = fs + 1;
        end

    end
    
    lengthID_schema = [lengthID_schema; schema]; 
    lengthID_fs = [lengthID_fs; fs]; 

    % remove trials with std > 3
    badTrials1 = find(nanstd(mat1)>3);
    badTrials2 = find(nanstd(mat2)>3);
    badTrials3 = find(nanstd(mat3)>3);
    badTrials4 = find(nanstd(mat4)>3);
    badTrials5 = find(nanstd(mat5)>3);
    badTrials6 = find(nanstd(mat6)>3);
    
    if ~isempty(badTrials1)        mat1(:,badTrials1) = []; end
    if ~isempty(badTrials2)        mat2(:,badTrials2) = []; end
    if ~isempty(badTrials3)        mat3(:,badTrials3) = []; end
    if ~isempty(badTrials4)        mat4(:,badTrials4) = []; end
    if ~isempty(badTrials5)        mat5(:,badTrials5) = []; end
    if ~isempty(badTrials6)        mat6(:,badTrials6) = []; end
    
    MAT1 = [MAT1 nanmean(mat1(:, end-fam+1+length(badTrials1):end), 2)];
    MAT2 = [MAT2 nanmean(mat2(:, end-foils+1+length(badTrials2):end), 2)];
    MAT3 = [MAT3 nanmean(mat3(:, end-sem+1+length(badTrials3):end), 2)];
    MAT4 = [MAT4 nanmean(mat4(:, end-schema+1+length(badTrials4):end), 2)];
    MAT5 = [MAT5 nanmean(mat5(:, end-wn+1+length(badTrials5):end), 2)];
    MAT6 = [MAT6 nanmean(mat6(:, end-fs+1+length(badTrials6):end), 2)];
    
    length(badTrials1)/fam
    length(badTrials2)/foils
    length(badTrials3)/sem
    length(badTrials4)/schema
    length(badTrials5)/wn
    length(badTrials6)/fs
 
end



%% extract peaks and latencies for each subject
[maxFam, lat1]= max (MAT1(50:150,:)); %0.5 to 2.5 sec 
[maxFoils, lat2] = max (MAT2(50:150,:));
[maxSem, lat3] = max (MAT3(50:150,:));
[maxSchema, lat4] = max (MAT4(50:150,:));
[maxWn, lat5]= max (MAT5(50:150,:));
[maxFs, lat6] = max (MAT6(50:150,:));
lat3 = (lat3'+25)/50;
lat4 = (lat4'+25)/50;
lat5 = (lat5'+25)/50;
lat6 = (lat6'+25)/50;
maxFam=maxFam';
maxSem = maxSem';
maxSchema = maxSchema';
maxWn = maxWn';
maxFoils=maxFoils';
maxFs=maxFs';

%extract auc for each subject
aucfam = trapz (MAT1(50:150,:))'; %0.5 to 2.5 sec
aucfoils = trapz (MAT2(50:150,:))'; %0.5 to 2.5 sec
aucsem = trapz (MAT3(50:150,:))'; %0.5 to 2.5 sec
aucschema= trapz (MAT4(50:150,:))'; %0.5 to 2.5 sec
aucwn = trapz (MAT5(50:150,:))'; %0.5 to 2.5 sec
aucfs = trapz (MAT6(50:150,:))'; %0.5 to 2.5 sec

%% plot all on same graphs with and without shaded errors
% % -0.5 to 3sec

figure();hold on;
% without error bars
x = 1:size(MAT1(1:end-1,:),1);
plot (x, nanmean(MAT1(1:end-1,:),2), '-^','MarkerIndices',1:5:225,'MarkerSize',4, 'Color', 'k'); hold on
plot (x, nanmean(MAT3(1:end-1,:),2), '-+','MarkerIndices',1:5:225,'MarkerSize',4, 'Color', 'r'); hold on
plot (x, nanmean(MAT4(1:end-1,:),2),  '-o','MarkerIndices',1:5:225,'MarkerSize',4, 'Color', 'b'); hold on
plot (x, nanmean(MAT5(1:end-1,:),2), '-x','MarkerIndices',1:5:225,'MarkerSize',4, 'Color', [.49 .18 .56]);

figure();hold on;
% with error bars
shadedErrorBar(1:size(MAT1(1:end-1,:),1),nanmean(MAT1(1:end-1,:),2),nanstd(MAT1(1:end-1,:),0,2)/sqrt(size(MAT1(1:end-1,:),2)),'-k',1); hold on;
shadedErrorBar(1:size(MAT3(1:end-1,:),1),nanmean(MAT3(1:end-1,:),2),nanstd(MAT3(1:end-1,:),0,2)/sqrt(size(MAT3(1:end-1,:) ,2)),'-r',1); hold on;
shadedErrorBar(1:size(MAT4(1:end-1,:),1),nanmean(MAT4(1:end-1,:),2),nanstd(MAT4(1:end-1,:),0,2)/sqrt(size(MAT4(1:end-1,:) ,2)),'-b',1); hold on;
shadedErrorBar(1:size(MAT5(1:end-1,:),1),nanmean(MAT5(1:end-1,:),2),nanstd(MAT5(1:end-1,:),0,2)/sqrt(size(MAT5(1:end-1,:) ,2)),{'-','Color',[.49 .18 .56],'markerfacecolor',[.49 .18 .56]},1);

title ('HC','FontSize',16);
set(gca,'XLim',[0 225],'XTick',[0:25:225],'FontSize',14)
set(gca,'XTickLabel',[-0.5:0.5:4])
ylabel ('pupil [zscore]', 'FontSize',16);

h=line([25 25], [-2 1.5]);
set (h, 'Color', 'black');
set (h, 'LineWidth', 1);
xlabel ('time (s)', 'FontSize',16)

legend ({'no deviant', 'semantic', 'syntactic', 'white-noise'},'FontSize',12);
legend BOXOFF



%% plot same conditions for different groups
figure()
shadedErrorBar(1:size(MAT3(1:end-1,:),1),nanmean(MAT3(1:end-1,:),2),nanstd(MAT3(1:end-1,:),0,2)/sqrt(size(MAT3(1:end-1,:),2)),'-k',1); hold on;
shadedErrorBar(1:size(MAT3SD(1:end-1,:),1),nanmean(MAT3SD(1:end-1,:),2),nanstd(MAT3SD(1:end-1,:),0,2)/sqrt(size(MAT3SD(1:end-1,:) ,2)),'-r',1); hold on;
shadedErrorBar(1:size(MAT3PNFA(1:end-1,:),1),nanmean(MAT3PNFA(1:end-1,:),2),nanstd(MAT3PNFA(1:end-1,:),0,2)/sqrt(size(MAT3PNFA(1:end-1,:) ,2)),'-b',1); hold on;
shadedErrorBar(1:size(MAT3BV(1:end-1,:),1),nanmean(MAT3BV(1:end-1,:),2),nanstd(MAT3BV(1:end-1,:),0,2)/sqrt(size(MAT3BV(1:end-1,:) ,2)),'-b',1); hold on;


%% syntactic in familiar and unfamiliar tunes

figure();hold on;
x = 1:size(MAT1(1:end-1,:),1);
plot (x, nanmean(MAT1(1:end-1,:),2), '-','LineWidth',2, 'Color', 'r'); hold on
plot (x, nanmean(MAT2(1:end-1,:),2), '-','LineWidth',2, 'Color', 'b'); hold on
plot (x, nanmean(MAT4(1:end-1,:),2),  '-.','LineWidth',2, 'Color', 'r'); hold on
plot (x, nanmean(MAT6(1:end-1,:),2), '-.','LineWidth',2, 'Color', 'b');
title ('svPPA','FontSize',16);
set(gca,'XLim',[0 225],'XTick',[0:25:225],'FontSize',12)
set(gca,'XTickLabel',[-0.5:0.5:4])
ylabel ('pupil [zscore]', 'FontSize',14);

h=line([25 25], [-0.8 0.8]);
set (h, 'Color', 'black');
xlabel ('time (s)', 'FontSize',14 )
legend ({'no deviant - familiar', 'no deviant - unfamiliar', 'syntax deviant - familiar ', 'syntax deviant - unfamiliar'},'FontSize',12);
legend BOXOFF



%% plot maximums (scatter plots)
L = length(find (groups==i));

X1=repmat(1.2, 1,L);
X2=repmat(2.2,1,L);
X3=repmat(3.2,1,L);
X4=repmat(4.2,1,L);

STD=[nanstd(maxFam)/sqrt(L) nanstd(maxSem)/sqrt(L) ...
    nanstd(maxSchema)/sqrt(L) nanstd(maxWn)/sqrt(L)];

STDDiff=[nanstd(maxSem-maxFam)/sqrt(L) ...
    nanstd(maxSchema-maxFam)/sqrt(L) nanstd(maxWn-maxFam)/sqrt(L)];

MedTot=[nanmedian(maxFam) nanmedian(maxSem)  nanmedian(maxSchema) nanmedian(maxWn)];
MeanTot=[nanmean(maxFam) nanmean(maxSem)  nanmean(maxSchema) nanmean(maxWn)];

MeanTotDiff=[nanmean(maxSem-maxFam)  nanmean(maxSchema-maxFam) nanmean(maxWn-maxFam)];


figure();
h=bar(MeanTot);

hold on;
e1=errorbar(MeanTot,STD,'.');

e1.Color=[0 0 0];
e1.LineWidth=1;

hold on;
h1=scatter(X1,maxFam,'o');
h2=scatter(X2,maxSem,'o');
h3=scatter(X3,maxSchema,'o');
h4=scatter(X4,maxWn,'o');

% h1Diff=scatter(X1,(maxSem-maxFam),'o');
% h2Diff=scatter(X2,(maxSchema-maxFam),'o');
% h3Diff=scatter(X3,(maxWn-maxFam),'o');



hold off;
aX.XAxis.TickValues=0:1:5; %change to 5
set(gca, 'YLim', [0 2]);%change to [0 2] (or [-1 1.5])
set(gca, 'XLim', [0 5]); %change to 5
set(gca, 'XTick', [1:1:4]); %change to 4
set(gca, 'XTickLabel', {'no deviant', 'semantic',  'syntactic', 'white-noise'});
set(h,'FaceColor',[.9 .9 .9]);
ylabel ('maximum pupil (deviant)');%'maximum PDR (deviant) - maximum PDR (no deviant) '
title ('AD');


%% plot unf fs
% TC
figure();hold on;
x = 1:size(foils(1:end-51,:),1);
p = plot (x, nanmean(foils(1:end-51,:),2), x, nanmean(fs(1:end-51,:),2));
title ('AD');
set(gca,'XLim',[0 175],'XTick',[0:25:175])
set(gca,'XTickLabel',[-0.5:0.5:3])
ylabel ('PDR [zscore]');
h=line([25 25], [-0.6 0.8]);
set (h, 'Color', 'black');
xlabel ('time (s)')

% max
L = length(find (groups==i));

X1=repmat(1.2, 1,L);
X2=repmat(2.2,1,L);

STD=[nanstd(maxFoils)/sqrt(L)  ...
    nanstd(maxFs)/sqrt(L)];
MeanTot=[nanmean(maxFoils) nanmean(maxFs)];

figure();
h=bar(MeanTot);

hold on;
e1=errorbar(MeanTot,STD,'.');

e1.Color=[0 0 0];
e1.LineWidth=1;

hold on;

h1=scatter(X1,maxFoils,'o');
h2=scatter(X2,maxFs,'o');

hold off;
aX.XAxis.TickValues=0:1:3; %change to 5
set(gca, 'YLim', [-1 1]);%change to [0 2]
set(gca, 'XLim', [0 3]); %change to 5
set(gca, 'XTick', [1:1:2]); %change to 4
set(gca, 'XTickLabel', {'foils',  'foils-schematic'});
set(h,'FaceColor',[.9 .9 .9]);
ylabel ('maximum PDR (deviant) - maximum PDR (no deviant) ');
title ('AD');


%% plot with shaded errors
figure();hold on;
shadedErrorBar(1:size(fam,1),mean(fam,2),std(fam,0,2)/sqrt(size(fam,2)),'-r',1);
shadedErrorBar(1:size(unf,1),mean(unf,2),std(unf,0,2)/sqrt(size(unf,2)),'-b',1);
title ('AD: PDR to melody familiarity');
set(gca,'XTickLabel',[-0.5 0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5])

figure();hold on;
shadedErrorBar(1:size(foils,1),nanmean(foils,2),nanstd(foils,0,2)/sqrt(size(foils,2)),'-b',1);
shadedErrorBar(1:size(fs,1),nanmean(fs,2),nanstd(fs,0,2)/sqrt(size(fs ,2)),'-r',1);
title ('AD: PDR to schematic deviants in unfamiliar melodies');
set(gca,'XTickLabel',[-0.5 0.5 1.5 2.5 3.5 4.5])

figure();hold on;
shadedErrorBar(1:size(famous,1),nanmean(famous,2),nanstd(famous,0,2)/sqrt(size(famous,2)),'-b',1);
shadedErrorBar(1:size(sem,1),nanmean(sem,2),nanstd(sem,0,2)/sqrt(size(sem ,2)),'-r',1);
title ('AD: PDR to semantic deviants in familiar melodies');
set(gca,'XTickLabel',[-0.5 0.5 1.5 2.5 3.5 4.5])

figure();hold on;
shadedErrorBar(1:size(famous,1),nanmean(famous,2),nanstd(famous,0,2)/sqrt(size(famous,2)),'-b',1);
shadedErrorBar(1:size(schema,1),nanmean(schema,2),nanstd(schema,0,2)/sqrt(size(schema ,2)),'-r',1);
title ('AD: PDR to schematic deviants in familiar melodies');
set(gca,'XTickLabel',[-0.5 0.5 1.5 2.5 3.5 4.5])

figure();hold on;
shadedErrorBar(1:size(famous,1),nanmean(famous,2),nanstd(famous,0,2)/sqrt(size(famous,2)),'-b',1);
shadedErrorBar(1:size(famW,1),nanmean(famW,2),nanstd(famW,0,2)/sqrt(size(famW ,2)),'-r',1);
title ('AD: PDR to white-noise deviants in familiar melodies');
set(gca,'XTickLabel',[-0.5 0.5 1.5 2.5 3.5 4.5])


