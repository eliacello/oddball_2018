%% Sijia Zhao (23 June 2020) sijia.zhao@psy.ox.ac.uk

% Run the time-series permutation test over five groups (with different
% sample sizes) and plot the data.

% close all; clear;
addpath /Users/Eliacello/Dropbox/study1_oddball/analysis_pupil_122019/permutationTest/
addpath(genpath('functions'));

% load('data.mat');

data{1} = MAT5HC;
data{2} = MAT5AD;
data{3} = MAT5BV;
data{4} = MAT5SD;
data{5} = MAT5PNFA;

% data{1} = MATAD;
% data{2} = MAT3BV;

% data{1} = MAT1HC;
% data{2} = MAT3HC;
% data{3} = MAT4HC;
% data{4} = MAT5HC;

groups = {'HC','AD','bvFTD','svPPA','nfvPPA'};
% groups = {'no deviant','semantic','syntactic','acoustic'};

%between_groups
stats_pair = [...
    1,2; 1,3; 1,4; 1,5;...
    2,3; 2,4; 2,5;...
    3,4; 3,5;
    4,5];

%within-groups
% stats_pair = [...
%     1,2; 1,3; 1,4; ...
%     2,3; 2,4; ...
%     3,4];

samfreq = 1/50;
timeaxis = -0.5:samfreq:4;

colourmap_group = cbrewer('qual', 'Dark2',5);

yl = [-1,0.5];

xl = [-0.5,4];

%% Plot all groups' PDR together
figure(1); clf;
hold on;
g = 1; 
   pdr = data{g}';
    nsub = size(pdr,1);
    a = nanmean(pdr,1);
    plot(timeaxis,a,'LineWidth',1.5,'Color','black')
for g = 2:numel(groups)
    pdr = data{g}';
    nsub = size(pdr,1);
    a = nanmean(pdr,1);
    plot(timeaxis,a,'LineWidth',1.5,'Color',colourmap_group(g,:));
end
% legend(groups, 'FontSize',14);
% legend BOXOFF

g = 1;
   pdr = data{g}';
    nsub = size(pdr,1);
    a = nanmean(pdr,1);
    b = nanstd(pdr)/sqrt(nsub);
    
    curve1 = a'+b';
    curve2 = flipud(a'-b');
    X = [timeaxis'; flipud(timeaxis')];
    Y = [curve1; curve2];
    figfill = fill(X',Y','black','edgecolor','none','facealpha',0.2);


for g = 2:numel(groups)
    pdr = data{g}';
    nsub = size(pdr,1);
    
    a = nanmean(pdr,1);
    b = nanstd(pdr)/sqrt(nsub);
    
    curve1 = a'+b';
    curve2 = flipud(a'-b');
    X = [timeaxis'; flipud(timeaxis')];
    Y = [curve1; curve2];
    figfill = fill(X',Y',colourmap_group(g,:),'edgecolor','none','facealpha',0.2);
end
% legend(groups,'Location','northwest');


% add diff lines: 
% syntactic: 
plot(timeaxis,diff(1,:),'LineWidth',16,'Color',colourmap_group(2,:)); % AD vs controls

plot(timeaxis,diff(2,:)*1.2,'LineWidth',16,'Color',colourmap_group(3,:)+0.1); % BV vs HC
% plot(timeaxis,diff(5,:)*1.4,'LineWidth',16,'Color',colourmap_group(3,:)+0.1); % BV vs AD
plot(timeaxis,diff(9,:)*1.4,'LineWidth',16,'Color',colourmap_group(3,:)+0.2); % BV vs PNFA

% plot(timeaxis,diff(3,:)*1.4,'LineWidth',16,'Color',colourmap_group(4,:)-0.1); % SD vs HC
plot(timeaxis,diff(6,:)*1.6,'LineWidth',16,'Color',[0.98    0.25   0.60]); % SD vs AD
% plot(timeaxis,diff(10,:)*1.6,'LineWidth',16,'Color',[0.98    0.50   0.70]); % SD vs PNFA

plot(timeaxis,diff(4,:)*1.8,'LineWidth',16,'Color',colourmap_group(5,:)); % PNFA vs HC



ylim([-1.8 1]);
xlim(xl);
xlabel('Time [s]');
ylabel('Z-score');

% title('group average, shaded area = 1 SEM');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5]);
saveas(gcf,'/Users/Eliacello/Dropbox/study1_oddball/oddball_Brain/permutation_pupil/semantic_all.png');

%% Run permutation test for each pair
figure(2);clf;
diff = [];
for i = 1:size(stats_pair,1)
    
    subplot(2,5,i); hold on;
% subplot(2,3,i); hold on;
    
    g1 = stats_pair(i,1);
    g2 = stats_pair(i,2);
    
    cond1 = data{g1}';
    cond2 = data{g2}';
    dataB = bootstrap(cond1) - bootstrap(cond2); %default number of iteractions = 10000   
    %%% If it's repeated measure you should do:
%         dataB = bootstrap(cond1-cond2);
    %%%
    % bonferroni correction: replace 0.05 by 0.05/n with n: number of
    % comparison 
     
    nb_comp = 10; % 10 for 5 groups - 6 for 4 conditions 
    s = findSigDiff(dataB, 0.05/nb_comp);   %dataB contains all p-values  
    statlimit = findStatLimit([-.5 0.25],s,timeaxis); % check baseline if there is any significant cluster
    [s_corrected,~] = bootstrapCluster(s,statlimit,timeaxis); % remove any cluster whose length < the sig. cluster during baseline
    s = s_corrected; 
    %     s = [s_corrected NaN(1,56)];
  
    pdr = cond1;
    nsub = size(pdr,1);
    a = nanmean(pdr,1);
    b = nanstd(pdr)/sqrt(nsub);
    plot(timeaxis,a,'LineWidth',1.5,'Color',colourmap_group(g1,:));
    
    %add error envelope: 
    curve1 = a'+b';
    curve2 = flipud(a'-b');
    X = [timeaxis'; flipud(timeaxis')];
    Y = [curve1; curve2];
    figfill = fill(X',Y',colourmap_group(g1,:),'edgecolor','none','facealpha',0.2);

    
    pdr = cond2;
    nsub = size(pdr,1);
    a = nanmean(pdr,1);
    b = nanstd(pdr)/sqrt(nsub);
    plot(timeaxis,a,'LineWidth',1.5,'Color',colourmap_group(g2,:));
    
    %add error envelope: 
    curve1 = a'+b';
    curve2 = flipud(a'-b');
    X = [timeaxis'; flipud(timeaxis')];
    Y = [curve1; curve2];
    figfill = fill(X',Y',colourmap_group(g2,:),'edgecolor','none','facealpha',0.2);

    
    ylim(yl);
    xlim(xl);
    
    ab = ylim;
    aa= 0.025*(max(ab)-min(ab));
    dist = 0.05*(max(ab)-min(ab));
    Ystat = ab(1) + aa;
    Ystat = Ystat + dist;
    plot(timeaxis,Ystat*abs(s),'LineWidth',2.5,'Color','k');
    
    legend({groups{g1}, groups{g2}},'Location','northoutside');
    
    % save Ystat*abs(s) in separate array to plot in main graph 
    diff = [diff; Ystat*abs(s)];
end
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8]);
% saveas(gcf,'/Users/Eliacello/Dropbox/study1_oddball/oddball_Brain/permutation_pupil/acoustic_between.png');
