function out=timecourseAUC_zscore(index_time_playsound, index_time_oddsound, pupil,pupilBasAvg,pupilBasStd, sr)

if nargin==3
    sr=50;
end

% timecourse 0.5s baseline before sound onset
for i=1:length(index_time_playsound)
    
        baseline{1,i} = mean(pupil(index_time_playsound(i)-0.5*sr...
            :index_time_playsound(i)));

        %zscore first
        zscore_pupil = (pupil(index_time_playsound(i)-0.5*sr...
             :index_time_oddsound(i))-pupilBasAvg)./(pupilBasStd);  %divide by sqrt(0.5*sr) if standard error 
        zscore_baseline = (baseline{1,i}-pupilBasAvg)./(pupilBasStd);
        %then baseline-correct 
         out{1,i} = zscore_pupil-zscore_baseline;
    
end




