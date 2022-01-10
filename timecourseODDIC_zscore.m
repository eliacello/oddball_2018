function out=timecourseODDIC_zscore(index_time_oddsound, pupil,pupilBasOddAvg,pupilBasOddStd, sr)

if nargin==3
    sr=50;
end

% timecourse 0.5s baseline before sound onset
for i=1:length(index_time_oddsound)
    
        baseline{1,i} = mean(pupil(index_time_oddsound(i)-0.5*sr...
            :index_time_oddsound(i)));
         if (index_time_oddsound(i)+4*sr) < length (pupil)
             %zscore first
             zscore_pupil = (pupil(index_time_oddsound(i)-0.5*sr...
                 :index_time_oddsound(i)+4*sr)-pupilBasOddAvg)./(pupilBasOddStd); %standard error: divide by sqrt(36)
             zscore_baseline = (baseline{1,i}-pupilBasOddAvg)./(pupilBasOddStd);
             %then baseline-correct \
             out{1,i} = zscore_pupil-zscore_baseline;

         else out{1,i} = [];
         end 

    
end