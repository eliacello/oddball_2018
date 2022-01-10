function out=timecourse_fam_zscore(type, index_time_playsound, pupil, dur,pupilBasAvg,pupilBasStd, sr )

if nargin==4
    sr=50;
end

% timecourse 0.5s baseline before sound onset
for i=1:length(type)
    if isnan (type (i)) 
        out{1,i}=[];
    
    elseif index_time_playsound(type(i))~= 0 

                   
        baseline{1,i} = mean(pupil(index_time_playsound(type(i))-0.5*sr...
            :index_time_playsound(type(i))));


        %zscore first
        zscore_pupil = (pupil(index_time_playsound(type(i))-0.5*sr...
             :index_time_playsound(type(i))+dur*0.001*sr)-pupilBasAvg)./(pupilBasStd/sqrt(0.5*sr)); 
        zscore_baseline = (baseline{1,i}-pupilBasAvg)./(pupilBasStd/sqrt(0.5*sr));
        %then baseline-correct 
        out{1,i} = zscore_pupil-zscore_baseline;
                    
        
        
    else
        out{1,i}=[]; 
    end
    
end

end


