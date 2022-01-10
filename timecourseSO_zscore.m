function out=timecourseSO_zscore(type, index_time_playsound, index_time_oddsound, pupil,pupilBasAvg,pupilBasStd, sr)

if nargin==3
    sr=50;
end

% timecourse 0.5s baseline before sound onset
for i=1:size(type,1)
    for j=1:size(type,2)
       if isnan(type(i,j))
           out{i,j}=[]; 
       
           
       elseif index_time_playsound(type(i,j))~= 0 
            baseline{i,j} = mean(pupil(index_time_playsound(type(i,j))-0.5*sr...
            :index_time_playsound(type(i,j))));
        
        
            %zscore first
            zscore_pupil = (pupil(index_time_playsound(type(i,j))-0.5*sr...
             :index_time_oddsound(type(i,j)))-pupilBasAvg)./(pupilBasStd/sqrt(0.5*sr)); 
            zscore_baseline = (baseline{i,j}-pupilBasAvg)./(pupilBasStd/sqrt(0.5*sr));
            %then baseline-correct 
            out{1,i} = zscore_pupil-zscore_baseline;

        
    else
        out{i,j}=[]; 
    end
    end
end

end




