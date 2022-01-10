function out=timecourse_zscore(type, position, index_time_oddball, pupil, pupilBasOddAvg, pupilBasOddStd, sr)

if nargin==3
    sr=50;
end

for j=1:length(position)
        type_position{1,j}=intersect(type, position{1,j});
        for i=1:length(type_position{1,j})
            if isempty(type)==0 
                baseline{i,j} = mean(pupil(index_time_oddball(type_position{1,j}(i))-0.5*sr...
                    :index_time_oddball(type_position{1,j}(i))));
                if (index_time_oddball(type_position{1,j}(i))+4*sr) < length (pupil)
                    
                    %zscore first (baseline as well) 
                    zscore_pupil = (pupil(index_time_oddball(type_position{1,j}(i))-0.5*sr...
                    :index_time_oddball(type_position{1,j}(i))+4*sr)-pupilBasOddAvg)./(pupilBasOddStd/sqrt(0.5*sr)); 
                    zscore_baseline = (baseline{i,j}-pupilBasOddAvg)./(pupilBasOddStd/sqrt(0.5*sr));
                    %then baseline-correct 
                out.position{i,j} = zscore_pupil-zscore_baseline;
                
                else out.position{i,j}=[];
                    
                end 
            else
                out.position{i,j}=[];
            end
        end

end


