function []=SetOpacity(es,ed,F,F2)


for j = 1:size(F,1)
    for i=1:size(F,2)
        if ~isempty(F(j,i))
            if i <= es.Value
                set(F{j,i},'facealpha',1)
                set(F2{j,i},'facealpha',1)
            else
                set(F{j,i},'facealpha',0,'EdgeAlpha',0.25)
                set(F2{j,i},'facealpha',0,'EdgeAlpha',0.25)
            end
        end
    end
end


end