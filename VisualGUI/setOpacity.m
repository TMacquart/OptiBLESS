function []=setOpacity(es,ed,Fig,F)

% keyboard
% 

for j = 1:size(F,1)
    for i=1:size(F,2)
        if ~isempty(F(j,i))
            if i <= es.Value
                set(F{j,i},'facealpha',1)
            elseif i < es.Value +1
                set(F{j,i},'facealpha',es.Value-floor(es.Value))
            else
                set(F{j,i},'facealpha',0)
            end
        end
    end
end

end