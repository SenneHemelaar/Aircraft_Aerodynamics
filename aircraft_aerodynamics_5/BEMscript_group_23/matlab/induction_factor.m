%% Find induction for a certain thrust coefficient 
function [a] = induction_factor(CT)
    a = zeros(size(CT));
    CT1=1.816;
    CT2=2*sqrt(CT1)-CT1;
    a(CT>=CT2) = 1 + (CT(CT>=CT2)-CT1)/(4*(sqrt(CT1)-1));
    a(CT<CT2) = 0.5-0.5*sqrt(1-CT(CT<CT2));
end %fun

