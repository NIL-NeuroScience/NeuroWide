function [h,p] = f_lme(MouseID,g1,g2,alpha)

N = numel(MouseID);

col_MouseID = categorical([MouseID,MouseID])';
col_Group = categorical([repmat({'g1'},1,N),repmat({'g2'},1,N)])';
col_Data = cat(1,g1,g2);

G = size(col_Data,2);

p = zeros(G,1);

for i = 1:G
    Data = col_Data(:,i);
    T = table(col_MouseID,col_Group,Data);
    
    lme = fitlme(T,'Data ~ col_Group + (1|col_MouseID)');
    coeff = lme.Coefficients;
    p(i) = coeff.pValue(strcmp(coeff.Name,'col_Group_g2'));
end
h = p < alpha;

end
