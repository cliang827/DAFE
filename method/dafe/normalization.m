function F = normalization(F, range, center, method)
% save('./temp/normalization.mat', 'F', 'range', 'center', 'method');
% purpose: normalize f into range=[-1,1] and center at center=0

% clear
% clc
% load('./temp/normalization.mat');

switch method
    case 'mean-priority'
        k = size(F,2);
        for i=1:k
            f = F(:,i);
            f = f - mean(f);
            f1 = abs(max(f));
            f2 = abs(min(f));
            sigma = 1/max(f1, f2);
            f = f*sigma;

            assert(max(f)<=range(2) && min(f)>=range(1) && norm(mean(f)-center)<1e-4)
            F(:,i) = f;
        end
        
    case 'range-priority'
        k=size(F,2);
        for i=1:k
            f=F(:,i);
            f = ((f-min(f))./(max(f)-min(f)));
            F(:,i)=f;
        end
end
