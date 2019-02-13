function r = regu(v,p,order,method)
% save('./temp/regu.mat', 'v', 'p', 'order', 'method');


% clear
% clc
% load('./temp/regu.mat');

v(v>1) = 1;
v(v<0) = 0;

switch method
    case 'cvpr07_spectral_matting'
        assert(length(p)<=2);
        if length(p)==1
            p0 = p; p1 = p;
        else
            p0 = p(1); p1 = p(2);
        end
        
        assert(0<=p0 && p0<=1);
        assert(0<=p1 && p1<=1);
        
        if order == 0
            r = norm(v,p0)^p0 + norm(1-v,p1)^p1;
        elseif order == 1
            r = p0*v.^(p0-1)-p1*(1-v).^(p1-1);
            r(r==inf) = 0;
            r(r==-inf) = 0;
        end
        assert(isreal(r));
    case 'cvpr07_spectral_matting+negative_1_norm_v2'
       
        r1 = regu(v,p,order,'cvpr07_spectral_matting');
        r2 = regu(1-v,1,order,'negative_p_norm');
        r = r1+0.1*r2;
        
    case 'cvpr07_spectral_matting+negative_2_norm'
       
        r1 = regu(v,p,order,'cvpr07_spectral_matting');
        r2 = regu(1-v,2,order,'negative_p_norm');
        r = r1+0.1*r2;    
        
    case 'cvpr07_spectral_matting+negative_1_norm'
        assert(length(p)<=2);
        if length(p)==1
            p0 = p; p1 = p;
        else
            p0 = p(1); p1 = p(2);
        end
        
        assert(0<=p0 && p0<=1);
        assert(0<=p1 && p1<=1);
        
        if order == 0
            v = 1 - v;     
            pp = 1;
            r1 = -norm(v,pp);

            v = 1 - v;
            r = norm(v,p0)^p0 + norm(1-v,p1)^p1+0.1*r1;
        elseif order == 1
            v = 1 - v;
            pp = 1;
            if norm(v)<1e-6
                n = length(v);
                r1 = zeros(n,1);
            else
                r1 = -1*v.^(pp-1)/(sum(v.^pp)^(1-1/pp));
            end
            
            v = 1 - v;
            r = p0*v.^(p0-1)-p1*(1-v).^(p1-1) + 0.1*r1;
            r(r==inf) = 0;
        end    
        assert(isreal(r));
    case 'negative_p_norm'
%         assert(1<=p)
        if order == 0
            r = -norm(v,p);
        elseif order == 1
            if norm(v)<1e-6
                n = length(v);
                r = zeros(n,1);
            else
                r = -1*v.^(p-1)/(sum(v.^p)^(1-1/p));
                r(r==inf) = 0;
                r(r==-inf) = 0;
            end
        end
        
    case 'positive_1_norm'
        p = 1;
        if order == 0
            r = norm(v,1);
        elseif order == 1
            n = length(v);
            if norm(v)<1e-6
                r = zeros(n,1);
            else
                r = ones(n,1);
            end
        end
        
    case 'norm1_norm2'
        if order == 0
            r = norm(v,1) - norm(v,2);
        elseif order == 1
            n = length(v);
            if norm(v)<1e-6
                r = ones(n,1);
            else
                r = ones(n,1) - v/norm(v,2);
            end
        end
        
end