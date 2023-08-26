function [coor1,vtol,all_violations] = adaptive_embedding_refiner(coor0,D,vtol)
% coor1 = adaptive_embedding_refiner(coor0,D,vtol)
%
% A set of coordinates obtained from embedding is iteratively changed until
% it conforms to given lower and upper distance bounds within tolerance
% vtol, the function increases and returns vtol, if required convergence
% cannot be achieved
%
% (c) G.Jeschke, 2023

frac_viol = 0.25;

lower_bounds = D;
upper_bounds = D;
[C,~] = size(D);
for c1 = 1:C-1
    for c2 = c1+1:C
        lower_bounds(c1,c2) = lower_bounds(c1,c2) - (1-frac_viol)*vtol;
        lower_bounds(c2,c1) = lower_bounds(c2,c1) - (1-frac_viol)*vtol;
        upper_bounds(c1,c2) = upper_bounds(c1,c2) + (1-frac_viol)*vtol;
        upper_bounds(c2,c1) = upper_bounds(c2,c1) + (1-frac_viol)*vtol;
    end
end
[lower_dev,upper_dev]=test_embed(coor0,lower_bounds,upper_bounds);
violation_sum = single(sum(sum(lower_dev.^2+upper_dev.^2)));
iter=round(violation_sum);
if iter<1000
    iter=1000;
end
if iter>50000
    iter=50000;
end

if ~exist('vtol','var') || isempty(vtol)
    vtol = 2.5;
end

lambda=linspace(1/iter,1,iter);

[~,n]=size(lower_bounds);

violation0 = sqrt(max(max(lower_dev.^2+upper_dev.^2)));
while violation0 > vtol
    coor=coor0;
    coor1=coor0;
    lower_bounds = D;
    upper_bounds = D;
    [C,~] = size(D);
    for c1 = 1:C-1
        for c2 = c1+1:C
            lower_bounds(c1,c2) = lower_bounds(c1,c2) - (1-frac_viol)*vtol;
            lower_bounds(c2,c1) = lower_bounds(c2,c1) - (1-frac_viol)*vtol;
            upper_bounds(c1,c2) = upper_bounds(c1,c2) + (1-frac_viol)*vtol;
            upper_bounds(c2,c1) = upper_bounds(c2,c1) + (1-frac_viol)*vtol;
        end
    end
    [lower_dev,upper_dev]=test_embed(coor0,lower_bounds,upper_bounds);
    violation0 = sqrt(max(max(lower_dev.^2+upper_dev.^2)));
    all_violations = zeros(iter,1);
    for k=1:iter
        if mod(k,100)==0
            [lower_dev,upper_dev]=test_embed(coor,lower_bounds,upper_bounds);
            violation = sqrt(max(max(lower_dev.^2+upper_dev.^2)));
            all_violations(k) = violation;
            if violation < vtol || violation > violation0
                k = k - 1; %#ok<FXSET>
                break;
            end
            violation0 = violation;
            % fprintf(1,'%s%i%s%5.3f\n','Iteration ',k,' Violation number: ',violations);
        end
        dmat=coor2dmat(coor);
        ndmat=dmat+(dmat==0);
        [lower_dev,upper_dev]=test_embed(coor,lower_bounds,upper_bounds);
        ndev=(lower_dev-upper_dev)./ndmat;
        if sum(sum(isnan(ndev)))
            k = k - 1; %#ok<FXSET>
            break
        end
        for p=1:n
            cvec=zeros(1,3);
            for pp=1:n
                cvec=cvec+ndev(p,pp)*(coor(p,:)-coor(pp,:));
            end
            coor1(p,:)=coor(p,:)+lambda(k)*cvec;
        end
        coor=coor1;
    end
    if violation0 > (1-frac_viol)*vtol
        vtol = vtol + 0.1;
    end
end
all_violations = all_violations(1:k);
all_violations = all_violations(all_violations > 0); 
