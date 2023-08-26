function [coor1,err,violations,all_violations] = embedding_refiner(coor0,lower_bounds,upper_bounds,iter,vtol)
% coor1 = embedding_refiner(coor0,lower_bounds,upper_bounds,iter)
%
% A set of coordinates obtained from embedding is iteratively changed until
% it conforms to given lower and upper distance bounds within tolerance
% vtol, the function returns after iter iterations, if required convergence
% cannot be achieved
%
% if the procedure does not converge, an error code 1 is given back,
% otherwise the error code is 0
%
% (c) G.Jeschke, 2008

err=0;

[lower_dev,upper_dev]=test_embed(coor0,lower_bounds,upper_bounds);
violations= single(sum(sum(lower_dev.^2+upper_dev.^2)));
% disp(sprintf('%s%5.3f','Bound refiner start. Violation number: ',violations));

if nargin<4
    iter=round(violations);
    if iter<1000
        iter=1000;
    end
    if iter>50000
        iter=50000;
    end
end

if nargin<5
    vtol=single(1e-6);
end

lambda=linspace(1/iter,1,iter);

[~,n]=size(lower_bounds);

coor=coor0;
coor1=coor0;
violations = 1e12;
all_violations = zeros(iter,1);
for k=1:iter
    if mod(k,100)==0
        [lower_dev,upper_dev]=test_embed(coor,lower_bounds,upper_bounds);
        violations0 = violations;
        violations=sum(sum(lower_dev.^2+upper_dev.^2));
        all_violations(k) = violations;
        if violations<vtol || violations > violations0
            k = k - 1; %#ok<FXSET>
            break;
        end
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
all_violations = all_violations(1:k);
if violations<vtol
    % disp(sprintf('%s%i%s','Refinement converged after ',k,' iterations.'));
else
    % disp(sprintf('%s%i%s','Refinement not converged after ',k,' iterations.'));
    err=1;
end
