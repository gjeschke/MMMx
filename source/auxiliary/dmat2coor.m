function [coor,sevals]=dmat2coor(dmat)
% Cartesian embedding of a distance matrix
coor=[];
n=size(dmat,1);
mmat=metric_matrix(dmat);
[evec,evalsm]=eig(mmat); %,'nobalance'); % diagonalization of the metric matrix
evec=evec';
evalsm=evalsm';
evals=diag(evalsm); % vector of eigenvalues
[sevals,index]=sort(evals,1,'descend'); % sorted eigenvalues
for k=1:3,
    if sevals(k)<=0,
        fprintf(1,'%s%i%s','### Error: eigenvalue ',k,' not positive ###');
        return
    end;
end;
% disp(sprintf('%s%6.4f','Ratio of 4th to 3rd eigenvalue of the metric matrix: ',sevals(4)/sevals(3)));
% disp(sqrt(sevals(1:3)));
% disp(evec(index(1),1:8));
% disp(evec(index(2),1:8));
% disp(evec(index(3),1:8));
coor=zeros(n,3);
for k=1:3,
    coor(:,k)=sqrt(sevals(k))*real(evec(index(k),:))';
end;
