function mmat=metric_matrix(dmat)
% Compute metric matrix from squared distance matrix
% (c) G. Jeschke, 2008


n=size(dmat,1);
dmat2=dmat.^2;
mvec=zeros(1,n);
hvec=ones(1,n);
bckg=sum(sum(dmat2))/(2*n^2);
for k=1:n,
    mvec(k)=sum(dmat2(k,:))/n-bckg;
end;
mmat=(kron(mvec',hvec)+kron(hvec',mvec)-dmat2)/2;
% mmat=zeros(n,n);
% for i=1:n,
%     for j=1:n,
%         mmat(i,j)=(mvec(i)+mvec(j)-dmat2(i,j))/2;
%     end;
% end;
% disp(sprintf('%s%6.4f','Maximum deviation between metric matrices: ',max(max(abs(mmat-mmat0)))));
