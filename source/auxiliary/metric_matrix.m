function mmat = metric_matrix(dmat)
% Compute metric matrix from squared distance matrix
% see also: G. Young, A. S. Householder, Psychometrika, 1938, 3(1), 19-21 
% (c) G. Jeschke, 2008


n=size(dmat,1);
dmat2=dmat.^2;
mvec=zeros(1,n);
hvec=ones(1,n);
bckg=sum(sum(dmat2))/(2*n^2);
for k=1:n
    mvec(k)=sum(dmat2(k,:))/n-bckg;
end
mmat=(kron(mvec',hvec)+kron(hvec',mvec)-dmat2)/2;
