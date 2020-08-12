function c=cross_rowvec(a,b)
% A fast cross product that works only for two three-element row vectors 

c = [a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)];