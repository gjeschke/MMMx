function [lower_bounds,upper_bounds,err]=triangle_bound_smoothing(lower_bounds,upper_bounds,silent)
% [lower_bounds,upper_bounds,err]=triangle(lower_bounds,upper_bounds)
%
% Triangle bound smoothing
% algorithm as in G. M. Crippen, T. F. Havel, Distance Geometry and
% Molecular Conformation, Research Studies Press, Taunton, 1988, p. 252 f.
%
% error output to the command window is suppressed if input argument silent 
% is provided and is not zero
%
% given initial bounds on distances lower_bound and upper_bounds, the
% function ensures that the smoothed bounds fulfil the triangle
% inequalities
%
% Error codes err:
% 0   no error
% 1   erroneous bounds 
% 2   at least one of the bound matrices is not square
% 3   bound matrices do not have the same size
%
% (c) G. Jeschke 2008

if ~exist('silent','var'), silent= true; end

err=0;
[m,n]=size(lower_bounds);
if m~=n,
    if ~silent, add_msg_box('ERROR: Lower bound matrix not square.'); end;
    err=2;
    return;
end;
[m2,n2]=size(upper_bounds);
if m2~=n2,
    if ~silent, add_msg_box('ERROR: Upper bound matrix not square.'); end;
    err=2;
    return;
end;
if n~=n2,
    if ~silent, add_msg_box('ERROR: Bound matrices have different size.'); end;
    err=3;
    return;
end;

for k=1:n,
    for kk=1:n-1,
        for kkk=kk+1:n,
            if upper_bounds(kk,kkk) > upper_bounds(kk,k)+upper_bounds(k,kkk)
                upper_bounds(kk,kkk)=upper_bounds(kk,k)+upper_bounds(k,kkk);
            end;
            if lower_bounds(kk,kkk) < lower_bounds(kk,k)-upper_bounds(k,kkk)
               lower_bounds(kk,kkk)=lower_bounds(kk,k)-upper_bounds(k,kkk);
            elseif lower_bounds(kk,kkk) < lower_bounds(kkk,k)-upper_bounds(k,kk)
                lower_bounds(kk,kkk) = lower_bounds(kkk,k)-upper_bounds(k,kk);
            end;
            if lower_bounds(kk,kkk)>upper_bounds(kk,kkk)
                if ~silent,
                    add_msg_box('ERROR: Erroneous bounds.');
                    add_msg_box(sprintf('%s%i%s%i%s%5.2f%s%5.2f','(',kk,',',kkk,'): Lower: ',lower_bounds(kk,kkk),', Upper: ',upper_bounds(kk,kkk)));
                end;
                err=1;
                return
            end;
        end;
    end;
end;

% Make matrix symmetric again
for k=1:n-1,
    for kk=k+1:n,
        lower_bounds(kk,k)=lower_bounds(k,kk);
        upper_bounds(kk,k)=upper_bounds(k,kk);
    end;
end;

