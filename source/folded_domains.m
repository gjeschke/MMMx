function domains = folded_domains(pae,options)

if ~exist('options','var') || isempty(options) 
    options.threshold = 5;
end

if ~isfield(options,'mindomain')
    options.mindomain = 25;
end

if ~isfield(options,'figures')
    options.figures = false;
end

[nres,~] = size(pae);
resax = 1:nres;

if options.figures
    figure; hold on
    plot(1,1,'k.');
    plot(nres,nres,'k.');
    image(resax,resax,pae,'CDataMapping','scaled','ButtonDownFcn',@(hObject,eventdata,handles)axes_rmsd_ButtonDownFcn);
    curr_axis = gca;
    curr_axis.YDir = 'normal';
    colorbar;
    axis tight
    xlabel('Residue');
    ylabel('Residue');
    title('PAE matrix');
    axis equal
end

% convert PAE matrix to binary disorder matrix 
disorder = pae;
disorder(pae <= options.threshold) = 0;
disorder(pae > options.threshold) = 1;

% if PAE for a pair indicates order for one of the two alignments, assume
% order
for r1 = 1:nres-1
    for r2 = r1+1:nres
        if disorder(r1,r2)*disorder(r2,r1) == 0
            disorder(r1,r2) = 0;
            disorder(r2,r1) = 0;
        end
    end
end

if options.figures
    figure; hold on
    plot(1,1,'k.');
    plot(nres,nres,'k.');
    image(resax,resax,disorder,'CDataMapping','scaled');
    curr_axis = gca;
    curr_axis.YDir = 'normal';
    colorbar;
    axis tight
    xlabel('Residue');
    ylabel('Residue');
    title('Disorder matrix');
    axis equal
end

% generate disorder profile
dprofile = zeros(1,nres);
for r = 1:nres
    check = nres - sum(disorder(r,1:nres));
    if check >= options.mindomain
        dprofile(r) = 1;
    end
end

% determine the residue ranges for the folded domains (if any)
ordered = false;
ndomains = 0;
ranges = zeros(1000,2);

for r = 1:nres
    if ~ordered && dprofile(r)
        ndomains = ndomains + 1;
        ranges(ndomains,1) = r;
    end
    if ordered && ~dprofile(r)
        ranges(ndomains,2) = r-1;
    end
    ordered = dprofile(r);
end

if ndomains == 0
    domains.ranges = [];
    pae0 = sum(sum(pae));
    % mean pae in intrinsically disordered protein, disregarding the zeros between a residue
    % and itself
    domains.disorder = pae0/(nres*(nres-1));
    return
end

if ranges(ndomains,2) == 0
    ranges(ndomains,2) = nres;
end
domains.ranges = ranges(1:ndomains,:);
domains.disorder = zeros(ndomains);

for d1 = 1:ndomains
    range1 = ranges(d1,:);
    pae0 = sum(sum(pae(range1(1):range1(2),range1(1):range1(2))));
    % mean pae in folded domain, disregarding the zeros between a residue
    % and itself
    domains.disorder(d1,d1) = pae0/((range1(2)-range1(1)+1)*(range1(2)-range1(1)));
    for d2 = d1+1:ndomains
        range2 = ranges(d2,:);
        pae1 = sum(sum(pae(range1(1):range1(2),range2(1):range2(2))));
        pae2 = sum(sum(pae(range2(1):range2(2),range1(1):range1(2))));
        domains.disorder(d1,d2) = pae1/((range1(2)-range1(1)+1)*(range2(2)-range2(1)+1));
        domains.disorder(d2,d1) = pae2/((range1(2)-range1(1)+1)*(range2(2)-range2(1)+1));
    end
end

