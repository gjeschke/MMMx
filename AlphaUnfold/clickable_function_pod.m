function clickable_function_pod(options)

if ~exist("options",'var') || ~isfield(options,'diffmode') || isempty(options.diffmode)
    options.diffmode = false;
end

if ~isfield(options,'clear') || isempty(options.clear)
    options.clear = true;
end

if ~isfield(options,'filtered') || isempty(options.filtered)
    options.filtered = 0;
end

if ~isfield(options,'threshold') || isempty(options.threshold)
    options.threshold = 30;
end


[filename, pathname] = uigetfile('function_*.csv', 'Select file');
if isequal(filename,0)
    disp('User selected Cancel');
else
    fullpath = fullfile(pathname, filename);
end

data = get_csv(fullpath,'data.cell');
[GO_ids,~] = size(data);

top_ten = -10*ones(10,3);
top_ten_GO = cell(10,3);
bottom_ten = 10*ones(10,3);
bottom_ten_GO = cell(10,3);

n_GO = 0;
n_stat = 0;

if options.clear
    figure; hold on
end

for p = 1:GO_ids
    GO_id = data{p,1};
    taxonomy = data{p,2};
    x = str2double(data{p,3});
    y = str2double(data{p,4});
    z = str2double(data{p,5});
    n = str2double(data{p,9});
    if options.diffmode
        xyz = [x,y,-z];
    else
        xyz = [x,y,1 - z];
    end
    stdx = str2double(data{p,6});
    stdy = str2double(data{p,7});
    stdz = str2double(data{p,8});
    stddev = [stdx,stdy,stdz];
    if options.filtered > 0 && max(stddev) > options.filtered % skip entries with broad distribution 
        continue
    end
    if n > options.threshold
        n_stat = n_stat + 1;
    end
    for k = 1:3
        if xyz(k) > min(top_ten(:,k))
            [top_ten(:,k),indices] = sort(top_ten(:,k),'descend');
            top_ten_GO(:,k) = top_ten_GO(indices,k);
            top_ten(10,k) = xyz(k);
            top_ten_GO{10,k} = GO_id;
        end
        if xyz(k) < max(bottom_ten(:,k))
            [bottom_ten(:,k),indices] = sort(bottom_ten(:,k),'ascen');
            bottom_ten_GO(:,k) = bottom_ten_GO(indices,k);
            bottom_ten(10,k) = xyz(k);
            bottom_ten_GO{10,k} = GO_id;
        end
    end
    plot_pod_ellipsoid(xyz,stddev,GO_id,taxonomy,n,options.filtered);
    n_GO = n_GO + 1;
end

fprintf(1,'%i gene ontology identifiers visualized\n',n_GO);
fprintf(1,'%i of these have at least %i proteins assigned\n',n_stat,options.threshold);

axis equal
grid on;
if ~options.diffmode
    axis([0,1,0,1,0,1]);
end

view(35,30);
lighting gouraud 
material shiny

if options.diffmode
    xlabel('\Delta f_{IDR}');
    ylabel('\Delta f_{fuzzy}');
    zlabel('\Delta f_{residual}');
else
    xlabel('f_{IDR}');
    ylabel('f_{fuzzy}');
    zlabel('1 - f_{residual}');
end

for k = 1:3
    [top_ten(:,k),indices] = sort(top_ten(:,k),'descend');
    top_ten_GO(:,k) = top_ten_GO(indices,k);
end
for k = 1:3
    [bottom_ten(:,k),indices] = sort(bottom_ten(:,k),'ascend');
    bottom_ten_GO(:,k) = bottom_ten_GO(indices,k);
end
