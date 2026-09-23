function plot_taxonomy(tree,leaf_nodes,taxonomy,feature,options)
% plot_taxonomy(tree,feature)
%
% Generates a circular plot of a taxonomy tree color-coded by a selected
% feature
%
% Input:
% tree          array of struct (one struct per node) with at least fields
%               .name           name of the node
%               .identifier     taxonomy identifier
%               .parent         parent in the tree
%               .(feature)      value of the selected feature
% leaf_nodes    taxonomy identifiers of the leaf nodes
% taxonomy      taxonomy index table for nodes
% feature       name of the feature for color-coding
% options       struct with options
%               .colors     color map, defaults to 150 shades in parula
%               .range      value range for color coding [min_v,max_v]
%               .coverage   optional value for minimum coverage, default
%                           0.95
%
% G. Jeschke, 2026

smallest_proteome = 62; % Vidania strain VFMALBOS, https://doi.org/10.1038/s41467-026-69238-x
minimum_coverage = 0.95;
excluded = 0;

if exist('options','var') && isfield(options,'coverage')
    minimum_coverage = options.coverage;
end

% make color map
if ~exist('options','var') || ~isfield(options,'colors')
    options.colors = parula(150);
end
[clength,~] = size(options.colors);

if ~isfield(options,'mark_homo_sapiens')
    options.mark_homo_sapiens = false;
end

% format feature string if it contains subscript
parts = split(feature,'_');
if length(parts) > 1
    feature_string = sprintf('%s_{%s}',parts{1},parts{2});
else
    feature_string = feature;
end

leaf_nodes = leaf_nodes(leaf_nodes ~= 1);
if ~isfield(options,'range') || isempty(options.range)
    % determine the range of values for the selected feature
    min_value = 1e12;
    max_value = -1e12;
    for n = 1:length(leaf_nodes)
        node = leaf_nodes(n);
        assigned = tree(node).proteins/tree(node).proteome_size;
        % skip proteomes with low AlphaFold coverage or with double coverage
        if assigned < minimum_coverage || assigned > 1
            continue
        end
        if tree(node).proteins < smallest_proteome
            continue
        end
        value = tree(node).(feature);
        if value < min_value
            min_value = value;
        end
        if value > max_value
            max_value = value;
        end
    end
else
    min_value = options.range(1);
    max_value = options.range(2);
end
fprintf(1,'Plot range: [%.3f,%.3f]\n',min_value,max_value);

proteins = zeros(length(leaf_nodes),1);
coverage = zeros(length(leaf_nodes),1);
values = zeros(length(leaf_nodes),1);
for n = 1:length(leaf_nodes)
    node = leaf_nodes(n);
    proteins(n) = tree(node).proteins;
    coverage(n) = tree(node).proteins/tree(node).proteome_size;
    if proteins(n) > 0 && (coverage(n) > 1 || coverage(n) < minimum_coverage)
        excluded = excluded + 1;
    end
    values(n) = tree(node).(feature);
end
relevant_nodes = leaf_nodes(coverage >= minimum_coverage & coverage <= 1);

values = values(coverage >= minimum_coverage & coverage <= 1);
covered = length(values);
[~,idx] = sort(values);
relevant_nodes = relevant_nodes(idx);
dphi = 2*pi/length(relevant_nodes);
node_list = zeros(length(taxonomy),6);
nodes = 0;
tree_depth = 0;
% determine maximum tree depth and make list of all (parent) nodes
for n = 1:length(relevant_nodes)
    node = relevant_nodes(n);
    if tree(node).identifier == 1
        continue
    end
    nodes = nodes + 1;
    phi = (n-1)*dphi;
    node_list(nodes,1) = node; % node identifier
    node_list(nodes,2) = tree(node).(feature); % associated value
    node_list(nodes,3) = 1; % multiplicity
    node_list(nodes,4) = phi; % angle in cicular plot
    node_list(nodes,5) = tree(node).depth; % depth
    if tree(node).depth > tree_depth
        tree_depth = tree(node).depth;
    end
    pnode = tree(node).parent;
    [isa,pnode] = ismember(pnode,taxonomy);
    if ~isa
        continue
    end
    node_list(nodes,6) = pnode; % parent node
    while ~isempty(tree(node).parent) && tree(node).parent ~= 0
        pnode = tree(node).parent;
        [isa,pnode] = ismember(pnode,taxonomy);
        if ~isa
            break
        end
        node = pnode; % parent node
        [isa,nl] = ismember(node,node_list(1:nodes,1));
        if isa
            node_list(nl,2) = (node_list(nl,3)*node_list(nl,2) + tree(node).(feature))/(node_list(nl,3)+1); % mean value
            node_list(nl,4) = (node_list(nl,3)*node_list(nl,4) + phi)/(node_list(nl,3)+1); % mean angle
            node_list(nl,3) = node_list(nl,3) + 1; % multiplicity
        else
            if ~isempty(tree(node).(feature))
                nodes = nodes + 1;
                node_list(nodes,1) = node; % node identifier
                node_list(nodes,2) = tree(node).(feature); % associated value
                node_list(nodes,4) = phi; % angle
                node_list(nodes,3) = 1; % multiplicity
                node_list(nodes,5) = tree(node).depth; % depth
            end
        end
        pnode = tree(node).parent;
        if isempty(pnode)
            break
        end
        [isa,pnode] = ismember(pnode,taxonomy);
        if ~isa
            break
        end
        node_list(nodes,6) = pnode; % parent
    end
end
node_list = node_list(1:nodes,:); 

figure; hold on;

outer_radius = 1; % 90*tree_depth; % length(leaf_nodes);
msize = 12;
% % plot the tree
% for n = 1:nodes
%     phi = node_list(n,4);
%     depth = node_list(n,5);
%     radius = depth*outer_radius/tree_depth;
%     x = radius*cos(phi);
%     y = radius*sin(phi);
%     if node_list(n,6) > 0
%         [isa,node] = ismember(node_list(n,6),node_list(:,1));
%         if isa
%             if node_list(node,4) ~= phi
%                 keyboard
%             end
%             phi = node_list(node,4);
%             depth = node_list(node,5);
%             radius = depth*outer_radius/tree_depth;
%             xp = radius*cos(phi);
%             yp = radius*sin(phi);
%             plot([x,xp],[y,yp],'Color',[0.75,0.75,0.75]);
%         else
%             continue
%         end
%     else
%         continue
%     end
% end
% plot the nodes
for n = 1:nodes
    node = node_list(n,1);
    if tree(node).proteins < smallest_proteome
        continue
    end
    phi = node_list(n,4);
    depth = node_list(n,5);
    radius = depth*outer_radius/tree_depth;
    x = radius*cos(phi);
    y = radius*sin(phi);
    if node_list(n,6) == 0
        color = [0,0,0];
    else
        color_idx = 1 + round((clength-1)*(node_list(n,2) - min_value)/(max_value-min_value));
        color = options.colors(color_idx,:);
    end
    if options.mark_homo_sapiens
        if contains(tree(node).name,'Homo sapiens')
            color = [0.7,0.1,0.2];
        end
    end
    obj = plot(x,y,'.','Color',color,'MarkerSize',msize);
    obj.UserData.name = tree(node).name;
    obj.UserData.TaxonId = tree(node).identifier;
    obj.UserData.value = tree(node).(feature);
    obj.UserData.proteome = tree(node).proteome;
    if isfield(tree(node),'proteins')
        obj.UserData.proteins = tree(node).proteins;
    else
        obj.UserData.proteins = NaN;
    end
    obj.ButtonDownFcn = @taxonomy_node_clicked;
    obj.DataTipTemplate.DataTipRows(1).Label = tree(node).name;
    obj.DataTipTemplate.DataTipRows(1).Value = tree(node).proteins;
    obj.DataTipTemplate.DataTipRows(2).Label = feature_string;
    obj.DataTipTemplate.DataTipRows(2).Value = tree(node).(feature);
end

axis equal
axis off
caxis([min_value max_value]); 
c = colorbar;
if strcmpi(feature_string,'RMSF')
    feature_string = [feature_string '(Å)'];
end
ylabel(c, feature_string, 'FontSize', 12);

fprintf(1,'Out of %i proteomes, %i are covered at minimum fraction of %.3f\n',length(leaf_nodes),covered,minimum_coverage);
fprintf(1,'%i proteomes were excluded for being outside the coverage range\n',excluded);