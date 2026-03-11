function [domain,error] = assign_domain(taxonomy_id)

persistent assignment
error = '';
domain = '';

if isnan(taxonomy_id) || taxonomy_id == 1
    error = 'Invalid taxonomy ID';
    return
end

if isempty(assignment)
    assignment = load('taxonomy_nodes.mat');
end
[diff,idx] = min(abs(assignment.nodes(:,1) - taxonomy_id));
if diff ~= 0
    domain = '';
    error = 'Invalid taxonomy identifier';
    return
end

depth = 100;
while assignment.nodes(idx,3) ~= 1 && depth > 0
    [diff,idx] = min(abs(assignment.nodes(:,1) - assignment.nodes(idx,2)));
    if diff ~= 0
        domain = '';
        error = 'Taxonomy hierarchy does not lead back to domain';
        return
    end
    depth = depth - 1;
end
    
if assignment.nodes(idx,3) ~= 1
    domain = '';
    error = 'Taxonomy hierarchy depth exceeded';
    return
end

[diff,idx2] = min(abs(assignment.decoder_ids - assignment.nodes(idx,1)));
if diff ~= 0
    domain = '';
    error = sprintf('Wrong taxonomy identifier %i at domain level',assignment.nodes(idx,1));
    return
end

domain = assignment.decoder(idx2).name;
