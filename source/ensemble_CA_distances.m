function [M,S] = ensemble_CA_distances(entity,conformers)

if ~exist('conformers','var')
    conformers = 1:length(entity.populations);
end

[backbones,W] = get_backbones_ensemble(entity,'A');
coor0 = backbones.A.bb{1};
[na,~] = size(coor0);
C = length(conformers);
D = zeros(C,na,na);
W = W(conformers);
ind = 0;
for c = conformers
    ind = ind + 1;
    coor = backbones.A.bb{c};
    D(ind,:,:) = coor2dmat(coor); % distance matrix for first conformer
end

% Reshape weights to match dimensions of D
W_reshaped = reshape(W, [length(W), 1, 1]);

% Compute weighted mean M
weighted_sum = sum(W_reshaped .* D, 1);  % Sum along C dimension
M = weighted_sum / sum(W);
M = squeeze(M);  % Remove singleton dimension

% Compute weighted variance V
deviations = D - reshape(M, [1, size(M, 1), size(M, 2)]);  % Broadcast subtraction
squared_deviations = deviations .^ 2;
weighted_sq_sum = sum(W_reshaped .* squared_deviations, 1);
V = weighted_sq_sum / sum(W);
V = squeeze(V);
S = sqrt(V);