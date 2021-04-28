function [coor,min_dist,cycles,max_shift] = clash_repair(coor,environ,fixed,options)
% [coor,min_dist,cycles,max_shift] = clash_repair(ecoor,environ,fixed)
%
% tries to deform the structure specified by coor, so that it does not
% clash strongly with the environment specified by environ, atoms with
% coordinates provided in array fixed are not moved
% the algorithm tries to minimize local deformation
%
% coor      coordinate array of the "soft" structure (n,3)
% environ   coordinate array of the "hard" environment (m,3)
% fixed     coordinate array of atoms in the "soft" structure that should
%           not move (f,3)
% options   structure with fields
%           .clash_thr  clash threshold, defaults to 1.4 Å
%           .max_shift  maximum allowed shift of a soft structure atom,
%                       defaults to 2 Å
%
% coor      update coordinate array of the soft structure
% min_dist  minimum distance for any pair of "soft" structure and
%           environment atoms after deformation, is 1.4 Å (minimum H bond 
%           length) on success, remove H atoms from coor if you really need
%           close packing or reduce clash_thr
% cycles    iteration cycles
% max_shift maximum shift of any atom in the "soft" structure
%
% G. Jeschke, 3.1.2018

maxcycles = 100;
add_safety = 5; % additional safety margin beyond shift and clash threshold

if ~exist('options','var') || ~isfield(options,'clash_thr') || isempty(options.clash_thr)
    clash_thr = 1.4;
else
    clash_thr = options.clash_thr;
end

if ~exist('options','var') || ~isfield(options,'max_shift')
    options.max_shift = 2.5;
end

% only environment atoms closer to any soft structure atom than this 
% safety margin at the beginning are considered in iterations
safety_margin = add_safety + options.max_shift + clash_thr;

[n,~] = size(coor);
[f,~] = size(fixed);
cycles = 0;
max_shift = 0;

coor0 = coor;

% determine closest approach and the corresponding atom pair
pair_dist = get_all_pair_dist(coor,environ);
mdv = min(pair_dist);
min_dist = min(mdv);


% check whether the structure is already non-clashing
if min_dist >= clash_thr
    return
end

[~,km] = find(pair_dist < safety_margin);
environ = environ(unique(km),:);
pair_dist = get_all_pair_dist(coor,environ);
[mdv,mdp] = min(pair_dist);
[min_dist,mdp2] = min(mdv);
mdp1 = mdp(mdp2);
mdvec = zeros(1,maxcycles+1);
mdvec(1) = min_dist;

while min_dist < clash_thr && cycles < maxcycles
    cycles = cycles + 1;
    if min_dist < 0.5*clash_thr % structures are likely to overlap
        [kn,km] = find(pair_dist < clash_thr);
        soft = coor(unique(kn),:);
        [ms,~] = size(soft);
        if ms > 1
            soft = mean(coor(unique(kn),:));
        end
        if length(km) > 1
            hard = mean(environ(unique(km),:));
        else
            hard =environ(unique(km),:);
        end
    else
        hard = environ(mdp2,:); % clashing environment atom
        soft = coor(mdp1,:); % clashing soft structure atom  
    end
    shift = clash_thr-min_dist; % required shift of the clashing atom 
    shiftvec = soft - hard; % vector along which the atoms are to be shifted
    shiftvec = shift*shiftvec/norm(shiftvec); % actual shift vector for the clashing atom
    if f > 0 % there are fixed atoms in the soft structure
        diff = repmat(soft,f,1) - fixed; % vectors between fixed atoms and clashing soft atom
        dist = sqrt(sum(diff.^2,2)); % length of these vectors
        decay = min(dist); % minimum distance of moving atom to a fixed atom
        if decay < 3*shift % if the moving atom is too close to a fixed atom, give up
            return
        end
    else
        decay = 1e10; % if there is no fixed atom, the whole soft structure is translated
    end
    diff = coor - repmat(soft,n,1); % vectors between soft structure atoms an moving atom
    shift_frac = 1-sqrt(sum(diff.^2,2))/decay; % fraction of shift vector that should be applied
    shift_frac(shift_frac < 0) = 0; % atoms further away than the decay distance should not be shifted
    coor = coor + shift_frac*shiftvec;
    % determine closest approach and the corresponding atom pair
    pair_dist = get_all_pair_dist(coor,environ);
    [mdv,mdp] = min(pair_dist);
    [min_dist,mdp2] = min(mdv);
    mdvec(cycles+1) = min_dist;
    if isnan(min_dist)
        min_dist = 0;
        max_shift = 100;
        return
    end
    mdp1 = mdp(mdp2);
    diff = coor - coor0; % coordinate changes as vectors
    dist = sqrt(sum(diff.^2,2)); % length of these vectors
    max_shift = max(dist);
    if max_shift > options.max_shift
        return
    end
end

% if min_dist < clash_thr
%     figure(1); clf;
%     plot(mdvec(1:cycles+1));
%     disp('Aber hallo!');
% end


