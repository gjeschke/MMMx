entity = get_pdb('1nti');

taur = 2e-9;
taui = 2e-10;
td = 10e-3;
R2_dia = 12.6;
larmor = 750;

coor = zeros(86,3);
valid = true(86,1);
H_index = entity.A.R1.H1.tab_indices(1);
coor(1,:) = entity.xyz(H_index,:);
for res = 2:86
    resfield = sprintf('R%i',res);
    if isfield(entity.A.(resfield),'H')
        H_index = entity.A.(resfield).H.tab_indices(1);
        coor(res,:) = entity.xyz(H_index,:);
    else
        valid(res) = false;
    end
end

[positions,entity] = get_label(entity,'mtsl','positions','(A)86');
label.coor = positions{1};
[populations,entity,exceptions] = get_label(entity,'mtsl','populations','(A)86');
label.pop = populations{1};

all_pre = pre(coor,label,td,taur,taui,R2_dia,larmor);

figure(5); clf;
plot(all_pre,'r');
title('I86C');
