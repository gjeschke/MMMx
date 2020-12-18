entity = get_pdb('1nti');

label = '(A)65';

taur = 2e-9;
taui = 2e-10;
td = 10e-3;
R2_dia = 12.6;
larmor = 750;

site_list(86).chain = 'A';
site_list(86).residue = 86;
for k = 1:86
    site_list(k).chain = 'A';
    site_list(k).residue = k;
end

[pre_list,exceptions] = pre_of_entity(entity,site_list,label,td,taur,taui,R2_dia,larmor);

figure(3); clf; hold on;
for k = 1:length(pre_list)
    plot(pre_list(k).residue,pre_list(k).pre,'k.','MarkerSize',14);
end
title(label);
set(gca,'FontSize',14);
xlabel('Residue');
ylabel('I_{para}/I_{dia}');
axis([1,86,0,1.05]);

return
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
