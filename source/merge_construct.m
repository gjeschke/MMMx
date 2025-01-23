function entity = merge_construct(entities,parts)

entity.name = 'MMMx';
entity.populations = 1;
entity.selected = 0;
entity.index_array = zeros(200000,5,'uint16');
xyz = zeros(200000,3);
elements = zeros(1,200000,'uint8');
occupancies = ones(200000,1,'uint8');
atnum = 0;
for p = 1:length(parts)
    ctag = char(double('A')+p-1); % chain tag for this part in new entity
    entity.(ctag).selected = false;
    entity.(ctag).index = p;
    conf = parts(p).conformer;
    chain = parts(p).chain;
    range = parts(p).range(1):parts(p).range(2);
    c_entity = entities{parts(p).entity};
    for r = range % loop over all selected residues
        residue = sprintf('R%i',r);
        entity.(ctag).(residue).index = 1 + r - range(1);
        entity.(ctag).(residue).selected = false;
        entity.(ctag).(residue).selected_rotamers = 1;
        entity.(ctag).(residue).name = c_entity.(chain).(residue).name;
        entity.(ctag).(residue).locations = ' ';
        entity.(ctag).(residue).populations = 1;
        atoms = fieldnames(c_entity.(chain).(residue));
        for ka = 1:length(atoms)
            atom = atoms{ka};
            if isstrprop(atom(1),'upper') % these are atom fields
                this_atom = c_entity.(chain).(residue).(atom);
                at_index = this_atom.tab_indices(conf);
                atnum = atnum + 1;
                this_atom.tab_indices = atnum;
                xyz(atnum,:) = c_entity.xyz(at_index,:);
                elements(atnum) = c_entity.elements(at_index);
                entity.(ctag).(residue).(atom) = this_atom;
            end
        end
    end
end
entity.xyz = xyz(1:atnum,:);
entity.elements = elements(1:atnum);
entity.occupancies = occupancies(1:atnum);