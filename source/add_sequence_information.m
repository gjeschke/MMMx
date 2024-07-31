function entity = add_sequence_information(entity,chain)

if ~isfield(entity,chain)
    entity.(chain) = [];
    return
end
residues = fieldnames(entity.(chain));
sequence = char(zeros(1,length(residues)));
resnum = zeros(1,length(residues));
spoi = 0;
for kr = 1:length(residues)
    residue = residues{kr};
    if strcmp(residue(1),'R') % these are residue fields
        slc = tlc2slc(entity.(chain).(residue).name);
        if ~isempty(slc)
            spoi = spoi + 1;
            sequence(spoi) = slc;
            resnum(spoi) = str2double(residue(2:end));
        end
    end
end
entity.(chain).sequence = sequence(1:spoi);
entity.(chain).resnum = resnum(1:spoi);