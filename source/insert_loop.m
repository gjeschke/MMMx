function insert_loop_ensemble(core,loop_ensemble,N_anchor_chain,C_anchor_chain,loop_chain)

function insert_loop(in_name,out_name,loop_name,N_anchor_chain,C_anchor_chain,loop_chain,res1,resend)

if isempty(N_anchor_chain)
    anchor_chain = C_anchor_chain;
else
    anchor_chain = N_anchor_chain;
end
m_entity = get_pdb(loop_name);
entity1 = get_pdb(in_name);
for resnum = res1:resend
    resname = sprintf('R%i',resnum);
    entity1 = add_residue(entity1,anchor_chain,resnum,m_entity.(loop_chain).(resname),m_entity.xyz);
end
if ~isempty(N_anchor_chain) && ~isempty(C_anchor_chain)...
        && ~strcmp(N_anchor_chain,C_anchor_chain) % chain must be stitched
    residues = fieldnames(entity1.(C_anchor_chain));
    % the following does not produce a consistent
    % MMMx entity, but entity1 can be correctly
    % saved as PDB file
    for kres = 1:length(residues)
        resname = residues{kres};
        if resname(1) == 'R'
            entity1.(N_anchor_chain).(resname) = entity1.(C_anchor_chain).(resname);
        end
    end
    entity1 = rmfield(entity1,C_anchor_chain);
end
put_pdb(entity1,out_name);