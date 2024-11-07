function AF3 = get_AF3_output(basname,seed)

fname = sprintf('fold_%s_%i.zip',basname,seed);
if exist(fname,'file')
    unzip(fname);
    mbasname = sprintf('fold_%s_%i_model_',basname,seed);
    for m = 1:5
        fname = sprintf('%s%i.cif',mbasname,m-1);
        entity = get_cif(fname);
        modnum = 5*(seed-1)+m;
        oname = sprintf('%s_m%i.pdb',basname,modnum);
        AF3.files{m} = oname;
        put_pdb(entity,oname);
        full_data = sprintf('fold_%s_%i_full_data_%i.json',basname,seed,m-1);
        entity.AF3_full_data = full_data;
        jsonStr = fileread(full_data);
        full_data = jsondecode(jsonStr);
        confidences = sprintf('fold_%s_%i_summary_confidences_%i.json',basname,seed,m-1);
        entity.AF3_confidences = confidences;
        matfname = sprintf('%s_m%i.mat',basname,modnum);
        save(matfname,'entity');
        jsonStr = fileread(confidences);
        confidences = jsondecode(jsonStr);
        AF3.atom_chain_ids{m} = full_data.atom_chain_ids;
        AF3.atom_plddts{m} = full_data.atom_plddts;
        AF3.contact_probs{m} = full_data.contact_probs;
        AF3.pae{m} = full_data.pae;
        AF3.token_chain_ids{m} = full_data.token_chain_ids;
        AF3.token_res_ids{m} = full_data.token_res_ids;
        AF3.chain_iptm{m} = confidences.chain_iptm;
        AF3.chain_pair_iptm{m} = confidences.chain_pair_iptm;
        AF3.chain_pair_pae_min{m} = confidences.chain_pair_pae_min;
        AF3.chain_ptm{m} = confidences.chain_ptm;
        AF3.fraction_disordered{m} = confidences.fraction_disordered;
        AF3.has_clash{m} = confidences.has_clash;
        AF3.iptm{m} = confidences.iptm;
        AF3.ptm{m} = confidences.ptm;
        AF3.num_recycles{m} = confidences.num_recycles;
        AF3.ranking_score{m} = confidences.ranking_score;
    end
end

