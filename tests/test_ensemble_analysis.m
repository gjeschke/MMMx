clear options
entity = get_pdb('2ad9');


[backbones,pop,exceptions] = get_backbones_ensemble(entity);

if ~isempty(exceptions{1})
    keyboard
end

[backbones2,pop2,exceptions] = get_backbones_ensemble('2adx');

if ~isempty(exceptions{1})
    keyboard
end

[common,dl1,dl2] = comp_struct(backbones,backbones2,2,0,1e-5);
