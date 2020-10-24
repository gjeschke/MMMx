clear options
entity1 = get_pdb('2lzm');
options.altlocs = false;
options.Bfactor = false;
options.element = false;
% options.MMMcol = false;
[id,entity2] = cx_get_pdb('2lzm',options);
% entity1.xyz(317,2) = 0;
% entity1.A.R142.N.selected = 1;
[common,dl1,dl2] = comp_struct(entity1,entity2,2,0,1e-5);