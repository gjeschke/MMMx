whichone = 'color';
value = 'steelblue';

clear options
options.altlocs = false;
options.Bfactor = false;
options.element = false;
% id = cx_get_pdb('1a6m',options);
% id = cx_get_pdb('2ad9',options);
% id = cx_get_pdb('2lzm',options);

% address = '(A)57';
address = '(A)15';
% address = '(A)154.FE';
% address = '(A)156.S';
attribute = cx_get_residue(id,address,whichone),

attribute = cx_set_residue(id,address,whichone,value);

attribute = cx_get_residue(id,address,whichone),
