whichone = 'transparency';
value = 10;

clear options
options.altlocs = false;
options.Bfactor = false;
options.element = false;
% id = cx_get_pdb('1a6m',options);

address = '(A)4.CG';
address = '(A)156.S';
address = '(A)154.FE';
address = '(B)';
% attribute = cx_get_atom(id,address,whichone),

attribute = cx_set_atom(id,address,whichone,value);

% attribute = cx_get_atom(id,address,whichone),
