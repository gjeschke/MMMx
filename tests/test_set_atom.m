whichone = 'coor';
value = 10;

clear options
options.altlocs = false;
options.Bfactor = false;
options.element = false;
% id = cx_get_pdb('1a6m',options);

% address = '(A)4.CG';
% address = '(A)156.S';
address = '(A)154.FE';
% address = '(B)';
xyz = cx_get_atom(id,address,whichone);
fprintf(1,'Coordinates of %s : %8.3f,%8.3f,%8.3f ?\n',address,xyz);
xyz = xyz + [1,-2,-4];

attribute = cx_set_atom(id,address,whichone,xyz);

% attribute = cx_get_atom(id,address,whichone),
