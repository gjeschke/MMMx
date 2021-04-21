function transmat = transrot2affine(trans,euler)
% Inverse transformation for affine2transrot
%
% G. Jeschke, 2017

transmat1 = affine('xyz',euler);
transmat2 = affine('translation',trans);
transmat = transmat2*transmat1;