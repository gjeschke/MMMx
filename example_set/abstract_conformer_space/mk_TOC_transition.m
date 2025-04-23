
entity1 = get_PED('PED00499','e002');
entity1 = inertia_frame(entity1);
entity2 = get_PED('PED00502','e002');
entity2 = inertia_frame(entity2);
options.visualization = 'GHSR1a_nP_P.mmm';
options.fname1 = 'PED499_nP.pdb';
options.fname2 = 'PED502_P.pdb';
options.figname = 'GHSR1a_transition_nP_P.pdf';

[clusters,D] = cluster_transition(entity1,entity2,options);
