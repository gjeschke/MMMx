function transmat=affine(mode,param1,param2,param3)
% function transmat=affine(mode,param1,param2,param3)
%
% Creates a 4x4 matrix for an affine coordinate transformation in 3D space
% see: M. Bender, M. Brill, Computergrafik, Hanser, M?nchen, 2. Aufl.,
%      2006, section 2.1
% use functions affine_trafo, affine_trafo_point, and affine_trafo_vector 
% for performing the actual transformations
% an empty matrix is returned when a non-existing mode is requisted
% function assumes that the correct number of correctly dimensioned
% parameters is provided by the client (caller)
%
% mode      transformation mode, case sensitive string
% param     parameters
%
% Implemented combinations of modes and parameters
% mode              explanation             parameters
% 'identity'        identity transformation (none)
% 'translation'     translation             coordinates [x,y,z] of new
%                                           origin
% 'Euler'           rotation                Euler angles [alpha,beta,gamma]
%                                           for rotations about z, y', and
%                                           z'' axes, in radians
%                                           WARNING: the objects are
%                                           rotated, not the frame, use
%                                           -alpha, -beta, -gamma for the
%                                           frame rotation more usual in
%                                           magnetic resonance
% 'xyz'             rotation                Euler angles for xyz convention
% 'rotx'            rotation about x axis   rotation angle phi (radians)
% 'roty'            rotation about y axis   rotation angle phi (radians)
% 'rotz'            rotation about z axis   rotation angle phi (radians)
% 'rotn'            rotation about an       rotation angle phi (radians)
%                   arbitrary axis          and vector representation
%                                           [nx,ny,nz] of the rotation axis
% 'invert'          inversion at origin     (none)
% 'reflectxy'       reflection at xy plane  (none) 
% 'reflectxz'       reflection at xz plane  (none) 
% 'reflectyz'       reflection at yz plane  (none) 
% 'scale'           scaling                 scaling factors [ax,ay,az] for
%                                           the three coordinates or
%                                           uniform scaling factor a for
%                                           all coordinates
% 'shear'           shearing                vector
%                                           [sxy,sxz,syx,syz,szx,szy] of
%                                           shearing coefficients
% 'screw'           rotation about an       rotation angle phi (radians),
%                   arbitrary axis          vector representation
%                   followed by trans-      [nx,ny,nz] of the rotation axis
%                   lation along the same   and shift dz along the axis
%                   axis
%
% G. Jeschke, 2009

transmat=eye(4); % initialize identity transformation

switch mode
    case 'identity' % needed to prevent otherwise clause for identity mode
    case 'translation'
        transmat(1:3,4)=param1';
    case 'Euler'
        ca=cos(param1(1)); sa=sin(-param1(1));
        cb=cos(param1(2)); sb=sin(-param1(2));
        cg=cos(param1(3)); sg=sin(-param1(3));
        transmat(1,1)=cg*cb*ca-sg*sa;
        transmat(1,2)=cg*cb*sa+sg*ca;
        transmat(1,3)=-cg*sb;
        transmat(2,1)=-sg*cb*ca-cg*sa;
        transmat(2,2)=-sg*cb*sa+cg*ca;
        transmat(2,3)=sg*sb;
        transmat(3,1)=sb*ca;
        transmat(3,2)=sb*sa;
        transmat(3,3)=cb;
    case 'xyz' % https://gamedev.stackexchange.com/questions/50963/how-to-extract-euler-angles-from-transformation-matrix
        cx=cos(param1(1)); sx=sin(param1(1));
        cy=cos(param1(2)); sy=sin(param1(2));
        cz=cos(param1(3)); sz=sin(param1(3));
        transmat(1,1)= cy*cz;
        transmat(1,2)= -cy*sz;
        transmat(1,3)= sy;
        transmat(2,1)= sx*sy*cz+cx*sz;
        transmat(2,2)= -sx*sy*sz+cx*cz;
        transmat(2,3)= -sx*cy;
        transmat(3,1)= -cx*sy*cz+sx*sz;
        transmat(3,2)= cx*sy*sz+sx*cz;
        transmat(3,3)= cx*cy;
    case 'rotx'
        cfi=cos(param1(1));
        sfi=sin(param1(1));
        transmat(2,2)=cfi;
        transmat(2,3)=-sfi;
        transmat(3,2)=sfi;
        transmat(3,3)=cfi;
    case 'roty'
        cfi=cos(param1(1));
        sfi=sin(param1(1));
        transmat(1,1)=cfi;
        transmat(1,3)=sfi;
        transmat(3,1)=-sfi;
        transmat(3,3)=cfi;
    case 'rotz'
        cfi=cos(param1(1));
        sfi=sin(param1(1));
        transmat(1,1)=cfi;
        transmat(1,2)=-sfi;
        transmat(2,1)=sfi;
        transmat(2,2)=cfi;
    case 'rotn'
        c=cos(param1(1));
        s=sin(param1(1));
        t=1-c;
        n=param2/norm(param2);
        nx=n(1);
        ny=n(2);
        nz=n(3);
        transmat(1,1)=t*nx^2+c;
        transmat(1,2)=t*nx*ny-s*nz;
        transmat(1,3)=t*nx*nz+s*ny;
        transmat(2,1)=t*nx*ny+s*nz;
        transmat(2,2)=t*ny^2+c;
        transmat(2,3)=t*ny*nz-s*nx;
        transmat(3,1)=t*nx*nz-s*ny;
        transmat(3,2)=t*ny*nz+s*nx;
        transmat(3,3)=t*nz^2+c;
    case 'invert'
        transmat(1,1)=-1;
        transmat(2,2)=-1;
        transmat(3,3)=-1;
    case 'reflectxy'
        transmat(3,3)=-1;
    case 'reflectxz'
        transmat(2,2)=-1;
    case 'reflectyz'
        transmat(1,1)=-1;
    case 'scale'
        if length(param1)==1, param1=[1,1,1]*param1; end
        transmat(1,1)=param1(1);
        transmat(2,2)=param1(2);
        transmat(3,3)=param1(3);
    case 'shear'
        transmat(2,1)=param1(1);
        transmat(3,1)=param1(2);
        transmat(1,2)=param1(3);
        transmat(3,2)=param1(4);
        transmat(1,3)=param1(5);
        transmat(2,3)=param1(6);
    case 'screw'
        c=cos(param1(1));
        s=sin(param1(1));
        t=1-c;
        n=param2/norm(param2);
        nx=n(1);
        ny=n(2);
        nz=n(3);
        transmat(1,1)=t*nx^2+c;
        transmat(1,2)=t*nx*ny-s*nz;
        transmat(1,3)=t*nx*nz+s*ny;
        transmat(2,1)=t*nx*ny+s*nz;
        transmat(2,2)=t*ny^2+c;
        transmat(2,3)=t*ny*nz-s*nx;
        transmat(3,1)=t*nx*nz-s*ny;
        transmat(3,2)=t*ny*nz+s*nx;
        transmat(3,3)=t*nz^2+c;
        transmat(1,4)=param3*param2(1);
        transmat(2,4)=param3*param2(2);
        transmat(3,4)=param3*param2(3);
    otherwise
        transmat=[];
end
        
       