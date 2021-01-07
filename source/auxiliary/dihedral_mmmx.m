function ang = dihedral_mmmx(c1,c2,c3,c4)
% function ang = DIHEDRAL_MMMX(c1,c2,c3,c4)
%
% Compute dihedral angle for four points with coordinates c1,c2,c3,c4
%
% ang   dihedral angle in units of radians
%
% (c) 2015 Gunnar Jeschke
%

vec1=c1-c2;
vec2=c4-c3;
vec3=c2-c3;
vec3=vec3/norm(vec3);
vec1=vec1-sum(vec1.*vec3)*vec3;
vec2=vec2-sum(vec2.*vec3)*vec3;
vec1=vec1/norm(vec1);
vec2=vec2/norm(vec2);
vec4 = [vec1(2)*vec2(3)-vec1(3)*vec2(2),vec1(3)*vec2(1)-vec1(1)*vec2(3),vec1(1)*vec2(2)-vec1(2)*vec2(1)];
if norm(vec4)<1e-6
    sign0=1;
else
	vec4=vec4/norm(vec4);
    sign0=-sum(vec4.*vec3);
end
ang=sign0*acos(sum(vec1.*vec2));
