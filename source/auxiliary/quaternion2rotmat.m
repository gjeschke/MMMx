function R = quaternion2rotmat(q)
% R = QUATERNION2ROTMAT(q)
% generates rotation matrix R from normalized quaternion q

q1q = q(2)^2;
q2q = q(3)^2;
q3q = q(4)^2;
q01 = q(1)*q(2);
q02 = q(1)*q(3);
q03 = q(1)*q(4);
q12 = q(2)*q(3);
q13 = q(2)*q(4);
q23 = q(3)*q(4);
R = [1-2*q2q-2*q3q,2*(q12+q03),2*(q13-q02);...
    2*(q12-q03),1-2*q1q-2*q3q,2*(q23+q01);...
    2*(q13+q02),2*(q23-q01),1-2*q1q-2*q2q];
