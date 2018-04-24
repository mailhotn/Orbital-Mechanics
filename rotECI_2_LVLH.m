function [ Rot ] = rotECI_2_LVLH( O,i,w,th )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
R1 = [ cos(O)*cos(w) - sin(O)*cos(i)*sin(w), - cos(O)*sin(w) - sin(O)*cos(i)*cos(w),  sin(O)*sin(i);
       sin(O)*cos(w) + cos(O)*cos(i)*sin(w),   cos(O)*cos(i)*cos(w) - sin(O)*sin(w), -cos(O)*sin(i);
                              sin(i)*sin(w),                          cos(w)*sin(i),         cos(i)].';
R2 = [cos(th), sin(th), 0;
     -sin(th), cos(th), 0;
            0        0, 1];
Rot = R2*R1;
end

