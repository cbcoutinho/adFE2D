function [ area ] = areaTriangle( X1, X2, X3 )
%AREATRIANGLE Determines the area of a triangle from 3 points in 2D
%   A triangle is formed by 3 points. The area of that triangle is
%   determined using the approach defined here:
%       http://mathworld.wolfram.com/TriangleArea.html
%       http://www.mathopenref.com/coordtrianglearea.html

matrix = [X1(1) X1(2) 1;
    X2(1) X2(2) 1;
    X3(1) X3(2) 1];

area = 0.5*det(matrix);

end

