/*
Gmsh .geo file to create mesh. Original file from:
https://openfoamwiki.net/index.php/2D_Mesh_Tutorial_using_GMSH

Still need to figure out a way to output boundary faces
*/

// Inputs
lenX = 2; //m
lenY = 2; //m

gridsize = 0.05; // lenX / 20;

// All numbering counterclockwise from bottom-left corner
Point(1) = {0, 0, 0, gridsize};
Point(2) = {lenX, 0, 0, gridsize};
Point(3) = {lenX, lenY, 0, gridsize};
Point(4) = {0, lenY, 0, gridsize};
Line(1) = {1, 2};               // bottom line
Line(2) = {2, 3};               // right line
Line(3) = {3, 4};               // top line
Line(4) = {4, 1};               // left line

//Physical Line(7) = {4};
//Physical Line(8) = {1};
//Physical Line(9) = {2};
//Physical Line(10) = {3};

Line Loop(5) = {1, 2, 3, 4};
// the order of lines in Line Loop is used again in surfaceVector[]
Plane Surface(6) = {5};

Transfinite Surface {6};
Recombine Surface {6};
