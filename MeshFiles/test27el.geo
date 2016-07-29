Point(1) = {0, 0, 0, 1.0};
Point(2) = {0, 1, 0, 1.0};
Point(3) = {1, 0, 0, 1.0};
Point(4) = {1, 1, 0, 1.0};
Point(5) = {0, 1, 1, 1.0};
Point(6) = {1, 0, 1, 1.0};
Point(7) = {1, 1, 1, 1.0};
Point(8) = {0, 0, 1, 1.0};
Line(1) = {5, 2};
Line(2) = {5, 8};
Line(3) = {8, 1};
Line(4) = {1, 2};
Line(5) = {8, 6};
Line(6) = {5, 7};
Line(7) = {7, 6};
Line(8) = {6, 3};
Line(9) = {3, 1};
Line(10) = {3, 4};
Line(11) = {4, 7};
Line(12) = {2, 4};
Line Loop(13) = {10, -12, -4, -9};
Plane Surface(14) = {13};
Line Loop(15) = {12, 11, -6, 1};
Plane Surface(16) = {15};
Line Loop(17) = {6, 7, -5, -2};
Plane Surface(18) = {17};
Line Loop(19) = {5, 8, 9, -3};
Plane Surface(20) = {19};
Line Loop(21) = {2, 3, 4, -1};
Plane Surface(22) = {21};
Line Loop(23) = {11, 7, 8, 10};
Plane Surface(24) = {23};
Surface Loop(25) = {20, 18, 16, 14, 24, 22};
Volume(26) = {25};

Transfinite Line "*" = 5;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";

Physical Surface("top") = {18};
Physical Surface("bottom") = {14};
Physical Surface("sides") = {20, 24, 16, 22};
Physical Volume("cube") = {26};
