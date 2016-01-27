cl__1 = 1;
cl_vol1 = 0.004;
cl_vol2 = 0.05;
Point(1) = {0, 0, 0, cl_vol1};
Point(2) = {0.1, 0, 0, cl_vol1};
Point(3) = {0.1, 0.1, 0, cl_vol1};
Point(4) = {0, 0.1, 0, cl_vol1};
Point(5) = {0, 0, 0.1, cl_vol1};
Point(6) = {0.1, 0, 0.1, cl_vol1};
Point(7) = {0.1, 0.1, 0.1, cl_vol1};
Point(8) = {0, 0.1, 0.1, cl_vol1};

Point(9) = {0.1, 0, 0, cl_vol1};
Point(10) = {0.1, 0.1, 0, cl_vol1};
Point(11) = {0, 0.1, 0, cl_vol1};
Point(12) = {0, 0, 0.1, cl_vol1};
Point(13) = {0.1, 0, 0.1, cl_vol1};
Point(14) = {0.1, 0.1, 0.1, cl_vol1};
Point(15) = {0, 0.1, 0.1, cl_vol1};

Point(16) = {0.0, 0.2, 0.0, cl_vol2};
Point(17) = {0.2, 0.2, 0.0, cl_vol2};
Point(18) = {0.2, 0.0, 0.0, cl_vol2};
Point(19) = {0.0, 0.2, 0.2, cl_vol2};
Point(20) = {0.2, 0.2, 0.2, cl_vol2};
Point(21) = {0.0, 0.0, 0.2, cl_vol2};
Point(22) = {0.2, 0.0, 0.2, cl_vol2};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};
Line(13) = {10, 11};
Line(14) = {11, 15};
Line(15) = {15, 14};
Line(16) = {14, 10};
Line(17) = {9, 10};
Line(18) = {14, 13};
Line(19) = {13, 9};
Line(20) = {15, 12};
Line(21) = {12, 13};
Line(22) = {11, 16};
Line(23) = {16, 17};
Line(24) = {17, 18};
Line(25) = {18, 9};
Line(26) = {16, 19};
Line(27) = {19, 21};
Line(28) = {21, 12};
Line(29) = {21, 22};
Line(30) = {22, 18};
Line(31) = {19, 20};
Line(32) = {20, 17};
Line(33) = {20, 22};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Line Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};
Line Loop(3) = {9, -8, -12, 4};
Plane Surface(3) = {3};
Line Loop(4) = {9, 5, -10, -1};
Plane Surface(4) = {4};
Line Loop(5) = {10, 6, -11, -2};
Plane Surface(5) = {5};
Line Loop(6) = {12, -7, -11, 3};
Plane Surface(6) = {6};
Line Loop(7) = {14, 15, 16, 13};
Plane Surface(7) = {7};
Line Loop(8) = {-19, -18, 16, -17};
Plane Surface(8) = {8};
Line Loop(9) = {21, -18, -15, 20};
Plane Surface(9) = {9};
Line Loop(10) = {-28, 29, 30, 25, -19, -21};
Plane Surface(10) = {10};
Line Loop(11) = {-24, -23, -22, -13, -17, -25};
Plane Surface(11) = {11};
Line Loop(12) = {22, 26, 27, 28, -20, -14};
Plane Surface(12) = {12};
Line Loop(13) = {23, -32, -31, -26};
Plane Surface(13) = {13};
Line Loop(14) = {31, 33, -29, -27};
Plane Surface(14) = {14};
Line Loop(15) = {24, -30, -33, 32};
Plane Surface(15) = {15};

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};
Surface Loop(2) = {7, 8, 9, 10, 11, 12, 13, 14, 15};
Volume(2) = {2};

Physical Volume("ODLEW") = {1};
Physical Volume("FORMA") = {2};

Physical Surface("WB_4R_ODLEW") = {2, 5, 6};
Physical Surface("WB_3R_FORMA") = {13, 14, 15};

Periodic Surface 6 {12, -7, -11, 3} = 7 {14, 15, 16, 13};
Periodic Surface 5 {10, 6, -11, -2} = 8 {-19, -18, 16, -17};
Periodic Surface 2 {5, 6, 7, 8} =  9 {21, -18, -15, 20};




