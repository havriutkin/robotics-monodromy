-*

    This file contains all testing function.
    Any new functionality should be implemented as separate function in this file, 
    and then be called in main.m2 file.

*-


-- Load numerical packages
loadPackage "NumericalAlgebraicGeometry";
loadPackage "MonodromySolver";

-- Load others files
load 'kinematics.m2';
load 'parametrisation.m2';
load 'monodromy.m2';

testMonodromy = () -> (
    DH := generateRandomDhParameters(6);
    equations := constructEquations(DH#0, DH#1, DH#2);
    seed := findSeed(DH#0, DH#1, DH#2, DH#3);
    system := polySystem(equations);
    J := sub(jacobian system, matrix seed#1);
    badRows := findBadRows(J, 0.01);
    << badRows << endl;
);