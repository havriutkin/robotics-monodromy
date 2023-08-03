-*

    This file contains all testing function.
    Any new functionality should be implemented as separate function in this file, 
    and then be called in main.m2 file.

*-


-- Load numerical packages
loadPackage "NumericalAlgebraicGeometry";
loadPackage "MonodromySolver";

-- Load others files
load "kinematics.m2";
load "parametrisation.m2";
load "monodromy.m2";


generateMonodromyGroup = () -> (
    -- Set up
    errorDepth = 0;
    setRandomSeed 0     -- Fix random seed

    -- Use KUKA KR-15/2 DH parameters
    dof := 6;
    alpha := {pi/2, 0.0, pi/2, -pi/2, pi/2, 0.0};
    r := {0.3, 0.65, 0.155, 0.0, 0.0, 0.0};
    d := {0.675, 0.0, 0.0, 0.0, 0.0, 0.140};
    theta := {0.0, -pi/2, 0.0, 0.0, 0.0, 0.0};

    -- Generate random thetas
    Dh=generateRandomDhParameters(6);
    DH := {alpha, r, d, Dh#3};

    -- Construct system of inverse kinematics equations
    G = constructEquations(DH#0, DH#1, DH#2);
    Sys=polySystem(G);

    -- Generate seed for monodromy
    Seed = findSeed(DH#0, DH#1, DH#2, DH#3);

    -- Find rows to elimenate in order to make system squared
    findBadRows(sub(jacobian Sys, matrix Seed#1), 50);
    G = G_{5,7,8,9,10,11}|take (G,-6);  -- TODO: Automate this part
    SquareSys = polySystem(G);

    -- Perform monodromy solve
    V = first monodromySolve(SquareSys, Seed#0, {Seed#1}, NumberOfEdges=>5, Verbose=>true, FilterCondition=>filterFunction);
    saturateEdges V.Graph

    -- Extract solutions from graph
    solutions = points V.PartialSols

    -- Generate monodromy group using graph
    G = monodromyGroup(V.Graph, FileName => "6RFinal.gp")
);

