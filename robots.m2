loadPackage "NumericalAlgebraicGeometry";
loadPackage "MonodromySolver";

R = CC[tx, ty, tz, ax, ay, az, cx, cy, cz, sx, sy, sz][ct1, ct2, ct3, ct4, ct5, ct6, st1, st2, st3, st4, st5, st6];

-- Function to generate random DH parameters
generateRandomDhParameters = dof -> (
    alpha := {};
    r := {};
    d := {};
    theta := {};

    -- Generate random values
    for i from 1 to dof do (
        alpha = append(alpha, random(-pi/2, pi/2));
        r = append(r, random(-10, 10));
        d = append(d, random(0, 10));
        theta = append(theta, random(-pi, +pi));
    );

    return {alpha, r, d, theta};
);

-- Forward kinematics function
forwardKinematics = (alpha, r, d, theta) -> (
  -- Define the transformation matrix for one joint
  T := i -> matrix{{cos(theta#i), -sin(theta#i)*cos(alpha#i),  sin(theta#i)*sin(alpha#i), r#i*cos(theta#i)},
                   {sin(theta#i),  cos(theta#i)*cos(alpha#i), -cos(theta#i)*sin(alpha#i), r#i*sin(theta#i)},
                   {0,             sin(alpha#i),               cos(alpha#i),              d#i},
                   {0,             0,                          0,                         1}};

  -- Compute the forward kinematics transformation matrix
  Tf := product(#d, i -> T(i));

  -- Return the result matrix
  return Tf;
);

-- Construct inverse kinematics equations
constructEquations = (A, alpha, r, d) -> (
    ctheta := {ct1, ct2, ct3, ct4, ct5, ct6};
    stheta := {st1, st2, st3, st4, st5, st6};
    
    -- Construct SE(3) with unknowns
    T := i -> matrix{{ctheta#i, -stheta#i*cos(alpha#i),  stheta#i*sin(alpha#i), r#i*ctheta#i},
                    {stheta#i,  ctheta#i*cos(alpha#i), -ctheta#i*sin(alpha#i), r#i*stheta#i},
                    {0,         sin(alpha#i),           cos(alpha#i),          d#i},
                    {0,         0,                      0,                     1}};
    
    Tf := product(#d, i -> T(i));

    -- Set correspond elements of A and T to be equal
    myEquations := {};
    for i from 0 to 2 do(
        for j from 0 to 3 do( 
            myEquations = append(myEquations, Tf_(i,j) - A_(i,j));
        );
    );
    
    --for i from 0 to #alpha-1 do(
    --    myEquations = append(myEquations, (ctheta#i)^2 + (stheta#i)^2 -1);
    --);
    
    return myEquations;
);

-- Parametrise SE(3) using c_i and s_i
constructParametriseSE = () -> (
    -- SE(3) along X
    transX := matrix{{1, 0,  0,   ax},
                    {0, cx, -sx, 0},
                    {0, sx, cx,  0},
                    {0, 0,  0,   1}};
    
    -- SE(3) along Y
    transY := matrix{{cy,  0, sy, 0},
                    {0,   1, 0,  ay},
                    {-sy, 0, cy, 0},
                    {0,   0, 0,  1}};
                
    -- SE(3) along Z            
    transZ := matrix{{cz, -sz, 0, 0},
                    {sz, cz,  0, 0},
                    {0,  0,   1, az},
                    {0,  0,   0, 1}};
    
    result := transX * transY * transZ;
    return result;
);

-- Parametrise SE(3) using rational parametrisation of cos and sin
constructRationalParametriseSE = () -> (
    -- SE(3) along X
    transX := matrix{{1, 0, 0, ax}, 
                    {0, (1-tx^2)/(1+tx^2), (-2*tx)/(1+tx^2)}, 
                    {0, (2*tx)/(1+tx^2), (1-tx^2)/(1+tx^2), 0}, 
                    {0, 0, 0, 1}};
    
    -- SE(3) along Y
    transY := matrix{{(1-ty^2)/(1+ty^2), 0, (2*ty)/(1+ty^2), 0}, 
                    {0, 1, 0, ay},
                    {(-2*ty)/(1+ty^2), 0, (1-ty^2)/(1+ty^2), 0}, 
                    {0, 0, 0, 1}};
    
    -- SE(3) along Z
    transZ := matrix{{(1-tz^2)/(1+tz^2), (-2*tz)/(1+tz^2), 0, 0}, 
                    {(2*tz)/(1+tz^2), (1-tz^2)/(1+tz^2), 0, 0},
                    {0, 0, 1, az}, 
                    {0, 0, 0, 1}};
                
    result := transX * transY * transZ;
    return result;
);

-- Extracts thetax
extractXRotation = (A) -> ( 
    B:=matrix{{0,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,0}};
    y:=vector{0_RR,1_RR,0_RR,0_RR};
    Projx:=(B*(A*(y)));
    cx:= (Projx_0*y_0)+(Projx_1*y_1)+(Projx_2*y_2)+(Projx_3*y_3);
    thetaX = acos(cx);
    return thetaX;
);

-- Extracts thetay
extractYRotation = (A) -> ( 
    B:=matrix{{1,0,0,0}, {0,0,0,0}, {0,0,1,0}, {0,0,0,0}};
    z:=vector{0_RR,0_RR,1_RR,0_RR};
    Projy:=(B*(A*(z)));
    cy:= (Projy_0*z_0)+(Projy_1*z_1)+(Projy_2*z_2)+(Projy_3*z_3);
    thetaY = acos(cy);
    return thetaY;
);

-- Extracts thetaz
extractZRotation = (A) -> ( 
    B:=matrix{{1,0,0,0}, {0,1,0,0}, {0,0,0,0}, {0,0,0,0}};
    x:=vector{1_RR,0_RR,0_RR,0_RR};
    Projz:=(B*(A*(x)));
    cz:= (Projz_0*x_0)+(Projz_1*x_1)+(Projz_2*x_2)+(Projz_3*x_3);
    thetaZ = acos(cz);
    return thetaZ;
);

-- Extracts ax, ay, az
extractTranslation = (A) -> ( 
    v:=vector{0,0,0,1};
    S:=A*v;
    T:=vector{S_0, S_1, S_2};
   
    cx:=cos(extractXRotation(A));
    cy:=cos(extractYRotation(A));
    cz:=cos(extractZRotation(A));
    sx:=sin(extractXRotation(A));
    sy:=sin(extractYRotation(A));
    sz:=sin(extractZRotation(A));
    
    M:=matrix{{1, 0, sy}, {0, cy, -cy*sy}, {0, sx, cx*cy}};
    
    C:=solve(M, T);
    
    return entries C;
);

constructProblemSeed = (alpha, r, d, theta) -> (
    forwardSE := forwardKinematics(alpha, r, d, theta);
    thetaX := extractXRotation(forwardSE);
    thetaY := extractYRotation(forwardSE);
    thetaZ := extractZRotation(forwardSE);
    trans := extractTranslation(forwardSE);
    
    problemSeed := point{{trans#0, trans#1, trans#2, cos(thetaX), cos(thetaY), cos(thetaZ),
        sin(thetaX), sin(thetaY), sin(thetaZ)}};
    
    return problemSeed;
);

-- Creates seed for monodromy solver, returns list {seedProblem, seedSolution}
createMonodromySeed = (alpha, r, d, theta) -> (
    -- Create seed problem using forward kinematics
    problemSeed = constructProblemSeed(alpha, r, d, theta);
    
    -- Create seed solution
    solutionSeed := {};
    
    -- Add cosines of thetas
    for i from 0 to #theta - 1 do(
        solutionSeed = append(solutionSeed, cos(theta#i));
    );

    -- Add sines of thetas
    for i from 0 to #theta - 1 do(
        solutionSeed = append(solutionSeed, sin(theta#i));
    );

    -- Make seed solution a point
    solutionSeed = point{solutionSeed};

    seed := {problemSeed, solutionSeed};
    return seed;
);

testPlanar = () -> (
    -- Define DH parameters randomly
    dof := 2;
    alpha := {0.0, 0.0};
    r := {1.0, 1,0};
    d := {0.0, 0.0};
    theta := {pi/4, pi/4};

    dhParams := {alpha, r, d, theta};

    -- Compute the forward kinematics
    result := forwardKinematics(dhParams#0, dhParams#1, dhParams#2, dhParams#3);

    -- Print the result
    << "Result matrix: " << result << endl << endl;
)

main = () -> (
    -- Define dh parameters
    dof := 2;
    alpha = {0.0_CC, 0.0_CC};
    r = {1.0_CC, 1.0_CC};
    d = {0.0_CC, 0.0_CC};
    theta = {pi/2_CC, -pi/2_CC};   
    dhParams := {alpha, r, d, theta};
    
    -- Create parametrised end effector position
    A = constructParametriseSE();
    
    -- Construct equations and system
    myEquations = constructEquations(A, dhParams#0, dhParams#1, dhParams#2);
    mySystem = polySystem(myEquations);
    
    -- Create seed
    mySeed = createMonodromySeed(alpha, r, d, theta);
    
    
    -- Use monodromySolve
    V = monodromySolve(mySystem, mySeed#0, {mySeed#1}, Verbose=>true);
);

main()

