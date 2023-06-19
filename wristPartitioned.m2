-- Function finds normalised vectors of axies along which joints rotate.
--          returns list of vectors 
findRotationalAxies = (alpha, r, d, theta) -> (
    e = {vector{0_RR, 0_RR, 1_RR}};
    for i from 0 to #alpha-1 do(
        rotationMatr = matrix{{1, 0, 0},
                              {0, cos(alpha#i), -sin(alpha#i)},
                              {0, sin(alpha#i), cos(alpha#i)}};
        e = append(e, rotationMatr * e#i);
    );
    
    return e;
);

-- Function returns vector of angular velocities
findAngularVelocities = (alpha, r, d, theta) ->(
    axies = findRotationalAxies(alpha, r, d, theta);
    w = {};
    for i from 0 to #theta - 1 do(
        w = append(w, theta#i * axies#i);
    );
    return w;
);
