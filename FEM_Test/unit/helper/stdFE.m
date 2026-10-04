function fE = stdFE(name)
    % stdFE: standard 3D finite elements used in tests, cached within one run of `runTests` (cleared there).
    % Construction of 3D elements (e.g. P2, P3) takes seconds; the same elements are used by several test classes.

    % Name   | Element
    % ------ | --------------------------------------------------------------
    % P0     | constant, DoF at barycenter of element
    % P1     | linear, DoFs at vertices
    % P2     | quadratic, DoFs at vertices and midpoints of edges
    % P3     | cubic, DoFs at vertices, two points on each edge, barycenters of faces
    % DG1    | linear, DoFs at vertices, not shared
    % CR     | linear, DoFs at barycenters of faces
    % RT0    | lowest order Raviart-Thomas, DoFs are normal fluxes through faces
    % BDM1   | first order Brezzi-Douglas-Marini, DoFs are moments of normal component against barycentric coordinates of faces
    % NED1   | lowest order Nedelec (first kind), DoFs are tangential moments on edges
    % TP1    | linear trace on faces (reference face "D3FR"), DoFs at face vertices, not shared

    arguments (Input)
        name (1, 1) string {mustBeMember(name, ["P0", "P1", "P2", "P3", "DG1", "CR", "RT0", "BDM1", "NED1", "TP1"])};
    end
    arguments (Output)
        fE FE;
    end
    persistent cache
    if isempty(cache)
        cache = containers.Map();
    end
    if isKey(cache, char(name))
        fE = cache(char(name));
        return;
    end
    D3TElem = MshEnt("D3T").msh;
    d0 = [0; 0; 0];
    P2FS = "[1,x,y,z,x^2,y^2,z^2,x*y,y*z,z*x]";
    P3FS = "[1,x,y,z,x^2,y^2,z^2,x*y,y*z,z*x,x^3,y^3,z^3,x^2*y,x^2*z,y^2*x,y^2*z,z^2*x,z^2*y,x*y*z]";
    switch name
        case "P0"
            fE = FE("D3T", "1", NdDoF("D3", D3TElem, 3, [1/4; 1/4; 1/4], d0), "map", "affine");
        case "P1"
            fE = FE("D3T", "[1,x,y,z]", NdDoF("D3", D3TElem, 0, [], d0), "map", "affine");
        case "P2"
            fE = FE("D3T", P2FS, [NdDoF("D3", D3TElem, 0, [], d0), NdDoF("D3", D3TElem, 1, 1/2, d0)], "map", "affine");
        case "P3"
            fE = FE("D3T", P3FS, [NdDoF("D3", D3TElem, 0, [], d0), NdDoF("D3", D3TElem, 1, [1/3, 2/3], d0), ...
                NdDoF("D3", D3TElem, 2, [1/3; 1/3], d0)], "map", "affine");
        case "DG1"
            fE = FE("D3T", "[1,x,y,z]", NdDoF("D3", D3TElem, 0, [], d0, "share", false), "map", "affine");
        case "CR"
            fE = FE("D3T", "[1,x,y,z]", NdDoF("D3", D3TElem, 2, [1/3; 1/3], d0), "map", "affine");
        case "RT0"
            fE = FE("D3T", "[1,0,0; 0,1,0; 0,0,1; x,y,z].'", MoDoF("D3", D3TElem, 2, Fcn.cst(1), zeros(3), ...
                "coef", MshEnt("D3F").UNV, "orien", true, "GInt", GInt("D3F", 1)), "map", "piolaDiv");
        case "BDM1"
            fE = FE("D3T", FE.repFS("[1,x,y,z]", [3, 1]), MoDoF("D3", D3TElem, 2, ...
                [Fcn("D3R2", "1 - s - t"), Fcn("D3R2", "s"), Fcn("D3R2", "t")], zeros(3), ...
                "coef", MshEnt("D3F").UNV, "orien", true, "GInt", GInt("D3F", 2)), "map", "piolaDiv");
        case "NED1"
            fE = FE("D3T", "[1,0,0; 0,1,0; 0,0,1; 0,-z,y; z,0,-x; -y,x,0].'", MoDoF("D3", D3TElem, 1, Fcn.cst(1), zeros(3), ...
                "coef", MshEnt("D3L").UTV, "orien", true, "GInt", GInt("D3L", 1)), "map", "piolaCurl");
        case "TP1"
            fE = FE("D3FR", "[1,s,t]", NdDoF("D3R2", MshEnt("D3FR").msh, 2, [0, 1, 0; 0, 0, 1], [0; 0], "share", false));
    end
    cache(char(name)) = fE;
end
