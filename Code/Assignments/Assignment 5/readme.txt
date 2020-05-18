Tatsuhiko Sweet 05/18/2020
The following files were created for CESG 506 Assignment 3.

Problem 5-3
    main5-3.py:
        "Study of the snap-though behavior of a shallow arch."
        Solves using arc-length method. A mesh refinement study was performed by
        running with different values of nEle, which represents the number of elements.
        Minor adjustments to arc-length initialization parameters may be needed for
        different meshes.

    Prob5-3_GammaVsDisp-nEle.png:
        Gamma vs Displacement curve plot output for problem 5-3. The middle node, and
        quarter node vertical displacements are plotted. The title contains the number
        of elements used for the mesh (from 4 to 64, at which point the result converge).

Problem 5-4
    main5-4-1.py:
        "Buckling of a beam." Computes the change in the det(Kt) as the axial load on
        a straight beam approaches the buckling load, Pcr. Part 1 of assignment problem 5-4.

    Prob5-4_GammaVsDetKt-StraightBeam.png:
        Plots the determinant of Kt of a straight beam (formulated using a shallow beam element
        with elevation function = 0) against gamma (gamma = 1 is when P = Euler Buckling Load)

    main5-4-2.py:
        Repeating problem 5-4-1 but with an imperfection using the elevation funcition

    Prob5-4_GammaVsDetKt_nEle.png:
        Plots the determinant of Kt of the imperfect shallow-curved beam. The title includes the
        number of elements used in the mesh. (Note, for meshes with high number of elements,
        the det(Kt) shot off to inf, thus some of the plots are empty).

    Prob5-4-2_GammaVsDisp_nEle.png:
        Gamma vs Displacement curve plot output for problem 5-4. The middle node, and
        quarter node vertical displacements are plotted. The title contains the number
        of elements used for the mesh.

curvedBeamclass.py:
    Python class for a curved beam that is used for both HW problems.

helperFunctions.py:
    A file full of helper functions including shape functions, Gaussian quadrature, and assembly.