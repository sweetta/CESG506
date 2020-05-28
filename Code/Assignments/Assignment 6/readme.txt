Tatsuhiko Sweet 05/27/2020
The following files were created for CESG 506 Assignment 6.

*Note: all files are IDENTICAL to Assignment 5 submission EXCEPT for the curved beam element class,
which is contained in 'curvedForceBeamClass.py'

Problem 6-3
    main5-3.py:
        "Study of the snap-though behavior of a shallow arch."
        Solves using arc-length method. A mesh refinement study was performed by
        running with different values of nEle, which represents the number of elements.
        Minor adjustments to arc-length initialization parameters may be needed for
        different meshes.

    Prob6-3_GammaVsDisp-nEle.png:
        Gamma vs Displacement curve plot output for problem 5-3. The middle node, and
        quarter node vertical displacements are plotted. The title contains the number
        of elements used for the mesh (from 4 to 64, at which point the result converge).

Problem 6-4
    main6-4-1.py:
        "Buckling of a beam." Computes the change in the det(Kt) as the axial load on
        a straight beam approaches the buckling load, Pcr. Part 1 of assignment problem 5-4.

    Prob6-4_GammaVsDetKt-StraightBeam.png:
        Plots the determinant of Kt of a straight beam (formulated using a shallow beam element
        with elevation function = 0) against gamma (gamma = 1 is when P = Euler Buckling Load)

    main6-4-2.py:
        Repeating problem 5-4-1 but with an imperfection using the elevation funcition

    Prob6-4_GammaVsDetKt_nEle.png:
        Plots the determinant of Kt of the imperfect shallow-curved beam. The title includes the
        number of elements used in the mesh. (Note, for meshes with high number of elements,
        the det(Kt) shot off to inf, thus some of the plots are empty).

    Prob6-4-2_GammaVsDisp_nEle.png:
        Gamma vs Displacement curve plot output for problem 5-4. The middle node, and
        quarter node vertical displacements are plotted. The title contains the number
        of elements used for the mesh.

curvedForceBeamclass.py:
    Python class for a curved beam that is used for both HW problems. Weak form of constitutive relation used
    for the axial portion to reduce/eliminate locking.

helperFunctions.py:
    A file full of helper functions including shape functions, Gaussian quadrature, and assembly.