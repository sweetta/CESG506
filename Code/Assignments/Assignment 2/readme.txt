Tatsuhiko Sweet 04/22/2020
The following files were created for CESG 506 Assignment 2.

Problem 2-1
    main.py:
        Solves the two-truss problem from last week's HW using using Newton-Raphson, Henkey strain,
        and displacement control. Utilizes trussClass but assembly is done manually (maybe next HW it'll be auto!).

    Prob2-1_GammaVsDisp.png:
        Gamma vs Displacement curve plot output for problem 2-1 part 3.

    Prob2-1_ContourPlot.png:
        u vs v (path of free node displacement) with internal force component overlay for problem 2-1 part 4.

Problem 2-2
    main.py:
        Solves a 3D truss using Newton-Raphson, Henkey strain, and displacement control.
        Utilizes trussClass but assembly is done manually.

    Prob2-2_GammaVsDisp.png:
        Gamma vs Displacement curve plot output for problem 2-2 part c.

    Prob2-2_EquilibriumPath.png:
        Planner view of equilibrium path for nodes 5 and 6 for problem 2-2 part d.

trussClass.py
    Python class for a truss that is used for both HW problems.



Extra (Not intended as part of Assignment Submission)
------------------------------------------------------------------------------------------------------------------------
Problem 2-1 Brute Force
    main.py:
        Solves the two-truss problem from last week's HW using using Newton-Raphson,
        Henkey strain, and displacement control. Does not use trussClass, and was a direct modificaiton
        of last week's submission. Included in directory just as a reference.

