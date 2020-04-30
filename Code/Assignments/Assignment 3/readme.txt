Tatsuhiko Sweet 04/29/2020
The following files were created for CESG 506 Assignment 3.

Problem 3-1
    main3-1.py:
        Solves the 1 free DOF system using the Arc-length method and creates plots.

    Prob3-1_GammaVsDisp.png:
        Gamma vs Displacement curve plot output for problem 3-1. 'a' refers to alpha.
        Letting a = 0 gave the best results. Results from HW 1 using a closed form solution
        are also included here.

    R_vs_IterationStep.png:
        Plotting |R| over each iteration on a semi-log plot.

    FromHW1
        Various files adapted from HW1 to complete the prompt, mainly for comparing results.

Problem 3-2
    main3-2.py:
        Solves the 2 free DOF system using Arc-length method and creates plots.

    Prob3-2_GammaVsDisp.png:
        Gamma vs Displacement curve plot output. Results from HW2 using displacement control are
        also included for comparison.

    oldMain.py: (Not intended as part of Assignment Submission)
        Solves the 2 free DOF system without using an assembler. Was used for constructing the assembler
        but is now outdated. Left in just in case it may be used as a reference someday.

    FromHW2:
        Various files adapted from HW2 (displacement control) used for comparison.

Problem 3-3
    main3-3.py:
        Solves the large 2D truss system using Arc-length Method and creates all resulting plots.

    Prob3-3_DeformedShape.png:
        Deformation of the truss system plotting on equal axis grid. Shows trajectory of each node as well.

    Prob3-3_GammaVsDisp.png:
        Plot of gamma vs displacement of the top right node.

    ElementsTable.txt:
        Summary of element numbering and node IDs.

trussClass.py:
    Python class for a truss that is used for both HW problems.

trussAssemble.py:
    A function that takes a list of trusses and a list of displacements and returns K_global and F_int.
