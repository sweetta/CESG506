Tatsuhiko Sweet 04/15/2020
The following files were created for CESG 506 Assignment 1.

Problem 1-1
    main.py:
        Creates a Force-Displacement curve for the simple 1-element truss system. By changing the
        variable 'whatToPlot' from 'P' to 'e' will change this to a Strain vs Disp. curve. The different strains
        are named with a subscript a, b, c, d, which follows the assignment prompt. Upper case 'N' signifies
        equilibrium was calculated in the undeformed configuration; lower case 'n' means equilibrium was calculated
        using the deformed configuration.

    force_vs_disp.png:
        Force vs Displacement curve plot output

    strain_vs_disp.png:
        Strain vs Displacement curve plot output

Problem 1-2
    main.py:
        Uses Newton-Raphson and Henkey strain to find displacements for a given P on the 2-element truss.
        The focus is on finding residuals at each iteration step for different ratios of Pcr (gamma),
        these results are saved in a txt file as well as plotted and saved as a png.
        The same code was used for finding Pcr. This was done by solving for displacements in a given bracket
        of P values, interpreting the results (the P value where snap-through occurs is Pcr), and refining the
        P bracket. The final Pcr value used for the assignment was -0.98171345 kN (only y component).

    helperFunctions.py:
        A small set of helper functions to make main.py easier to read.

    R_vs_IterationSteps.png:
        Answer to part 5 of assignment. Plot of residuals at each iteration step vs cumulative number of steps taken.
        Residuals are plotted on a log scale. The tolerance is indicated with the horizontal black bar.

    residuals.txt:
        Used for part 4 & 5 of assignment. Output for tabulating residuals.
        Formatting is later done in excel before compiling submission pdf.