# Parellel & High Performance Computing - Project Choice - Louis Stenger

I will focus on the second proposition. The $N$-Body Gravity Problem will be
solved numerically in 2D using the Fast Multipole Method (FFM, on a
hiearchically decomposed domain, 2x2 splitting per level). Once a funtional
baseline program is available, the solver will be optimized and parallelized
with the MPI library.

## Notes & Questions

I would appreciate a (quick if possible) feedback on the following items (at least the
first) below.

* Searching information online, I saw most references ([here][1], [there][2] or
    even [here][3]) compute the gravitational potential instead of the force
    directly. Is this a possibility? The force would be computed with finite
    differences each time a force kick is required by temporal evolution, and
    this without increasing numerical complexity (if I'm not mistaken; the FMM
    has complexity $\mathbb{O}(n)$).
* I think the projects were handed in WAY too late. Why wait until the lecture
    on profiling was given? I could really have used some time during the
    holidays just to get a baseline code running, and _advanced_ knowledge in
    profiling would likely not have been needed.
* Also concerning the projects, how were the 50 hours of work estimated? I've
    been spending at least 5 to get familiar with the different algorithms
    available for problem 2, will likely need 10h minimum to write up my final
    report. I don't really see how I should get an optimized parallelized code
    running in 35h... (At least, not within the very bsuy end of the semester
    that is likely to begin.)

Best, Louis (<louis.stenger@epfl.ch>)

[1]: https://obswww.unige.ch/lastro/conferences/sf2013/pdf/lecture2.pdf
[2]: https://math.nyu.edu/faculty/greengar/shortcourse_fmm.pdf
[3]: https://www.youtube.com/watch?v=qMLIyZi8Sz0
