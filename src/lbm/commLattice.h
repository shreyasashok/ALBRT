#ifndef COMM_LATTICE_H
#define COMM_LATTICE_H

namespace ALBRT {

/**
 * Block fringe communication relationships for face neighbors.
 * Corresponds to p4est z-ordered numbering.
*/
struct FaceCommLattice {
    const static int d = 3; ///< 3-dimensional (ALBRT is always 3D)
    const static int q = 6; ///< 6 directions for face communication

    static constexpr int c[q][d] = { ///< Communication direction
        {-1, 0, 0},
        { 1, 0, 0},
        { 0,-1, 0},
        { 0, 1, 0},
        { 0, 0,-1},
        { 0, 0, 1}
    };

    static constexpr int opp[q] = { ///< Index of opposite communication direction
        1, 0, 3, 2, 5, 4
    };
};

/**
 * Block fringe communication relationships for edge neighbors.
 * Corresponds to p4est z-ordered numbering.
*/
struct EdgeCommLattice {
    const static int d = 3; ///< 3-dimensional (ALBRT is always 3D)
    const static int q = 12; ///< 12 directions for edge communication

    static constexpr int c[q][d] = { ///< Communication direction
        { 0,-1,-1},
        { 0, 1,-1},
        { 0,-1, 1},
        { 0, 1, 1},

        {-1, 0,-1},
        { 1, 0,-1},
        {-1, 0, 1},
        { 1, 0, 1},

        {-1,-1, 0},
        { 1,-1, 0},
        {-1, 1, 0},
        { 1, 1, 0}
    };

    static constexpr int opp[q] = { ///< Index of opposite communication direction
        3, 2, 1, 0, 
        7, 6, 5, 4,
        11, 10, 9, 8
    };
};

/**
 * Block fringe communication relationships for corner neighbors.
 * Corresponds to p4est z-ordered numbering.
*/
struct CornerCommLattice {
    const static int d = 3; ///< 3-dimensional (ALBRT is always 3D)
    const static int q = 8; ///< 8 directions for corner communication

    static constexpr int c[q][d] = { ///< Communication direction
        {-1,-1,-1},
        { 1,-1,-1},

        {-1, 1,-1},
        { 1, 1,-1},

        {-1,-1, 1},
        { 1,-1, 1},

        {-1, 1, 1},
        { 1, 1, 1}
    };

    static constexpr int opp[q] = { ///< Index of opposite communication direction
        7, 6, 5, 4, 3, 2, 1, 0
    };
};




} //end namespace

#endif