#ifndef D3Q19_H
#define D3Q19_H

namespace ALBRT {

/**
 * D3Q19 Lattice
*/
template <typename T>
struct D3Q19
{

    const static int d = 3;                    ///< Lattice dimension (always 3 for ALBRT)
    const static int q = 19;                   ///< Number of populations
    const static int dataSize = q + d + 1 + d; ///< Data size per cell (populations + density + velocity + force terms)

    const static int rhoIndex = q;           ///< Density index
    const static int uIndex = q + 1;         ///< Velocity index
    const static int forceIndex = q + 1 + d; ///< Force index

    static constexpr int c[q][d] = {
        {0, 0, 0},

        {-1, 0, 0},
        {0, -1, 0},
        {0, 0, -1},
        {-1, -1, 0},
        {-1, 1, 0},
        {-1, 0, -1},
        {-1, 0, 1},
        {0, -1, -1},
        {0, -1, 1},

        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
        {1, 1, 0},
        {1, -1, 0},
        {1, 0, 1},
        {1, 0, -1},
        {0, 1, 1},
        {0, 1, -1}};    ///< Velocity set vectors

    static const int opp[dataSize] = {
        0, 10, 11, 12, 13, 14, 15, 16, 17, 18, 1, 2, 3, 4, 5, 6, 7, 8, 9, 
        19, 20, 21, 22, 23, 24, 25};    ///< Index of opposite vector in velocity set; for macroscopic quantities, the opposite is considered itself.

    static const T t[q] = {
        (T)1 / (T)3,

        (T)1 / (T)18, (T)1 / (T)18, (T)1 / (T)18,
        (T)1 / (T)36, (T)1 / (T)36, (T)1 / (T)36,
        (T)1 / (T)36, (T)1 / (T)36, (T)1 / (T)36,

        (T)1 / (T)18, (T)1 / (T)18, (T)1 / (T)18,
        (T)1 / (T)36, (T)1 / (T)36, (T)1 / (T)36,
        (T)1 / (T)36, (T)1 / (T)36, (T)1 / (T)36};  ///< Lattice weights

    static constexpr T invCs2 = (T)3; ///< Inverse square of lattice speed of sound
};

} //end namespace
#endif
