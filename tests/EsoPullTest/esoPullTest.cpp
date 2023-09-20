#include <iostream>
#include <gtest/gtest.h>

#include "ALBRT.h"

using namespace ALBRT;

int main(int argc, char** argv){
    ::testing::InitGoogleTest(&argc, argv);
    std::cout << "Hello, from ALBRT EsoPullTest!\n";
    return RUN_ALL_TESTS();
}

TEST(EsoPullTestCase, D3Q19Double) {
    using T = double;
    using Lattice = D3Q19<T>;
    using Octant_T = Octant<T, Lattice, 20, 24, 28, 2>;
    Octant_T octant;

    int testCell[3] = {5,5,5};
    int sourceCellIndex = Octant_T::getCellIndex(testCell[0], testCell[1], testCell[2]);

    for (int iPop = 0; iPop < Lattice::q; iPop++) {
        int destCell[3] = {testCell[0] + Lattice::c[iPop][0],
                           testCell[1] + Lattice::c[iPop][1],
                           testCell[2] + Lattice::c[iPop][2]};

        int destCellIndex = Octant_T::getCellIndex(destCell[0], destCell[1], destCell[2]);

        T * popPointer = &octant.getEvenData()[Lattice::opp[iPop]][sourceCellIndex];
        T * popOppPointer = &octant.getOddData()[iPop][destCellIndex];
        EXPECT_EQ(popPointer, popOppPointer) << "iPop " << iPop << " pointers don't match." << std::endl;
    }
}

TEST(EsoPullTestCase, D3Q19Float) {
    using T = float;
    using Lattice = D3Q19<T>;
    using Octant_T = Octant<T, Lattice, 20, 24, 28, 2>;
    Octant_T octant;

    int testCell[3] = {5,5,5};
    int sourceCellIndex = Octant_T::getCellIndex(testCell[0], testCell[1], testCell[2]);

    for (int iPop = 0; iPop < Lattice::q; iPop++) {
        int destCell[3] = {testCell[0] + Lattice::c[iPop][0],
                           testCell[1] + Lattice::c[iPop][1],
                           testCell[2] + Lattice::c[iPop][2]};

        int destCellIndex = Octant_T::getCellIndex(destCell[0], destCell[1], destCell[2]);

        T * popPointer = &octant.getEvenData()[Lattice::opp[iPop]][sourceCellIndex];
        T * popOppPointer = &octant.getOddData()[iPop][destCellIndex];
        EXPECT_EQ(popPointer, popOppPointer) << "iPop " << iPop << " pointers don't match." << std::endl;
    }
}

TEST(EsoPullTestCase, D3Q19Mixed) {
    using T = float;
    using S = double;
    using Lattice = D3Q19<S>;
    using Octant_T = Octant<T, Lattice, 20, 24, 28, 2>;
    Octant_T octant;

    int testCell[3] = {5,5,5};
    int sourceCellIndex = Octant_T::getCellIndex(testCell[0], testCell[1], testCell[2]);

    for (int iPop = 0; iPop < Lattice::q; iPop++) {
        int destCell[3] = {testCell[0] + Lattice::c[iPop][0],
                           testCell[1] + Lattice::c[iPop][1],
                           testCell[2] + Lattice::c[iPop][2]};

        int destCellIndex = Octant_T::getCellIndex(destCell[0], destCell[1], destCell[2]);

        T * popPointer = &octant.getEvenData()[Lattice::opp[iPop]][sourceCellIndex];
        T * popOppPointer = &octant.getOddData()[iPop][destCellIndex];
        EXPECT_EQ(popPointer, popOppPointer) << "iPop " << iPop << " pointers don't match." << std::endl;
    }
}