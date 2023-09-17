#include <iostream>
#include <gtest/gtest.h>

#include "ALBRT.h"

int main(int, char**){
    std::cout << "Hello, from ALBRT EsoPullTest!\n";
}

TEST(TestTest, Case1) {
    EXPECT_EQ(1.0, 1.0);
}