#include <iostream>
#include "jacobi.hh"

int main()
{
    ice::jacobi<> jb{128, 64 , 64};
    jb.compute(20);
}

