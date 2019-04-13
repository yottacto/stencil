#include <iostream>
#include "jacobi.hh"

int main()
{
    ice::jacobi<> jb{129, 65 , 65};
    jb.compute(200);
}

