#include <iostream>
#include "jacobi.hh"
#include "config.hh"

int main()
{
    auto config{ice::config{}};
    auto nround = config.front().nround;
    auto len_i = config.front().len_i;
    auto len_j = config.front().len_j;
    auto len_k = config.front().len_k;

    ice::jacobi<> jb{len_i, len_j , len_k};

    // warmup
    for (auto i = 0; i < 4; i++)
        jb.compute(nround);

    jb.compute(nround, true);
}

