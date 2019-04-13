#pragma once
#include <algorithm>
#include <array>
#include <vector>
#include <mpi.h>

#include <iostream>

namespace ice
{

struct point
{
    int x;
    int y;
    int z;
};

template <class T = float>
struct jacobi
{
    using value_type = T;

    jacobi() = default;

    jacobi(int len_i, int len_j, int len_k) :
        len_i(len_i), len_j(len_j), len_k(len_k),
        n(len_i * len_j * len_k)
    {
        if (!MPI::Is_initialized())
            MPI::Init();
        size = MPI::COMM_WORLD.Get_size();
        rank = MPI::COMM_WORLD.Get_rank();

        // block_size = n / size;
        // auto extra = n % size;
        // start      = rank       * block_size + std::min(rank,     extra);
        // end        = (rank + 1) * block_size + std::min(rank + 1, extra);
        block_size = (n + size - 1) / size;
        start      = rank       * block_size;
        end        = (rank + 1) * block_size;

        value.resize(block_size * size);
        next_ss.resize(block_size * size);

        for (auto i = 0; i < len_i; i++)
            for (auto j = 0; j < len_j; j++)
                for (auto k = 0; k < len_k; k++)
                    value[id(i, j, k)] = static_cast<value_type>(k * k)
                        / static_cast<value_type>((len_k-1) * (len_k-1));
    }

    ~jacobi()
    {
        if (!MPI::Is_finalized())
            MPI::Finalize();
    }

    auto id(int i, int j, int k) const -> int
    {
        return i * (len_j * len_k) + j * len_k + k;
    }

    auto coordinate(int id) const -> point
    {
        auto i = id / (len_j * len_k);
        id %= (len_j * len_k);
        auto j = id / len_k;
        auto k = id % len_k;
        return point{i, j, k};
    }

    auto on_border(point const& p) const -> bool
    {
        return (p.x == 0 || p.x == len_i - 1)
            || (p.y == 0 || p.y == len_j - 1)
            || (p.z == 0 || p.z == len_k - 1);
    }

    void compute(int nround)
    {
        value_type gosa{0.0};
        for (auto r = 0; r < nround; r++) {
            for (auto i = start; i < std::min(end, n); i++) {
                auto p = coordinate(i);
                if (on_border(p)) {
                    // FIXME
                    continue;
                }
                auto x = p.x;
                auto y = p.y;
                auto z = p.z;
                auto s0 = a[0] * value[id(x + 1, y, z)]
                    + a[1] * value[id(x, y + 1, z)]
                    + a[2] * value[id(x, y, z+1)]
                    + b[0] * (value[id(x+1, y+1, z  )] - value[id(x+1, y-1, z  )]
                        - value[id(x-1, y+1, z  )] + value[id(x-1, y-1, z  )])
                    + b[1] * (value[id(x,   y+1, z+1)] - value[id(x,   y-1, z+1)]
                        - value[id(x,   y+1, z-1)] + value[id(x,   y-1, z-1)])
                    + b[2] * (value[id(x+1, y,   z+1)] - value[id(x-1, y,   z+1)]
                        - value[id(x+1, y,   z-1)] + value[id(x-1, y,   z-1)])
                    + c[0] * value[id(x-1, y,   z  )]
                    + c[1] * value[id(x,   y-1, z  )]
                    + c[2] * value[id(x,   y,   z-1)]
                    + wrk1;

                auto ss = (s0 * a[3] - value[id(x, y, z)]) * bnd;
                next_ss[i] = ss;
            }

            // TODO value_type mapto MPI::DATATYPE
            MPI::COMM_WORLD.Allgather(
                MPI::IN_PLACE, 0, MPI::DATATYPE_NULL,
                next_ss.data(), block_size, MPI::FLOAT
            );

            gosa = 0.0;
            for (auto i = 1; i < len_i - 1; i++)
                for (auto j = 1; j < len_j - 1; j++)
                    for (auto k = 1; k < len_k - 1; k++) {
                        auto ss = next_ss[id(i, j, k)];
                        gosa += ss*ss;
                        value[id(i, j, k)] = value[id(i, j, k)] + omega * ss;
                    }

            if (!rank)
                std::cerr << "iter " << r << ", gosa: " << gosa << "\n";
        }
    }

    // constant
    std::array<value_type, 4> a{1.0, 1.0, 1.0, 1.0/6.0};
    std::array<value_type, 3> b{0.0, 0.0, 0.0};
    std::array<value_type, 3> c{1.0, 1.0, 1.0};
    value_type wrk1{0.0};
    value_type bnd{1.0};
    value_type omega{0.8};

    int size;
    int rank;
    int start;
    int end;
    int block_size;
    int len_i;
    int len_j;
    int len_k;
    int n;

    std::vector<value_type> value;
    std::vector<value_type> next_ss;
};

} // namespace ice

