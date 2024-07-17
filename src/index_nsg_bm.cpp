#include "efanna2e/index_nsg_bm.h"

#include <omp.h>
#include <bitset>
#include <chrono>
#include <cmath>
#include <boost/dynamic_bitset.hpp>
// @CS0522
#include <iomanip>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <fstream>

#include <efanna2e/exceptions.h>
#include <efanna2e/parameters.h>

namespace efanna2e
{
    IndexNSG_BM::IndexNSG_BM(const size_t dimension, const size_t n, Metric m,
                             Index *initializer)
        : IndexNSG(dimension, n, m, initializer)
    {
    }

    IndexNSG_BM::~IndexNSG_BM()
    {
        if (distance_ != nullptr)
        {
            delete distance_;
            distance_ = nullptr;
        }
        if (initializer_ != nullptr)
        {
            delete initializer_;
            initializer_ = nullptr;
        }
    }

    // mmap 映射读取 nsg 图
    void IndexNSG_BM::Load(const char *filename, void *&big_mem_ptr)
    {
        // final_graph 建立在大内存中
        // final_graph_ptr 的初始化暂时放在了 Load 函数中
        final_graph_ptr = new (big_mem_ptr) std::vector<std::vector<unsigned>>();
        // get file descriptor
        int fd = open(filename, O_RDWR);
        // Open file failed
        if (fd == -1)
        {
            perror(("Open \"" + (std::string)filename + "\" failed. Reason").c_str());
            exit(EXIT_FAILURE);
        }
        // get file info
        struct stat fs;
        // get info failed
        if (fstat(fd, &fs) == -1)
        {
            perror("Get file information failed. Reason");
            close(fd);
            exit(EXIT_FAILURE);
        }
        // file size
        off_t graph_file_size = fs.st_size;
        std::cout << "Graph " << (std::string)filename << " size: " << graph_file_size << std::endl;

        // mmap mapping
        unsigned *graph_file_mapping_ptr = (unsigned *)mmap(nullptr, graph_file_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        if (graph_file_mapping_ptr == MAP_FAILED)
        {
            perror("Memory mapping failed. Reason");
            close(fd);
            exit(EXIT_FAILURE);
        }

        // read graph data
        long long byte_count = 0LL;
        std::cout << "Loading graph..." << std::endl;
        // 先读取 width 和 ep_
        this->width = *graph_file_mapping_ptr;
        this->ep_ = *(graph_file_mapping_ptr + 1);
        byte_count += 2 * sizeof(unsigned);
        unsigned cc = 0;
        unsigned *tmp_ptr = (graph_file_mapping_ptr + 2);
        while (byte_count < graph_file_size)
        {
            unsigned k;
            k = *tmp_ptr;
            byte_count += sizeof(unsigned);
            tmp_ptr += 1;
            if (byte_count >= graph_file_size)
                break;
            cc += k;
            std::vector<unsigned> tmp(k);
            for (int i = 0; i < k; ++i)
            {
                tmp[i] = *(tmp_ptr + i);
            }
            // store in the final_graph_ptr
            final_graph_ptr->emplace_back(tmp);

            // forward
            byte_count += k * sizeof(unsigned);
            tmp_ptr += k;
        }
        // 平均出度（？）
        cc /= nd_;
        // std::cout << cc << std::endl;

        std::cout << "Loading graph completed. " << std::endl;
        close(fd);
        munmap(graph_file_mapping_ptr, graph_file_size);
    }

    void IndexNSG_BM::Search(const float *query, const float *x, size_t K,
                             const Parameters &parameters, unsigned *indices)
    {
        // NSG 随机选点
        if (NSG_RANDOM)
        {
            const unsigned L = parameters.Get<unsigned>("L_search");
            data_ = x;
            std::vector<Neighbor> retset(L + 1);
            std::vector<unsigned> init_ids(L);
            boost::dynamic_bitset<> flags{nd_, 0};
            // std::mt19937 rng(rand());
            // GenRandom(rng, init_ids.data(), L, (unsigned) nd_);

            // 将导航点的全部邻居放入init_ids
            unsigned tmp_l = 0;
            // final_graph
            std::cout << "final graph size: " << final_graph_ptr->size() << std::endl;
            for (; tmp_l < L && tmp_l < (*final_graph_ptr)[get_ep_()].size(); tmp_l++)
            {
                init_ids[tmp_l] = (*final_graph_ptr)[get_ep_()][tmp_l];
                flags[init_ids[tmp_l]] = true;
            }
            // tmp_l = 50

            // 导航点邻居不足L个则随机选取节点，直至init_ids包括L个节点
            while (tmp_l < L)
            {
                unsigned id = rand() % nd_;
                if (flags[id])
                    continue;
                flags[id] = true;
                init_ids[tmp_l] = id;
                tmp_l++;
            }

            // 将init_ids中的节点放入retset作为候选节点集
            for (unsigned i = 0; i < init_ids.size(); i++)
            {
                unsigned id = init_ids[i];
                float dist =
                    distance_->compare(data_ + dimension_ * id, query, (unsigned)dimension_);
                retset[i] = Neighbor(id, dist, true);
                // flags[id] = true;
            }

            std::sort(retset.begin(), retset.begin() + L);

            // greedy search
            int k = 0;
            while (k < (int)L)
            {
                int nk = L;

                if (retset[k].flag)
                {
                    retset[k].flag = false;
                    unsigned n = retset[k].id;

                    for (unsigned m = 0; m < (*final_graph_ptr)[n].size(); ++m)
                    {
                        unsigned id = (*final_graph_ptr)[n][m];

                        if (flags[id])
                            continue;
                        flags[id] = 1;
                        float dist =
                            distance_->compare(query, data_ + dimension_ * id, (unsigned)dimension_);

                        if (dist >= retset[L - 1].distance)
                            continue;
                        Neighbor nn(id, dist, true);
                        int r = InsertIntoPool(retset.data(), L, nn);
                        // TODO 修改随机补点的 InsertIntoPool
                        // int r = InsertIntoPool(retset, L, nn);

                        if (r < nk)
                            nk = r;
                    }
                }
                if (nk <= k)
                    k = nk;
                else
                    ++k;
            }
            // 记录搜索结果
            for (size_t i = 0; i < K; i++)
            {
                indices[i] = retset[i].id;
            }
        }

        // NSG 取消随机选点
        else
        {
            const unsigned L = parameters.Get<unsigned>("L_search");
            data_ = x;
            std::vector<Neighbor> retset;
            std::vector<unsigned> init_ids;
            boost::dynamic_bitset<> flags{nd_, 0};
            // std::mt19937 rng(rand());
            // GenRandom(rng, init_ids.data(), L, (unsigned) nd_);

            // 将导航点的全部邻居放入init_ids
            unsigned tmp_l = 0;
            for (; tmp_l < L && tmp_l < (*final_graph_ptr)[get_ep_()].size(); tmp_l++)
            {
                // init_ids[tmp_l] = (*final_graph_ptr)[get_ep_()][tmp_l];
                // flags[init_ids[tmp_l]] = true;
                init_ids.push_back((*final_graph_ptr)[get_ep_()][tmp_l]);
                flags[init_ids[tmp_l]] = true;
            }
            // tmp_l = 50

            // 不进行随机补点

            // 将init_ids中的节点放入retset作为候选节点集
            for (unsigned i = 0; i < init_ids.size(); i++)
            {
                unsigned id = init_ids[i];
                float dist =
                    distance_->compare(data_ + dimension_ * id, query, (unsigned)dimension_);
                // retset[i] = Neighbor(id, dist, true);
                // flags[id] = true;
                retset.push_back(Neighbor(id, dist, true));
            }

            // std::sort(retset.begin(), retset.begin() + L);
            std::sort(retset.begin(), retset.end());

            // greedy search
            int k = 0;
            while (k < (int)L)
            {
                int nk = (L < retset.size()) ? L : retset.size();

                if (retset[k].flag)
                {
                    retset[k].flag = false;
                    unsigned n = retset[k].id;

                    for (unsigned m = 0; m < (*final_graph_ptr)[n].size(); ++m)
                    {
                        unsigned id = (*final_graph_ptr)[n][m];

                        if (flags[id])
                            continue;
                        flags[id] = 1;
                        float dist =
                            distance_->compare(query, data_ + dimension_ * id, (unsigned)dimension_);

                        if (dist >= retset[(L < retset.size() ? L : retset.size()) - 1].distance)
                            continue;
                        Neighbor nn(id, dist, true);
                        // int r = InsertIntoPool(retset.data(), L, nn);
                        auto r = InsertIntoPool(retset, (L < retset.size()) ? L : retset.size(), nn);

                        if (r < nk)
                            nk = r;
                    }
                }
                if (nk <= k)
                    k = nk;
                else
                    ++k;
            }
            // 记录搜索结果
            for (size_t i = 0; i < K; i++)
            {
                indices[i] = retset[i].id;
            }
        }
    }

    void IndexNSG_BM::print_info(unsigned query_num, unsigned L, unsigned K)
    {
        std::cout << "\n===== NSG =====\n"
                  << std::endl;
        std::cout << "Queries: " << query_num << std::endl;
        // print hyper parameters
        std::cout << "Hyperparameters: "
                  << "Q (search length) = " << L
                  << ", K (nearest neighbor) = " << K << std::endl;
        std::cout << "Random: " << (NSG_RANDOM ? "enabled" : "disabled") << std::endl;

        std::cout << "\n===== END =====\n";
    }

    // 修改 final_graph_ 为 final_graph_ptr
    void IndexNSG_BM::Save(const char *filename)
    {
        std::ofstream out(filename, std::ios::binary | std::ios::out);
        assert(final_graph_ptr->size() == nd_);

        out.write((char *)&width, sizeof(unsigned));
        out.write((char *)&(this->ep_), sizeof(unsigned));
        for (unsigned i = 0; i < nd_; i++)
        {
            unsigned GK = (unsigned)(*final_graph_ptr)[i].size();
            out.write((char *)&GK, sizeof(unsigned));
            out.write((char *)(*final_graph_ptr)[i].data(), GK * sizeof(unsigned));
        }
        out.close();
    }

    // 修改为大内存
    void IndexNSG_BM::Build(size_t n, const float *data, const Parameters &parameters,
                            mi_heap_t *&heap_for_index)
    {
        std::string nn_graph_path = parameters.Get<std::string>("nn_graph_path");
        unsigned range = parameters.Get<unsigned>("R");
        // 修改为大内存
        // 堆上分配大内存
        auto big_mem_for_graph_ptr = mi_heap_malloc(heap_for_index, BIG_MEM_FOR_GRAPH_SIZE);
        // 在这个函数中进行 final_graph_ptr 初始化
        Load_nn_graph(nn_graph_path.c_str(), big_mem_for_graph_ptr);
        data_ = data;
        // 修改为大内存
        init_graph(parameters);
        // 堆上分配大内存
        auto big_mem_for_cut_graph_ptr = mi_heap_malloc(heap_for_index, BIG_MEM_FOR_GRAPH_SIZE);
        SimpleNeighbor *cut_graph_ = new (big_mem_for_cut_graph_ptr) SimpleNeighbor[nd_ * (size_t)range];
        Link(parameters, cut_graph_);
        final_graph_ptr->resize(nd_);

        for (size_t i = 0; i < nd_; i++)
        {
            SimpleNeighbor *pool = cut_graph_ + i * (size_t)range;
            unsigned pool_size = 0;
            for (unsigned j = 0; j < range; j++)
            {
                if (pool[j].distance == -1)
                    break;
                pool_size = j;
            }
            pool_size++;
            (*final_graph_ptr)[i].resize(pool_size);
            for (unsigned j = 0; j < pool_size; j++)
            {
                (*final_graph_ptr)[i][j] = pool[j].id;
            }
        }

        // 修改为大内存
        tree_grow(parameters);

        unsigned max = 0, min = 1e6, avg = 0;
        for (size_t i = 0; i < nd_; i++)
        {
            auto size = (*final_graph_ptr)[i].size();
            max = max < size ? size : max;
            min = min > size ? size : min;
            avg += size;
        }
        avg /= 1.0 * nd_;
        printf("Degree Statistics: Max = %d, Min = %d, Avg = %d\n", max, min, avg);

        has_built = true;
        delete cut_graph_;
        // 释放大内存
        mi_free(big_mem_for_cut_graph_ptr);
        mi_free(big_mem_for_graph_ptr);
    }

    // mmap 映射读取 nsg 图
    void IndexNSG_BM::Load_nn_graph(const char *filename, void *&big_mem_ptr)
    {
        // final_graph_ptr 建立在大内存中
        // 考虑是否将初始化放到 IndexNSG_BM 的构造函数中
        final_graph_ptr = new (big_mem_ptr) std::vector<std::vector<unsigned>>();
        // get file descriptor
        int fd = open(filename, O_RDWR);
        // open file failed
        if (fd == -1)
        {
            perror(("Open \"" + (std::string)filename + "\" failed. Reason").c_str());
            exit(EXIT_FAILURE);
        }
        // get file info
        struct stat fs;
        // get info failed
        if (fstat(fd, &fs) == -1)
        {
            perror("Get file information failed. Reason");
            close(fd);
            exit(EXIT_FAILURE);
        }
        // file size
        off_t graph_file_size = fs.st_size;
        std::cout << "Graph " << (std::string)filename << " size: " << graph_file_size << std::endl;

        // mmap mapping
        unsigned *graph_file_mapping_ptr = (unsigned *)mmap(nullptr, graph_file_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        if (graph_file_mapping_ptr == MAP_FAILED)
        {
            perror("Memory mapping failed. Reason");
            close(fd);
            exit(EXIT_FAILURE);
        }

        // read graph data
        std::cout << "Loading nn graph..." << std::endl;
        unsigned k = *graph_file_mapping_ptr;
        size_t num = (unsigned)(graph_file_size / (k + 1) / 4);

        final_graph_ptr->resize(num);
        final_graph_ptr->reserve(num);
        unsigned kk = (k + 3) / 4 * 4;
        auto tmp_ptr = graph_file_mapping_ptr;
        for (size_t i = 0; i < num; i++)
        {
            // 跳过一个 k (一个 unsigned)
            tmp_ptr += 1;
            (*final_graph_ptr)[i].resize(k);
            (*final_graph_ptr)[i].reserve(kk);
            for (int j = 0; j < k; ++j)
            {
                (*final_graph_ptr)[i][j] = *(tmp_ptr + j);
            }
            // forward
            tmp_ptr += k;
        }
        // read over
        std::cout << "Loading nn graph completed. " << std::endl;
        close(fd);
        munmap(graph_file_mapping_ptr, graph_file_size);
    }

    void IndexNSG_BM::init_graph(const Parameters &parameters)
    {
        float *center = new float[dimension_];
        for (unsigned j = 0; j < dimension_; j++)
            center[j] = 0;
        for (unsigned i = 0; i < nd_; i++)
        {
            for (unsigned j = 0; j < dimension_; j++)
            {
                center[j] += data_[i * dimension_ + j];
            }
        }
        for (unsigned j = 0; j < dimension_; j++)
        {
            center[j] /= nd_;
        }
        std::vector<Neighbor> tmp, pool;
        this->ep_ = rand() % nd_; // random initialize navigating point
        get_neighbors(center, parameters, tmp, pool);
        this->ep_ = tmp[0].id;
        delete center;
    }

    void IndexNSG_BM::get_neighbors(const float *query, const Parameters &parameter,
                                 std::vector<Neighbor> &retset,
                                 std::vector<Neighbor> &fullset)
    {
        unsigned L = parameter.Get<unsigned>("L");

        retset.resize(L + 1);
        std::vector<unsigned> init_ids(L);
        // initializer_->Search(query, nullptr, L, parameter, init_ids.data());

        boost::dynamic_bitset<> flags{nd_, 0};
        L = 0;
        for (unsigned i = 0; i < init_ids.size() && i < (*final_graph_ptr)[this->ep_].size(); i++)
        {
            init_ids[i] = (*final_graph_ptr)[this->ep_][i];
            flags[init_ids[i]] = true;
            L++;
        }
        while (L < init_ids.size())
        {
            unsigned id = rand() % nd_;
            if (flags[id])
                continue;
            init_ids[L] = id;
            L++;
            flags[id] = true;
        }

        L = 0;
        for (unsigned i = 0; i < init_ids.size(); i++)
        {
            unsigned id = init_ids[i];
            if (id >= nd_)
                continue;
            // std::cout<<id<<std::endl;
            float dist = distance_->compare(data_ + dimension_ * (size_t)id, query,
                                            (unsigned)dimension_);
            retset[i] = Neighbor(id, dist, true);
            // flags[id] = 1;
            L++;
        }

        std::sort(retset.begin(), retset.begin() + L);
        int k = 0;
        while (k < (int)L)
        {
            int nk = L;

            if (retset[k].flag)
            {
                retset[k].flag = false;
                unsigned n = retset[k].id;

                for (unsigned m = 0; m < (*final_graph_ptr)[n].size(); ++m)
                {
                    unsigned id = (*final_graph_ptr)[n][m];
                    if (flags[id])
                        continue;
                    flags[id] = 1;

                    float dist = distance_->compare(query, data_ + dimension_ * (size_t)id,
                                                    (unsigned)dimension_);
                    Neighbor nn(id, dist, true);
                    fullset.push_back(nn);
                    if (dist >= retset[L - 1].distance)
                        continue;
                    int r = InsertIntoPool(retset.data(), L, nn);

                    if (L + 1 < retset.size())
                        ++L;
                    if (r < nk)
                        nk = r;
                }
            }
            if (nk <= k)
                k = nk;
            else
                ++k;
        }
    }

    void IndexNSG_BM::get_neighbors(
        const float *query,
        const Parameters &parameter,
        boost::dynamic_bitset<> &flags,
        std::vector<Neighbor> &retset,
        std::vector<Neighbor> &fullset)
    {
        unsigned L = parameter.Get<unsigned>("L");

        retset.resize(L + 1);
        std::vector<unsigned> init_ids(L);
        // initializer_->Search(query, nullptr, L, parameter, init_ids.data());

        L = 0;
        for (unsigned i = 0; i < init_ids.size() && i < (*final_graph_ptr)[this->ep_].size(); i++)
        {
            init_ids[i] = (*final_graph_ptr)[this->ep_][i];
            flags[init_ids[i]] = true;
            L++;
        }
        while (L < init_ids.size())
        {
            unsigned id = rand() % nd_;
            if (flags[id])
                continue;
            init_ids[L] = id;
            L++;
            flags[id] = true;
        }

        L = 0;
        for (unsigned i = 0; i < init_ids.size(); i++)
        {
            unsigned id = init_ids[i];
            if (id >= nd_)
                continue;
            // std::cout<<id<<std::endl;
            float dist = distance_->compare(data_ + dimension_ * (size_t)id, query,
                                            (unsigned)dimension_);
            retset[i] = Neighbor(id, dist, true);
            fullset.push_back(retset[i]);
            // flags[id] = 1;
            L++;
        }

        std::sort(retset.begin(), retset.begin() + L);
        int k = 0;
        while (k < (int)L)
        {
            int nk = L;

            if (retset[k].flag)
            {
                retset[k].flag = false;
                unsigned n = retset[k].id;

                for (unsigned m = 0; m < (*final_graph_ptr)[n].size(); ++m)
                {
                    unsigned id = (*final_graph_ptr)[n][m];
                    if (flags[id])
                        continue;
                    flags[id] = 1;

                    float dist = distance_->compare(query, data_ + dimension_ * (size_t)id,
                                                    (unsigned)dimension_);
                    Neighbor nn(id, dist, true);
                    fullset.push_back(nn);
                    if (dist >= retset[L - 1].distance)
                        continue;
                    int r = InsertIntoPool(retset.data(), L, nn);

                    if (L + 1 < retset.size())
                        ++L;
                    if (r < nk)
                        nk = r;
                }
            }
            if (nk <= k)
                k = nk;
            else
                ++k;
        }
    }

    void IndexNSG_BM::tree_grow(const Parameters &parameter)
    {
        unsigned root = this->ep_;
        boost::dynamic_bitset<> flags{nd_, 0};
        unsigned unlinked_cnt = 0;
        while (unlinked_cnt < nd_)
        {
            DFS(flags, root, unlinked_cnt);
            // std::cout << unlinked_cnt << '\n';
            if (unlinked_cnt >= nd_)
                break;
            findroot(flags, root, parameter);
            // std::cout << "new root"<<":"<<root << '\n';
        }
        for (size_t i = 0; i < nd_; ++i)
        {
            if ((*final_graph_ptr)[i].size() > width)
            {
                width = (*final_graph_ptr)[i].size();
            }
        }
    }

    void IndexNSG_BM::DFS(boost::dynamic_bitset<> &flag, unsigned root, unsigned &cnt)
    {
        unsigned tmp = root;
        std::stack<unsigned> s;
        s.push(root);
        if (!flag[root])
            cnt++;
        flag[root] = true;
        while (!s.empty())
        {
            unsigned next = nd_ + 1;
            for (unsigned i = 0; i < (*final_graph_ptr)[tmp].size(); i++)
            {
                if (flag[(*final_graph_ptr)[tmp][i]] == false)
                {
                    next = (*final_graph_ptr)[tmp][i];
                    break;
                }
            }
            // std::cout << next <<":"<<cnt <<":"<<tmp <<":"<<s.size()<< '\n';
            if (next == (nd_ + 1))
            {
                s.pop();
                if (s.empty())
                    break;
                tmp = s.top();
                continue;
            }
            tmp = next;
            flag[tmp] = true;
            s.push(tmp);
            cnt++;
        }
    }

    void IndexNSG_BM::findroot(boost::dynamic_bitset<> &flag, unsigned &root,
                               const Parameters &parameter)
    {
        unsigned id = nd_;
        for (unsigned i = 0; i < nd_; i++)
        {
            if (flag[i] == false)
            {
                id = i;
                break;
            }
        }

        if (id == nd_)
            return; // No Unlinked Node

        std::vector<Neighbor> tmp, pool;
        get_neighbors(data_ + dimension_ * id, parameter, tmp, pool);
        std::sort(pool.begin(), pool.end());

        unsigned found = 0;
        for (unsigned i = 0; i < pool.size(); i++)
        {
            if (flag[pool[i].id])
            {
                // std::cout << pool[i].id << '\n';
                root = pool[i].id;
                found = 1;
                break;
            }
        }
        if (found == 0)
        {
            while (true)
            {
                unsigned rid = rand() % nd_;
                if (flag[rid])
                {
                    root = rid;
                    break;
                }
            }
        }
        (*final_graph_ptr)[root].push_back(id);
    }

}
