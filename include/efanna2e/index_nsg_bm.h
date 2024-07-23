/*
 * @Author: Chen Shi
 * @Date: 2024-07-10 11:28:27
 * @Description: 重写 index_nsg，基于虚拟大内存实现
 */
#ifndef EFANNA2E_INDEX_NSG_BM_H
#define EFANNA2E_INDEX_NSG_BM_H

#include "efanna2e/index_nsg.h"
#include <mimalloc-1.8/mimalloc.h>

namespace efanna2e
{
    class IndexNSG_BM : public IndexNSG
    {
    public:
        explicit IndexNSG_BM(const size_t dimension, const size_t n, Metric m,
                             Index *initializer, mi_heap_t *&heap_for_index);

        virtual ~IndexNSG_BM();

        unsigned get_ep_() { return this->ep_; }

        // mmap 映射读取 nsg 图
        virtual void Load(const char *filename, void *&big_mem_ptr);

        virtual void Search(const float *query, const float *x, size_t K,
                            const Parameters &parameters, unsigned *indices);

        virtual void print_info(unsigned, unsigned, unsigned);

        virtual void Build(size_t n, const float *data, const Parameters &parameters);

        virtual void Save(const char *filename);

    protected:
        // 相比于原来的 final_graph_，修改为指针，方便指定初始化地址
        CompactGraph *final_graph_ptr;

        void Load_nn_graph(const char *filename, void *&big_mem_ptr);

        void init_graph(const Parameters &parameters);

        void get_neighbors(
            const float *query,
            const Parameters &parameter,
            std::vector<Neighbor> &retset,
            std::vector<Neighbor> &fullset);
        void get_neighbors(
            const float *query,
            const Parameters &parameter,
            boost::dynamic_bitset<> &flags,
            std::vector<Neighbor> &retset,
            std::vector<Neighbor> &fullset);

        void tree_grow(const Parameters &parameter);

        void DFS(boost::dynamic_bitset<> &flag, unsigned root, unsigned &cnt);
        void findroot(boost::dynamic_bitset<> &flag, unsigned &root, const Parameters &parameter);

        void InterInsert(unsigned n, unsigned range, std::vector<std::mutex> &locks, SimpleNeighbor *cut_graph_);
        void sync_prune(unsigned q, std::vector<Neighbor> &pool, const Parameters &parameter, boost::dynamic_bitset<> &flags, SimpleNeighbor *cut_graph_);
        void Link(const Parameters &parameters, SimpleNeighbor *cut_graph_);

    private:
        unsigned width;
        unsigned ep_;

        mi_heap_t *heap_for_index;
    };
}

#endif