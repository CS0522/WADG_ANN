#ifndef EFANNA2E_INDEX_WADG_BM_H
#define EFANNA2E_INDEX_WADG_BM_H

// BM: big memory

// @CS0522
#include "efanna2e/lru_cache.h"
#include "efanna2e/index_nsg_bm.h"

namespace efanna2e
{
    class IndexWADG_BM : public IndexNSG_BM
    {
    public:
        explicit IndexWADG_BM(const size_t dimension, const size_t n, Metric m,
                           const unsigned K, Index *initializer);

        virtual ~IndexWADG_BM();

        // 记录历史查询请求并实现热点的更新
        // 与NSG相比，使用parameters传递k_search，并且将_data的初始化放在了Set_data函数中
        virtual void Search(
            const float *query,
            const Parameters &parameters,
            unsigned *&indices);

        virtual void Set_data(const float *x);

        // @CS0522
        // 需要在读取 final_graph 之后才能把导航点全部邻居放入 LRU 中
        virtual void Set_lru();

        // 打印信息
        virtual void print_info(unsigned, unsigned, unsigned);

    protected:
        // 记录搜索请求
        // 若请求记录窗口已满则进行热点更新，并删除旧的记录
        virtual void record_query(const float *query);
        // 更新热点
        // virtual void update_hot_points(std::vector<std::vector<unsigned> > &search_res);
        virtual void update_hot_point(unsigned hot_point_id);

    // private:
        unsigned max_hot_points_num;          // 最大热点数
        unsigned window_size;                 // 搜索请求记录窗口大小，远大于 K
        unsigned cluster_num;                 // 聚类中心数
        // std::vector<unsigned> hot_points_lru; // 包含全部有效热点id的LRU队列
        // @CS0522
        LRUCache *hot_points_lru;             // 包含全部有效热点id的LRU队列
        // @CS0522
        // "float *" -> "const float *"
        std::vector<const float *> query_list;      // 窗口内搜索请求记录

        // 记录更新热点次数
        int update_hot_points_count = 0;
    };
}
#endif
