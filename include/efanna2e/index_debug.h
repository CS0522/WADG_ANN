/*
 * @Author: Chen Shi
 * @Date: 2024-07-05 15:40:06
 * @Description: Index 的 DEBUG 基类
 */
#ifndef EFANNA2E_INDEX_DEBUG_H
#define EFANNA2E_INDEX_DEBUG_H

#include "efanna2e/index.h"

namespace efanna2e
{
    class Index_DEBUG
    {
    protected:
        // 存储搜索时间
        // pair: (time before greedy search, time of greedy search)
        std::vector<std::pair<double, double> > search_times;
        // 存储起点的距离
        // pair: (id, distance)
        std::vector<std::pair<unsigned, float> > start_points;

        // 记录每次 Search 的检索点数量
        std::vector<int> search_points_counts;
        // 记录每次 Search 的最长搜索路径
        std::vector<int> max_search_lengths;
        // 记录前驱的 pre 数组
        // 这里想初始化值为 -1，所以用 int。int 范围应该是够 sift 的
        int *pre;
        // 记录最长搜索路径的 mlen 数组
        int *mlen;

        // @CS0522
        virtual std::vector<std::pair<double, double> > get_search_times() { return search_times; }
        virtual std::vector<std::pair<unsigned, float> > get_start_points() { return start_points; }
        virtual std::vector<int> get_search_points_counts() { return search_points_counts; }
        virtual std::vector<int> get_max_search_lengths() { return max_search_lengths; }

    public:
        // constructor
        explicit Index_DEBUG()
        {
            // 显式执行了 IndexWADG 初始化
            // 初始化 pre 数组
            pre = (int *)malloc(sizeof(int) * 1000000);
            // 值为 -1
            memset(pre, 0b11111111, sizeof(int) * 1000000);

            // 初始化 mlen 数组
            mlen = (int *)malloc(sizeof(int) * 1000000);
            memset(mlen, 0b00000000, sizeof(int) * 1000000);
        }

        virtual ~Index_DEBUG() { }
    };
}

#endif