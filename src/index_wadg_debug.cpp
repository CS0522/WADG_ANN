/*
 * @Author: Chen Shi
 * @Date: 2024-07-05 16:00:31
 * @Description:
 */

#include <map>
#include <cstring>
#include "efanna2e/index_wadg_debug.h"
#include <vector>
#include <algorithm>
#include <iomanip>

#include <omp.h>
#include <bitset>
#include <chrono>
#include <cmath>
#include <boost/dynamic_bitset.hpp>

#include <efanna2e/exceptions.h>
#include <efanna2e/parameters.h>

namespace efanna2e
{
    IndexWADG_DEBUG::IndexWADG_DEBUG(const size_t dimension, size_t n, Metric m, const unsigned K,
                                     Index *initializer)
        : IndexWADG(dimension, n, m, K, initializer), Index_DEBUG()
    {
    }

    IndexWADG_DEBUG::~IndexWADG_DEBUG()
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

    void IndexWADG_DEBUG::Set_lru()
    {
        // 初始化 LRU 队列
        // (key, value): (下标, 有效热点 id)
        hot_points_lru = new LRUCache(max_hot_points_num);
        // 将导航点的全部邻居放入 LRU 队列
        unsigned ep = get_ep_();

        // print
        std::cout << "\n===== DEBUG: LRU Init =====\n"
                  << std::endl;
        std::cout << "Navigate node: " << std::endl;
        std::cout << "id: " << ep << std::endl
                  << std::endl;
        std::cout << "Neighbor points of navigate node: " << std::endl;
        for (unsigned i = 0; i < max_hot_points_num && i < final_graph_[ep].size(); i++)
        {
            hot_points_lru->put(final_graph_[ep][i]);

            std::cout << final_graph_[ep][i] << " ";
            // pre 数组中这些点的前驱都设为导航点
            pre[final_graph_[ep][i]] = ep;
        }
        // print lru
        std::cout << std::endl;
        std::cout << "\nInitial LRU: " << std::endl;
        hot_points_lru->print_lru_cache();
        // 与 nsg 中 init_ids 加入的导航点一致
        std::cout << std::endl;
    }

    // DEBUG version
    void IndexWADG_DEBUG::Search(const float *query,
                                 const Parameters &parameters,
                                 unsigned *&indices)
    {
        // 计时点
        auto search_begin = std::chrono::high_resolution_clock::now();

        auto L = parameters.Get<unsigned>("L_search");
        auto K = parameters.Get<unsigned>("K_search");

        std::vector<Neighbor> retset;
        boost::dynamic_bitset<> flags{nd_, 0};
        int lru_size = hot_points_lru->get_size();

        unsigned init_size = (L < lru_size) ? L : lru_size;

        // LRU 全部加入到 retset 中
        for (int i = 0; i < lru_size; ++i)
        {
            unsigned id = hot_points_lru->visit(i);
            float dist = distance_->compare(data_ + dimension_ * id, query, (unsigned)dimension_);

            retset.emplace_back(id, dist, true);
            // set flags
            flags[id] = true;
        }
        // 对 retset 从小到大排序
        // 函数指针
        auto cmp = [](Neighbor &a, Neighbor &b) -> bool
        {
            return a.distance < b.distance;
        };
        std::sort(retset.begin(), retset.end(), cmp);

        // 距离 query 最近的 id 放到 LRU 缓存头部
        hot_points_lru->put(retset[0].id);

        // 计时点
        auto before_greedy_search = std::chrono::high_resolution_clock::now();

        // 打印信息
        std::cout << std::endl
                  << "====== DEBUG: Search for query "
                  << this->search_points_counts.size() << " =====" << std::endl;
        // hot points info
        int num = hot_points_lru->get_size();
        std::cout << "\nHot points LRU (length = " << num << "): " << std::endl;
        for (int i = 0; i < num; ++i)
        {
            unsigned id = hot_points_lru->visit(i);
            float dist = distance_->compare(data_ + dimension_ * id, query, (unsigned)dimension_);
            std::cout << "id: " << std::setw(6) << id << ", dis: " << std::setw(6) << dist << std::endl;
        }
        std::cout << std::endl;
        // initial retset info
        std::cout << "Initial retset (length = " << retset.size() << "): " << std::endl;
        for (int i = 0; i < retset.size(); ++i)
        {
            std::cout << "id: " << std::setw(6) << retset[i].id << ", dis: " << std::setw(6) << retset[i].distance << std::endl;
            // retset 中的点无需更新前驱，因为 LRU 中的点的前驱已经设置
            // 路径为 1
            // pre[retset[i].id] = get_ep_();
            mlen[retset[i].id] = 1;
        }
        // 打印 greedy search 的起点信息，即 retset 的第一个点
        std::cout << "\nStart point (retset first point): " << std::endl;
        std::cout << "id: " << std::setw(6) << retset[0].id << ", dis: " << std::setw(6) << retset[0].distance << std::endl;
        // 存储搜索的起始点
        this->start_points.emplace_back(retset[0].id, retset[0].distance);
        // subtitle
        std::cout << std::endl
                  << "===== DEBUG: Greedy search =====\n"
                  << std::endl;

        // 计时点
        auto greedy_search_begin = std::chrono::high_resolution_clock::now();

        // greedy search
        int k = 0;
        // 检索点数量
        // LRU 都计算了一遍
        int search_points_count = retset.size();
        while (k < (int)L)
        {
            int len = (L < retset.size()) ? L : retset.size();
            int nk = len;

            if (retset[k].flag)
            {
                retset[k].flag = false;
                unsigned n = retset[k].id;

                // 打印当前检索的点
                float curr_dist = distance_->compare(query, data_ + dimension_ * n, (unsigned)dimension_);
                std::cout << "Level: " << std::setw(2) << mlen[n] << " - "
                          << "id: " << std::setw(6) << n << ", dis: "
                          << std::setw(6) << curr_dist << ", pre: " << std::setw(6) << pre[n] << " " << std::endl;

                for (unsigned m = 0; m < final_graph_[n].size(); ++m)
                {
                    unsigned id = final_graph_[n][m];

                    // 更新每个点的最长搜索路径
                    mlen[id] = std::max(mlen[id], mlen[n] + 1);

                    if (flags[id])
                    {
                        continue;
                    }
                    flags[id] = 1;
                    float dist =
                        distance_->compare(query, data_ + dimension_ * id, (unsigned)dimension_);

                    // 统计检索点数量
                    ++search_points_count;

                    if (dist >= retset[len - 1].distance)
                    {
                        continue;
                    }

                    Neighbor nn(id, dist, true);
                    auto r = InsertIntoPool(retset, len, nn);

                    // 更新前驱
                    pre[id] = n;
                    std::cout << "id: " << id << ", dis: " << dist << ", pre: " << pre[id] << ", 插入位置: " << r << ", 插入后 retset 长度: " << retset.size() << std::endl;

                    // update len
                    len = (L < retset.size()) ? L : retset.size();

                    if (r < nk)
                    {
                        // 在 r 位置插入了一个 nn
                        // 回溯
                        nk = r;
                    }
                }
            }
            // 回溯
            if (nk <= k)
            {
                k = nk;
            }
            else
            {
                ++k;
            }
        }

        // 计时点
        auto greedy_search_end = std::chrono::high_resolution_clock::now();

        // 记录搜索结果
        for (auto i = 0; i < K; i++)
        {
            indices[i] = retset[i].id;
        }

        // 记录本次 Search 的最长搜索路径
        // 统计检索点数量
        // subtitle
        std::cout << std::endl
                  << "===== DEBUG: Statistics =====\n";

        auto max_len = std::max_element(mlen, mlen + 1000000);
        std::cout << std::endl
                  << "Max search length of current query: " << *max_len << std::endl;
        // std::cout << std::endl
        //           << *max_len << std::endl;
        this->max_search_lengths.push_back(*max_len);
        // 还原 mlen
        memset(mlen, 0b00000000, sizeof(int) * 1000000);

        std::cout << "Search points count of current query: " << search_points_count << std::endl;
        // std::cout << search_points_count << std::endl;
        this->search_points_counts.push_back(search_points_count);

        // 搜索时间打印
        // Eplased time before greedy search
        auto time_before_greedy_search = (std::chrono::duration<double>(before_greedy_search - search_begin)).count();
        std::cout << "\nEplased time before greedy search: " << time_before_greedy_search << std::endl;
        // Eplased time of greedy search
        auto time_of_greedy_search = (std::chrono::duration<double>(greedy_search_end - greedy_search_begin)).count();
        std::cout << "Eplased time before greedy search: " << time_of_greedy_search << std::endl;

        // 加入到 search times 中
        this->search_times.emplace_back(time_before_greedy_search, time_of_greedy_search);

        // 把 query 结果中，离这个 query 最近的点，当作热点加入到热点队列中
        update_hot_point(indices[0]);
        // 打印该点信息
        std::cout << std::endl
                  << "New hot point id: " << indices[0] << ", dis: " << distance_->compare(data_ + indices[0] * dimension_, query, (unsigned)dimension_)
                  << std::endl;
    }

    void IndexWADG_DEBUG::print_info(unsigned query_num, unsigned L, unsigned K)
    {
        std::cout << "\n===== WADG =====\n"
                  << std::endl;
        std::cout << "Queries: " << query_num << std::endl;
        std::cout << "Update hot points times: " << this->update_hot_points_count << std::endl;

        // print hyper parameters
        std::cout << "Hyperparameters: "
                  << "Q (search length) = " << L
                  << ", K (nearest neighbors) = " << K
                  << ", L (cache length) = " << this->max_hot_points_num << std::endl;


        std::cout << "\n===== DEBUG =====\n"
                  << std::endl;
        // 每次搜索的点数
        std::cout << "Search points count in Search: " << std::endl;
        auto counts = this->search_points_counts;
        int total_counts = 0;
        // 每个 search 的搜索点数
        std::cout << "Each query (total " << counts.size() << " queries): " << std::endl;
        for (int i = 0; i < counts.size(); i++)
        {
            // calculate total counts
            total_counts += counts[i];
            // print each count
            std::cout << counts[i] << ", "[i == (counts.size() - 1)];
        }
        // 总的搜索点数，可以看平均值
        std::cout << "\nTotal (for " << counts.size() << " queries): " << total_counts << std::endl;
        // 最大值
        std::cout << "Max count: " << *std::max_element(counts.begin(), counts.end()) << std::endl
                  << std::endl;

        // 每次搜索的最长路径
        std::cout << "Max search length in Search: " << std::endl;
        auto lengths = this->max_search_lengths;
        int total_lengths = 0;
        // 每个 search 的最长搜索路径
        std::cout << "Each query (total " << lengths.size() << " queries): " << std::endl;
        for (int i = 0; i < lengths.size(); i++)
        {
            total_lengths += lengths[i];
            std::cout << lengths[i] << ", "[i == (lengths.size() - 1)];
        }
        // 总的最长搜索路径，可以看平均值
        std::cout << "\nTotal (for " << lengths.size() << " queries): " << total_lengths << std::endl;
        // 最大值
        std::cout << "Max length: " << *std::max_element(lengths.begin(), lengths.end()) << std::endl
                  << std::endl;

        // 每次搜索的起始点的信息
        std::cout << "Start point of each query: " << std::endl;
        auto points = this->start_points;
        for (int i = 0; i < points.size(); i++)
        {
            std::cout << "id: " << std::setw(6) << points[i].first << ", dis: " << std::setw(6) << points[i].second << std::endl;
        }

        // 分部分搜索时间
        std::cout << "\nPart eplased time of each search (time before greedy search, time of greedy search): " << std::endl;
        auto times = this->search_times;
        for (int i = 0; i < times.size(); ++i)
        {
            std::cout << "(" << times[i].first << ", " << times[i].second << ")\n";
        }

        std::cout << "\n===== END =====\n";
    }
}