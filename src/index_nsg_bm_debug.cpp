/*
 * @Author: Chen Shi
 * @Date: 2024-07-05 16:00:24
 * @Description:
 */

#include "efanna2e/index_nsg_bm_debug.h"

#include <omp.h>
#include <bitset>
#include <chrono>
#include <cmath>
#include <boost/dynamic_bitset.hpp>
// @CS0522
#include <iomanip>

#include <efanna2e/exceptions.h>
#include <efanna2e/parameters.h>

namespace efanna2e
{
    IndexNSG_BM_DEBUG::IndexNSG_BM_DEBUG(const size_t dimension, const size_t n, Metric m,
                                   Index *initializer, mi_heap_t *&heap_for_index)
        : IndexNSG_BM(dimension, n, m, initializer, heap_for_index), Index_DEBUG()
    {

    }

    IndexNSG_BM_DEBUG::~IndexNSG_BM_DEBUG()
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

    void IndexNSG_BM_DEBUG::Search(const float *query, const float *x, size_t K,
                                const Parameters &parameters, unsigned *indices)
    {
        // 计时点
        auto search_begin = std::chrono::high_resolution_clock::now();

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

            // 打印导航点
            std::cout << std::endl
                      << "===== DEBUG: Search for query " << search_points_counts.size() << " =====\n"
                      << std::endl;
            std::cout << "Navigate node: " << std::endl;
            std::cout << "id: " << std::setw(6) << get_ep_() << ", dis: " << std::setw(6)
                      << distance_->compare(query, data_ + dimension_ * get_ep_(), (unsigned)dimension_) << std::endl
                      << std::endl;
            std::cout << "Neighbor points of navigate node: " << std::endl;

            // 将导航点的全部邻居放入init_ids
            unsigned tmp_l = 0;
            for (; tmp_l < L && tmp_l < final_graph_[get_ep_()].size(); tmp_l++)
            {
                init_ids[tmp_l] = final_graph_[get_ep_()][tmp_l];
                flags[init_ids[tmp_l]] = true;

                // 放入 init_ids 中的点的 id 和 dist
                auto id = final_graph_[get_ep_()][tmp_l];
                std::cout << "id: " << std::setw(6) << id
                          << ", dis: " << std::setw(6) << distance_->compare(query, data_ + dimension_ * id, (unsigned)dimension_)
                          << std::endl;
            }
            // tmp_l = 50

            // 打印 init_ids 随机选点前的信息
            std::cout << std::endl
                      << "Initial init_ids before random (length = " << tmp_l << "): " << std::endl;
            for (int i = 0; i < tmp_l; ++i)
            {
                std::cout << "id: " << std::setw(6) << init_ids[i] << ", dis: "
                          << std::setw(6) << distance_->compare(query, data_ + dimension_ * init_ids[i], (unsigned)dimension_) << std::endl;
            }

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

            // 计时点
            auto before_greedy_search = std::chrono::high_resolution_clock::now();

            // 打印 greedy search 前的 retset 信息
            std::cout << std::endl
                      << "Initial retset after random and sort (length = " << retset.size() << "): " << std::endl;
            for (int i = 0; i < init_ids.size(); ++i)
            {
                std::cout << "id: " << std::setw(6) << retset[i].id << ", dis: " << std::setw(6) << retset[i].distance << std::endl;
                // retset 中的点的前驱为导航点，路径为 1
                pre[retset[i].id] = get_ep_();
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
            // 统计检索点的数量
            int search_points_count = init_ids.size();
            while (k < (int)L)
            {
                int nk = L;

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
                            continue;
                        flags[id] = 1;
                        float dist =
                            distance_->compare(query, data_ + dimension_ * id, (unsigned)dimension_);

                        // 统计检索点数量
                        ++search_points_count;

                        if (dist >= retset[L - 1].distance)
                            continue;
                        Neighbor nn(id, dist, true);
                        // 修改随机补点的 InsertIntoPool
                        int r = InsertIntoPool(retset.data(), L, nn);
                        // int r = InsertIntoPool(retset, L, nn);

                        // 更新前驱
                        pre[id] = n;
                        std::cout << "id: " << id << ", dis: " << dist << ", pre: " << pre[id] << ", 插入位置: " << r << ", 插入后 retset 长度: " << retset.size() << std::endl;

                        if (r < nk)
                            nk = r;
                    }
                }
                if (nk <= k)
                    k = nk;
                else
                    ++k;
            }

            // 计时点
            auto greedy_search_end = std::chrono::high_resolution_clock::now();

            // 记录搜索结果
            for (size_t i = 0; i < K; i++)
            {
                indices[i] = retset[i].id;
            }

            // 记录最长搜索路径和检索点数量
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

            // 打印导航点
            std::cout << std::endl
                      << "===== DEBUG: Search for query " << search_points_counts.size() << " =====\n"
                      << std::endl;
            std::cout << "Navigate node: " << std::endl;
            std::cout << "id: " << std::setw(6) << get_ep_() << ", dis: " << std::setw(6)
                      << distance_->compare(query, data_ + dimension_ * get_ep_(), (unsigned)dimension_) << std::endl
                      << std::endl;
            std::cout << "Neighbor points of navigate node: " << std::endl;

            // 将导航点的全部邻居放入init_ids
            unsigned tmp_l = 0;
            for (; tmp_l < L && tmp_l < final_graph_[get_ep_()].size(); tmp_l++)
            {
                // init_ids[tmp_l] = final_graph_[get_ep_()][tmp_l];
                // flags[init_ids[tmp_l]] = true;
                init_ids.push_back(final_graph_[get_ep_()][tmp_l]);
                flags[init_ids[tmp_l]] = true;

                // 放入 init_ids 中的点的 id 和 dist
                auto id = final_graph_[get_ep_()][tmp_l];
                std::cout << "id: " << std::setw(6) << id
                          << ", dis: " << std::setw(6) << distance_->compare(query, data_ + dimension_ * id, (unsigned)dimension_)
                          << std::endl;
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

            // 计时点
            auto before_greedy_search = std::chrono::high_resolution_clock::now();

            // 打印 greedy search 前的 retset 信息
            std::cout << std::endl
                      << "Initial retset without random after sort (length = " << retset.size() << "): " << std::endl;
            for (int i = 0; i < retset.size(); ++i)
            {
                std::cout << "id: " << std::setw(6) << retset[i].id << ", dis: " << std::setw(6) << retset[i].distance << std::endl;
                // retset 中的点的前驱为导航点，路径为 1
                pre[retset[i].id] = get_ep_();
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
            // 统计检索点的数量
            int search_points_count = retset.size();
            while (k < (int)L)
            {
                int nk = (L < retset.size()) ? L : retset.size();

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
                            continue;
                        flags[id] = 1;
                        float dist =
                            distance_->compare(query, data_ + dimension_ * id, (unsigned)dimension_);

                        // 统计检索点数量
                        ++search_points_count;

                        if (dist >= retset[(L < retset.size() ? L : retset.size()) - 1].distance)
                            continue;
                        Neighbor nn(id, dist, true);
                        // int r = InsertIntoPool(retset.data(), L, nn);
                        auto r = InsertIntoPool(retset, (L < retset.size()) ? L : retset.size(), nn);

                        // 更新前驱
                        pre[id] = n;
                        std::cout << "id: " << id << ", dis: " << dist << ", pre: " << pre[id] << ", 插入位置: " << r << ", 插入后 retset 长度: " << retset.size() << std::endl;

                        if (r < nk)
                            nk = r;
                    }
                }
                if (nk <= k)
                    k = nk;
                else
                    ++k;
            }

            // 计时点
            auto greedy_search_end = std::chrono::high_resolution_clock::now();

            // 记录搜索结果
            for (size_t i = 0; i < K; i++)
            {
                indices[i] = retset[i].id;
            }

            // 记录最长搜索路径和检索点数量
            // subtitle
            std::cout << std::endl
                      << "===== DEBUG: Statistics =====\n";

            auto max_len = std::max_element(mlen, mlen + 1000000);
            // std::cout << std::endl << "Max search length of current query: " << *max_len << std::endl;
            std::cout << std::endl
                      << *max_len << std::endl;
            this->max_search_lengths.push_back(*max_len);
            // 还原 mlen
            memset(mlen, 0b00000000, sizeof(int) * 1000000);

            // std::cout << "Search points count of current query: " << search_points_count << std::endl;
            std::cout << search_points_count << std::endl;
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
        }
    }

    void IndexNSG_BM_DEBUG::print_info(unsigned query_num, unsigned L, unsigned K)
    {
        std::cout << "\n===== NSG_BM =====\n"
                  << std::endl;
        std::cout << "Queries: " << query_num << std::endl;
        // print hyper parameters
        std::cout << "Hyperparameters: "
                  << "Q (search length) = " << L
                  << ", K (nearest neighbor) = " << K << std::endl;
        std::cout << "Random: " << (NSG_RANDOM ? "enabled" : "disabled") << std::endl;


        std::cout << "\n===== DEBUG: " << (NSG_RANDOM ? "random" : "no random") << " =====\n"
                  << std::endl;
        // 每次搜索的点数
        std::cout << "Search points in Search: " << std::endl;
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
        std::cout << "\nTotal (for " << lengths.size() << " queries): " << total_lengths << std::endl
                  << std::endl;
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