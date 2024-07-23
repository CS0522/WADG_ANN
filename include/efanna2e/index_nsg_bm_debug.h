/*
 * @Author: Chen Shi
 * @Date: 2024-07-05 13:02:01
 * @Description: 继承 index_wadg，重写 Search 函数，加入 DEBUG 功能
 */
#ifndef EFANNA2E_INDEX_NSG_BM_DEBUG_H
#define EFANNA2E_INDEX_NSG_BM_DEBUG_H

#include "efanna2e/index_nsg_bm.h"
#include "efanna2e/index_debug.h"

namespace efanna2e
{
    class IndexNSG_BM_DEBUG : public IndexNSG_BM, public Index_DEBUG
    {
    public:
        explicit IndexNSG_BM_DEBUG(const size_t dimension, const size_t n, Metric m,
                                   Index *initializer, mi_heap_t *&heap_for_index);

        virtual ~IndexNSG_BM_DEBUG();

        // TODO 需要重写其他函数
        virtual void Search(const float *query, const float *x, size_t K,
                            const Parameters &parameters, unsigned *indices);

        virtual void print_info(unsigned query_num, unsigned L, unsigned K);
    };
}

#endif