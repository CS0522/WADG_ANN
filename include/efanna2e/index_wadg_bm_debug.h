/*
 * @Author: Chen Shi
 * @Date: 2024-07-05 13:02:01
 * @Description: 继承 index_wadg，重写 Search 函数，加入 DEBUG 功能
 */
#ifndef EFANNA2E_INDEX_WADG_BM_DEBUG_H
#define EFANNA2E_INDEX_WADG_BM_DEBUG_H

#include "efanna2e/index_wadg_bm.h"
#include "efanna2e/index_debug.h"

namespace efanna2e
{
    class IndexWADG_BM_DEBUG : public IndexWADG_BM, public Index_DEBUG
    {
    public:
        explicit IndexWADG_BM_DEBUG(const size_t dimension, const size_t n, Metric m, 
                        const unsigned K, Index *initializer);

        virtual ~IndexWADG_BM_DEBUG();

        // TODO 需要重写其他函数
        virtual void Set_lru();

        virtual void Search(const float *query,
                                 const Parameters &parameters,
                                 unsigned *&indices);

        virtual void print_info(unsigned query_num, unsigned L, unsigned K);
    };
}

#endif