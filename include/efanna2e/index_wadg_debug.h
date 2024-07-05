/*
 * @Author: Chen Shi
 * @Date: 2024-07-05 13:02:01
 * @Description: 继承 index_wadg，重写 Search 函数，加入 DEBUG 功能
 */
#ifndef EFANNA2E_INDEX_WADG_DEBUG_H
#define EFANNA2E_INDEX_WADG_DEBUG_H

#include "efanna2e/index_wadg.h"
#include "efanna2e/index_debug.h"

namespace efanna2e
{
    class IndexWADG_DEBUG : public IndexWADG, public Index_DEBUG
    {
    public:
        explicit IndexWADG_DEBUG(const size_t dimension, const size_t n, Metric m, 
                        const unsigned K, Index *initializer);

        virtual ~IndexWADG_DEBUG();

        // 重写 set_lru 函数，添加 DEBUG 打印信息
        virtual void Set_lru();

        // 重写带有 DEBUG 功能的 Search 函数
        virtual void Search(
            const float *query,
            const Parameters &parameters,
            unsigned* &indices);

        // 重写打印信息
        // TODO
        virtual void print_info(unsigned, unsigned, unsigned);
    };
}

#endif