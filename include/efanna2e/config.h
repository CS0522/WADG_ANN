#ifndef EFANNA2E_CONFIG_H
#define EFANNA2E_CONFIG_H

namespace efanna2e
{
// @CS0522
#define NSG_RANDOM true
#define DEBUG false

    // 分配未对齐内存
#define BIG_MEM_FOR_DATA_SIZE 1024 * 1024 * 1024 * 65LL
#define BIG_MEM_FOR_QUERY_SIZE 1024 * 1024 * 1024 * 1LL
#define BIG_MEM_FOR_RES_SIZE 1024 * 1024 * 1024 * 10LL
#define BIG_MEM_FOR_GRAPH_SIZE 1024 * 1024 * 1024 * 10LL
}

#endif