//
// Created by 付聪 on 2017/6/21.
//
#include "efanna2e/index_nsg_bm_debug.h"
#include <efanna2e/util.h>
// @CS0522
#include <iomanip>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <cstring>
#include <mimalloc-1.8/mimalloc.h>

// TODO 修改为大内存

void load_data(char *filename, float *&data, unsigned &num,
               unsigned &dim)
{ // load data with sift10K pattern
  std::ifstream in(filename, std::ios::binary);
  if (!in.is_open())
  {
    std::cout << filename << " open file error" << std::endl;
    exit(-1);
  }
  in.read((char *)&dim, 4);
  std::cout << "data dimension: " << dim << std::endl;
  in.seekg(0, std::ios::end);
  std::ios::pos_type ss = in.tellg();
  size_t fsize = (size_t)ss;
  num = (unsigned)(fsize / (dim + 1) / 4);
  data = new float[num * dim * sizeof(float)];

  in.seekg(0, std::ios::beg);
  for (size_t i = 0; i < num; i++)
  {
    in.seekg(4, std::ios::cur);
    in.read((char *)(data + i * dim), dim * 4);
  }
  in.close();
}

// 使用 mmap 读取大数据集
void load_data_mmap(char *filename, float *&data, unsigned &num, unsigned &dim)
{
  // 获取文件描述符
  int fd = open((const char*)filename, O_RDWR);
  // 打开失败
  if (fd == -1)
  {
    perror(("Open \"" + (std::string)filename + "\" failed. Reason").c_str());
    exit(EXIT_FAILURE);
  }
  // 获取文件信息
  struct stat fs;
  // 获取失败
  if (fstat(fd, &fs) == -1)
  {
    perror("Get file information failed. Reason");
    close(fd);
    exit(EXIT_FAILURE);
  }

  // 文件大小
  off_t file_size = fs.st_size;
  std::cout << "Dataset " << (std::string)filename << " size: " << file_size << std::endl;

  // mmap 映射，获得指针
  u_char *file_mapping_ptr = (u_char *)mmap(nullptr, file_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
  // 映射失败
  if (file_mapping_ptr == MAP_FAILED)
  {
    perror("Memory mapping failed. Reason");
    close(fd);
    exit(EXIT_FAILURE);
  }

  // 读取数据集并保存数据
  long long byte_count = 0LL;
  // 每个向量开头第一个 int 为向量维度
  dim = *(int *)file_mapping_ptr;
  std::cout << "data dimension: " << dim << std::endl;
  // 初始 num 设为 0
  num = 0;
  std::cout << "Reading dataset..." << std::endl;
  while (byte_count < file_size)
  {
    u_char *q = file_mapping_ptr + byte_count;
    // 跳过一个 int 的大小
    q += sizeof(int);
    for (int i = 0; i < dim; i++)
    {
      // 存储到 data 中
      data[num * dim + i] = (float)*(q + i);
    }
    // forward
    byte_count += sizeof(int) + dim * sizeof(u_char);
    num += 1;
  }
  // Completed
  std::cout << "Reading dataset completed. " << std::endl;
  std::cout << "data num: " << num << std::endl;

  // 关闭文件和 mmap 映射
  close(fd);
  munmap(file_mapping_ptr, file_size);
}


void save_result(char *filename, std::vector<std::vector<unsigned>> &results)
{
  std::ofstream out(filename, std::ios::binary | std::ios::out);

  for (unsigned i = 0; i < results.size(); i++)
  {
    unsigned GK = (unsigned)results[i].size();
    out.write((char *)&GK, sizeof(unsigned));
    out.write((char *)results[i].data(), GK * sizeof(unsigned));
  }
  out.close();
}

int main(int argc, char **argv)
{
  if (argc != 7)
  {
    std::cout << argv[0]
              << " data_file query_file nsg_path search_L search_K result_path"
              << std::endl;
    exit(-1);
  }

  // @CS0522
  // Heap 内存池
  mi_heap_t *heap_for_search = mi_heap_new();
  // 堆中分配大内存
  auto big_mem_for_data_ptr = mi_heap_malloc(heap_for_search, BIG_MEM_FOR_DATA_SIZE);
  auto big_mem_for_query_ptr = mi_heap_malloc(heap_for_search, BIG_MEM_FOR_QUERY_SIZE);

  float *data_load = new (big_mem_for_data_ptr) float;
  unsigned points_num, dim;
  // load_data(argv[1], data_load, points_num, dim);
  load_data_mmap(argv[1], data_load, points_num, dim);
  
  float *query_load = new (big_mem_for_query_ptr) float;
  unsigned query_num, query_dim;
  // load_data(argv[2], query_load, query_num, query_dim);
  load_data_mmap(argv[2], query_load, query_num, query_dim);
  assert(dim == query_dim);

  unsigned L = (unsigned)atoi(argv[4]);
  unsigned K = (unsigned)atoi(argv[5]);

  if (L < K)
  {
    std::cout << "search_L cannot be smaller than search_K!" << std::endl;
    exit(-1);
  }

  // data_load = efanna2e::data_align(data_load, points_num, dim);//one must
  // align the data before build query_load = efanna2e::data_align(query_load,
  // query_num, query_dim);
  // efanna2e::L2确定距离比较器，默认欧氏距离平方
  efanna2e::IndexNSG_BM* index = new efanna2e::IndexNSG_BM(dim, points_num, efanna2e::L2, nullptr);
  // DEBUG
  if (DEBUG)
  {
    index = new efanna2e::IndexNSG_BM_DEBUG(dim, points_num, efanna2e::L2, nullptr);
  }
  // 修改为大内存
  // 堆中分配大内存
  auto big_mem_for_graph_ptr = mi_heap_malloc(heap_for_search, BIG_MEM_FOR_GRAPH_SIZE);
  index->Load(argv[3], big_mem_for_graph_ptr);

  efanna2e::Parameters paras;
  paras.Set<unsigned>("L_search", L);
  paras.Set<unsigned>("P_search", L); // 未找到使用P_search的地方

  auto s = std::chrono::high_resolution_clock::now();
  
  // 堆中分配大内存
  auto big_mem_for_res_ptr = mi_heap_malloc(heap_for_search, BIG_MEM_FOR_RES_SIZE);
  std::vector<std::vector<unsigned>> *res = new (big_mem_for_res_ptr) std::vector<std::vector<unsigned>>();

  std::string file_path = NSG_RANDOM ? "nsg_random" : "nsg_no_random";
  std::streambuf *cout_bak;
  std::ofstream *fout;

  for (unsigned i = 0; i < query_num; i++)
  {
    // DEBUG
    // redirect I/O stream
    if (DEBUG)
    {
      fout = new std::ofstream("./anals/" + file_path + "/nsg_query_" + std::to_string(i) + ".txt");
      // rdbuf() 重新定向，返回旧缓冲区指针
      cout_bak = std::cout.rdbuf(fout->rdbuf());
    }

    std::vector<unsigned> tmp(K);
    index->Search(query_load + i * dim, data_load, K, paras, tmp.data());
    res->push_back(tmp);

    // DEBUG
    // recover I/O stream
    if (DEBUG)
    {
      std::cout.rdbuf(cout_bak);
      fout->close();
    }
  }
  auto e = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = e - s;
  // std::cout << "search time: " << diff.count() << "\n";

  // DEBUG
  // redirect I/O stream
  if (DEBUG)
  {
    fout = new std::ofstream("./anals/" + file_path + "/nsg_result.txt");
    // rdbuf() 重新定向，返回旧缓冲区指针
    cout_bak = std::cout.rdbuf(fout->rdbuf());
  }

  std::cout << "search time: " << diff.count() << "\n";

  // print info
  index->print_info(query_num, L, K);

  save_result(argv[6], (*res));

  // 释放内存
  mi_free(big_mem_for_data_ptr);
  mi_free(big_mem_for_query_ptr);
  mi_free(big_mem_for_graph_ptr);
  mi_free(big_mem_for_res_ptr);
  mi_heap_delete(heap_for_search);
  mi_stats_print(NULL);
  // 手动触发一次垃圾回收
  mi_collect(true);

  // DEBUG
  // recover I/O stream
  if (DEBUG)
  {
    std::cout.rdbuf(cout_bak);
    fout->close();
  }

  return 0;
}
