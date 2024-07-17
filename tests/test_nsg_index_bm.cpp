//
// Created by 付聪 on 2017/6/21.
//

#include "efanna2e/index_nsg_bm.h"
#include <efanna2e/util.h>
// @CS0522
#include <iomanip>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <mimalloc-1.8/mimalloc.h>

// TODO 修改为大内存

void load_data(char* filename, float*& data, unsigned& num,
               unsigned& dim) {  // load data with sift10K pattern
  std::ifstream in(filename, std::ios::binary);
  if (!in.is_open()) {
    std::cout << "open file error" << std::endl;
    exit(-1);
  }
  in.read((char*)&dim, 4);
  in.seekg(0, std::ios::end);
  std::ios::pos_type ss = in.tellg();
  size_t fsize = (size_t)ss;
  num = (unsigned)(fsize / (dim + 1) / 4);
  data = new float[(size_t)num * (size_t)dim];

  in.seekg(0, std::ios::beg);
  for (size_t i = 0; i < num; i++) {
    in.seekg(4, std::ios::cur);
    in.read((char*)(data + i * dim), dim * 4);
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


int main(int argc, char** argv) {
  if (argc != 7) {
    std::cout << argv[0] << " data_file nn_graph_path L R C save_graph_file"
              << std::endl;
    exit(-1);
  }

  // @CS0522
  // Heap 内存池
  mi_heap_t *heap_for_index = mi_heap_new();
  // 堆中分配大内存，返回内存指针
  auto big_mem_for_data_ptr = mi_heap_malloc(heap_for_index, BIG_MEM_FOR_DATA_SIZE);

  float* data_load = new (big_mem_for_data_ptr) float;
  unsigned points_num, dim;
  // 用于最后关闭映射
  u_char *data_file_mapping_ptr = nullptr;
  off_t data_file_size;
  // load_data(argv[1], data_load, points_num, dim);
  load_data_mmap(argv[1], data_load, points_num, dim);

  std::string nn_graph_path(argv[2]);
  unsigned L = (unsigned)atoi(argv[3]);
  unsigned R = (unsigned)atoi(argv[4]);
  unsigned C = (unsigned)atoi(argv[5]);

  // data_load = efanna2e::data_align(data_load, points_num, dim);//one must
  // align the data before build
  // efanna2e::IndexNSG index(dim, points_num, efanna2e::L2, nullptr);
  efanna2e::IndexNSG_BM *index = new efanna2e::IndexNSG_BM(dim, points_num, efanna2e::L2, nullptr);

  auto s = std::chrono::high_resolution_clock::now();
  efanna2e::Parameters paras;
  paras.Set<unsigned>("L", L);
  paras.Set<unsigned>("R", R);
  paras.Set<unsigned>("C", C);
  paras.Set<std::string>("nn_graph_path", nn_graph_path);

  // 修改为大内存
  // 堆中分配大内存
  auto big_mem_for_graph_ptr = mi_heap_malloc(heap_for_index, BIG_MEM_FOR_GRAPH_SIZE);
  index->Build(points_num, data_load, paras, heap_for_index);
  auto e = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = e - s;

  std::cout << "indexing time: " << diff.count() << "\n";
  index->Save(argv[6]);

  return 0;
}
