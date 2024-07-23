#include <fstream>
#include <vector>
#include <iostream>
#include <algorithm>

// 读取 ivecs 文件
void read_ivecs(const std::string &filename, std::vector<std::vector<int>> &data)
{
    std::ifstream input(filename, std::ios::binary);
    if (!input.is_open())
    {
        std::cerr << "无法打开文件: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    while (!input.eof())
    {
        int dim;
        input.read(reinterpret_cast<char *>(&dim), sizeof(int));
        std::vector<int> vec(dim);
        input.read(reinterpret_cast<char *>(vec.data()), dim * sizeof(int));
        if (input)
        {
            data.push_back(std::move(vec));
        }
    }
}

// 读取 bvecs 文件
void read_bvecs(const std::string &filename, std::vector<std::vector<int>> &data)
{
    std::ifstream input(filename, std::ios::binary);
    if (!input.is_open())
    {
        std::cerr << "Open " << filename << " failed. " << std::endl;
        exit(EXIT_FAILURE);
    }

    while (!input.eof())
    {
        int dim;
        input.read((char *)&dim, sizeof(int));
        if (input.eof())
            break;
        
        std::vector<int> vec(dim);
        std::vector<u_char> buff(dim);
        input.read((char *)buff.data(), dim * sizeof(u_char));

        // u_char to float
        for (int i = 0; i < dim; i++)
        {
            vec[i] = (int)(buff[i]);
        }
        if (input)
            data.push_back(std::move(vec));
    }
}

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        std::cerr << "test_cal_top_K_precision {search result path} {ground truth path} {K}" << std::endl;
        return 1;
    }

    std::string result_file = argv[1];
    std::string gt_file = argv[2];
    int K = std::stoi(argv[3]);

    std::vector<std::vector<int>> result_data, gt_data;
    read_ivecs(result_file, result_data);
    read_ivecs(gt_file, gt_data);

    // if (result_data.size() != gt_data.size())
    // {
    //     std::cout << "两个文件的行数不同。" << std::endl;
    //     return 0;
    // }

    int overlap = 0;
    for (size_t i = 0; i < result_data.size(); ++i)
    {
        if (result_data[i].size() < K || gt_data[i].size() < K)
        {
            std::cout << "文件1向量数量为: " << result_data[i].size() << ", 文件2向量数量为: " << gt_data[i].size() << std::endl;
            std::cout << "在第 " << i + 1 << " 行，文件 " << (int)((result_data[i].size() < K) ? 1 : 2) << " 中的向量数量为 " << result_data[i].size() << " 小于K。" << std::endl;
            return -1;
        }
        for (int j = 0; j < K; ++j)
        {
            if (std::find(gt_data[i].begin(), gt_data[i].begin() + K, result_data[i][j]) != gt_data[i].begin() + K)
            {
                ++overlap;
            }
        }
    }

    double accuracy = static_cast<double>(overlap) / (result_data.size() * K);
    std::cout << "Top-" << K << " 准确度: " << accuracy << std::endl;

    return 0;
}
