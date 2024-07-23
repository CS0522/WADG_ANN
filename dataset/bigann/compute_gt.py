import numpy as np
import struct
import faiss

# read bvecs file
def read_bvecs(filename: str, first_n = -1):
    data = []
    num = 0
    with open(filename, 'rb') as f:
        while True:
            if (first_n != -1 and num >= first_n):
                break
            # read dimension
            dim_bytes = f.read(4)
            if not dim_bytes:
                break
            dim = struct.unpack('i', dim_bytes)[0]
            vec = np.fromfile(f, dtype=np.uint8, count=dim)
            data.append(vec)
            num += 1
    return np.array(data), dim, num


# compute groundtruth
def compute_gt_by_faiss(base_data: np.ndarray, query_data:np.ndarray, dim, gt_filename: str, top_k = 100):
    # 创建 faiss 索引
    index = faiss.IndexFlatL2(dim)
    index.add(base_data)

    # 计算 groundtruth
    D, I = index.search(query_data, top_k)

    # 保存 groundtruth
    groundtruth = I.astype(np.int32)
    print(f"groundtruth: num = {len(groundtruth)}, k = {len(groundtruth[0])}")
    # 写入 groundtruth: k(int) + k * index(int)
    # 添加一列
    dim_col = np.full(len(groundtruth), len(groundtruth[0]), dtype=np.int32)
    groundtruth = np.insert(groundtruth, 0, dim_col, axis=1)
    groundtruth.tofile(gt_filename)


if __name__ == "__main__":
    bigann_learn_small_data, dim, num = read_bvecs('./dataset/bigann/bigann_learn_small.bvecs')
    print(f"base dataset: dim = {dim}, num = {num}")
    bigann_query_data, query_dim, query_num = read_bvecs('./dataset/bigann/bigann_query.bvecs')
    print(f"queries: dim = {query_dim}, num = {query_num}")
    assert (dim == query_dim), "Dimensions of base data and queries are different. "
    compute_gt_by_faiss(bigann_learn_small_data, bigann_query_data, dim, 
                        './dataset/bigann/bigann_learn_small_groundtruth_faiss.ivecs')