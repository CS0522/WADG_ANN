'''
Author: Chen Shi
Date: 2024-03-21 14:28:44
Description: Plot class. Draw figures.
'''

import matplotlib.pyplot as plt
# 平滑曲线图
from scipy.interpolate import make_interp_spline
import numpy as np
import csv

'''
绘制 QPS、PRE
'''
class Plot1:
    def __init__(self, query_num:int, wadg_parmeters:str, top_k:int,
                 nsg_csv_file_path, wadg_csv_file_path):
        self.nsg_pre = []
        self.wadg_pre = []

        self.nsg_qps = []
        self.wadg_qps = []
        
        self.query_num = query_num
        self.top_k = top_k
        self.wadg_parmeters = wadg_parmeters

        # read
        self.read_csv(nsg_csv_file_path, 'nsg')
        self.read_csv(wadg_csv_file_path, 'wadg')

    # read csv file
    def read_csv(self, csv_file_path, nsg_or_wadg):
        with open(csv_file_path, 'r') as csv_file:
            reader = csv.reader(csv_file)
            # header is no need
            next(reader)
            for row in reader:
                if (nsg_or_wadg == 'nsg'):
                    self.nsg_pre.append(float(row[0]))
                    self.nsg_qps.append(self.query_num / float(row[1]))
                elif (nsg_or_wadg == 'wadg'):
                    self.wadg_pre.append(float(row[0]))
                    self.wadg_qps.append(self.query_num / float(row[1]))
        
        # print and check
        if (nsg_or_wadg == 'nsg'):
            print("nsg_pre:", self.nsg_pre)
            print("nsg_qps:", self.nsg_qps)
        elif (nsg_or_wadg == 'wadg'):
            print("wadg_pre:", self.wadg_pre)
            print("wadg_qps:", self.wadg_qps)


    # 对比平滑曲线图
    def comparison_curve_chart(self):
        x1 = self.nsg_pre
        x2 = self.wadg_pre

        y1 = self.nsg_qps
        y2 = self.wadg_qps

        x1_array = np.array(x1)
        x2_array = np.array(x2)

        y1_array = np.array(y1)
        y2_array = np.array(y2)

        x1_smooth = np.linspace(x1_array.min(), x1_array.max(), 300)
        x2_smooth = np.linspace(x2_array.min(), x2_array.max(), 300)

        y1_smooth = make_interp_spline(x1_array, y1_array)(x1_smooth)
        y2_smooth = make_interp_spline(x2_array, y2_array)(x2_smooth)

        # 子绘图 1
        # subplot(rows, cols, index)
        # 绘图一行两个
        # plt.subplot(1, 2, 1) 
        # 设置字体
        # plt.rcParams['font.sans-serif'] = ['Times New Roman']
        # x轴标题
        plt.xlabel('Precision@100')
        plt.xlim(0.96, 1)
        # y轴标题
        plt.ylabel('Queries per Second')
        # 绘制平滑曲线图，添加数据点，设置点的大小
        plt.plot(x1_smooth, y1_smooth)  
        plt.plot(x2_smooth, y2_smooth)
        # 曲线上打点
        plt.scatter(x1_array, y1_array, marker = 's')
        plt.scatter(x2_array, y2_array, marker = 's')
        # 设置曲线名称
        plt.legend(['NSG', 'WADG'])
        plt.grid(linestyle = '--', alpha = 0.5)
        # plt.title('SIFT, query_num=' + str(self.query_num) + ', wadg_args=' + self.wadg_parmeters +
        #           ', Top-K=' + str(self.top_k))
        plt.title(f'SIFT, query_num={self.query_num}, wadg_args={self.wadg_parmeters}, Top-K={self.top_k}')

        plt.show()


    # 对比折线图
    def comparison_line_chart(self):
        x1 = self.nsg_pre
        x2 = self.wadg_pre

        y1 = self.nsg_qps
        y2 = self.wadg_qps

        # x轴标题
        plt.xlabel('Precision@100')
        plt.xlim(0.96, 1)
        # y轴标题
        plt.ylabel('Queries per Second')
        # 绘制平滑曲线图，添加数据点，设置点的大小
        plt.plot(x1, y1, marker = 'o')  
        plt.plot(x2, y2, marker = 's')
        # 设置曲线名称
        plt.legend(['NSG', 'WADG'])
        plt.grid(linestyle = '--', alpha = 0.5)
        plt.title(f'SIFT, query_num={self.query_num}, wadg_args={self.wadg_parmeters}, Top-K={self.top_k}')

        plt.show()


    # 对比柱状图
    def comparision_histogram(self):
        base = ['1', '2', '3']
        x1 = [5, 10, 20]
        x2 = [6, 8, 15]

        # 设置字体
        plt.rcParams['font.sans-serif'] = ['Times New Roman']
        # 绘制柱状图
        # x轴标题
        # plt.xlabel('SIFT')  
        # y轴标题
        plt.ylabel('Avg Search Time') 
        base_ticks = range(len(base))
        plt.bar(base_ticks, x1, width = 0.1, label = 'NSG')
        plt.bar([i + 0.1 for i in base_ticks], x2, width = 0.1, label = 'WADG')
        plt.legend(['NSG', 'WADG'])
        # 修改 base 刻度
        plt.xticks(base_ticks, base)
        # 添加网格显示
        plt.grid(linestyle = '--', alpha = 0.5)
        #5、标题
        plt.title("SIFT Average Search Time Comparision")
        plt.show()


'''
根据 debug_output 绘制 检索起点距离、搜索时间对比
'''
class Plot2:
    def __init__(self, exprs: dict, save_figs_base_path: str, num = -1):
        self.save_figs_base_path = save_figs_base_path
        self.exprs = exprs
        self.num = num

        # 检索起点距离
        self.wadg_start_point_dist = []
        self.nsg_random_start_point_dist = []
        self.nsg_no_random_start_point_dist = []

        # 检索点个数
        self.wadg_search_points_count = []
        self.nsg_random_search_points_count = []
        self.nsg_no_random_search_points_count = []

        # 最长搜索路径
        self.wadg_max_search_length = []
        self.nsg_random_max_search_length = []
        self.nsg_no_random_max_search_length = []

        # 搜索时间
        self.wadg_search_time_before_greedy = []
        self.nsg_random_search_time_before_greedy = []
        self.nsg_no_random_search_time_before_greedy = []
        self.wadg_search_time_of_greedy = []
        self.nsg_random_search_time_of_greedy = []
        self.nsg_no_random_search_time_of_greedy = []
        
        # read data and plot
        for expr, file_path in self.exprs.items():
            # read
            self.read_csv(file_path, expr, self.num)
            # plot
            if ('start_point_dist' == expr):
                self.comparison_line_chart(self.wadg_start_point_dist, 
                                           self.nsg_random_start_point_dist, 
                                           self.nsg_no_random_start_point_dist,
                                           expr, self.save_figs_base_path)
                self.diff_line_chart(self.wadg_start_point_dist, 
                                     self.nsg_random_start_point_dist, 
                                     self.nsg_no_random_start_point_dist,
                                     expr, self.save_figs_base_path)
            elif ('search_points_count' == expr):
                self.comparison_line_chart(self.wadg_search_points_count,
                                           self.nsg_random_search_points_count,
                                           self.nsg_no_random_search_points_count,
                                           expr, self.save_figs_base_path)
                self.diff_line_chart(self.wadg_search_points_count,
                                     self.nsg_random_search_points_count,
                                     self.nsg_no_random_search_points_count,
                                     expr, self.save_figs_base_path)
            elif ('max_search_length' == expr):
                self.comparison_line_chart(self.wadg_max_search_length,
                                           self.nsg_random_max_search_length,
                                           self.nsg_no_random_max_search_length,
                                           expr, self.save_figs_base_path)
                self.diff_line_chart(self.wadg_max_search_length,
                                     self.nsg_random_max_search_length,
                                     self.nsg_no_random_max_search_length,
                                     expr, self.save_figs_base_path)
            elif ('search_time_before' == expr):
                    self.comparison_line_chart(self.wadg_search_time_before_greedy,
                                           self.nsg_random_search_time_before_greedy,
                                           self.nsg_no_random_search_time_before_greedy,
                                           expr, self.save_figs_base_path)
                    self.diff_line_chart(self.wadg_search_time_before_greedy,
                                         self.nsg_random_search_time_before_greedy,
                                         self.nsg_no_random_search_time_before_greedy,
                                         expr, self.save_figs_base_path)
            elif ('search_time_of' == expr):
                    self.comparison_line_chart(self.wadg_search_time_of_greedy,
                                           self.nsg_random_search_time_of_greedy,
                                           self.nsg_no_random_search_time_of_greedy,
                                           expr, self.save_figs_base_path)
                    self.diff_line_chart(self.wadg_search_time_of_greedy,
                                         self.nsg_random_search_time_of_greedy,
                                         self.nsg_no_random_search_time_of_greedy,
                                         expr, self.save_figs_base_path)
                    


    '''
    expr: 起点距离、检索点数、最长搜索路径
    num: 前 num 个 query 的运行结果
    '''
    def read_csv(self, file_path: str, expr: str, num = -1):
        with open(file_path, 'r') as csv_file:
            reader = csv.reader(csv_file)
            # head is no need
            next(reader)
            row_count = 0
            for row in reader:
                if (num != -1) and (row_count >= num):
                    break

                if ('start_point_dist' == expr):
                    self.wadg_start_point_dist.append(int(row[0]))
                    self.nsg_random_start_point_dist.append(int(row[1]))
                    self.nsg_no_random_start_point_dist.append(int(row[2]))
                elif ('search_points_count' == expr):
                    self.wadg_search_points_count.append(int(row[0]))
                    self.nsg_random_search_points_count.append(int(row[1]))
                    self.nsg_no_random_search_points_count.append(int(row[2]))
                elif ('max_search_length' == expr):
                    self.wadg_max_search_length.append(int(row[0]))
                    self.nsg_random_max_search_length.append(int(row[1]))
                    self.nsg_no_random_max_search_length.append(int(row[2]))
                elif ('search_time_before' == expr):
                    self.wadg_search_time_before_greedy.append(float(row[0]))
                    self.nsg_random_search_time_before_greedy.append(float(row[1]))
                    self.nsg_no_random_search_time_before_greedy.append(float(row[2]))
                elif ('search_time_of' == expr):
                    self.wadg_search_time_of_greedy.append(float(row[0]))
                    self.nsg_random_search_time_of_greedy.append(float(row[1]))
                    self.nsg_no_random_search_time_of_greedy.append(float(row[2]))
            
                row_count += 1


    '''
    expr: 起点距离、检索点数、最长搜索路径
    '''
    # 对比折线图
    def comparison_line_chart(self, wadg: list, nsg_random: list, nsg_no_random: list, expr: str, save_figs_base_path: str):
        x_coords = list(range(0, len(wadg)))
        
        plt.figure(figsize=(10, 6))
        # x轴标题
        plt.xlabel('Queries')
        # 绘制平滑曲线图
        plt.plot(x_coords, wadg, label = 'wadg')  
        plt.plot(x_coords, nsg_random, label = 'nsg_random')
        plt.legend()
        plt.grid(linestyle = '--', alpha = 0.3)
        plt.title(expr + ': wadg vs. nsg_random')
        # plt.show()
        plt.savefig(save_figs_base_path + f'{expr}_wadg_nsg_random.png')
        plt.clf()

        plt.figure(figsize=(10, 6))
        # x轴标题
        plt.xlabel('Queries')
        # 绘制平滑曲线图
        plt.plot(x_coords, wadg, label = 'wadg')  
        plt.plot(x_coords, nsg_no_random, label = 'nsg_no_random')
        plt.legend()
        plt.grid(linestyle = '--', alpha = 0.3)
        plt.title(expr + ': wadg vs. nsg_no_random')
        # plt.show()
        plt.savefig(save_figs_base_path + f'{expr}_wadg_nsg_no_random.png')
        plt.clf()


    '''
    expr: 起点距离、检索点数、最长搜索路径
    '''
    # 差值折线图（wadg 分别与 nsg_random、nsg_no_random 差值）
    def diff_line_chart(self, wadg: list, nsg_random: list, nsg_no_random: list, expr: str, save_figs_base_path: str):
        diff_wadg_nsg_random = [wadg[i] - nsg_random[i] for i in range(len(wadg))]
        diff_wadg_nsg_no_random = [wadg[i] - nsg_no_random[i] for i in range(len(wadg))]

        x_coords = list(range(0, len(wadg)))
        
        plt.figure(figsize=(10, 6))
        plt.xlabel('Queries')
        plt.plot(x_coords, diff_wadg_nsg_random, label = 'diff_wadg_nsg_random')
        plt.legend(['diff_wadg_nsg_random'])
        plt.grid(linestyle = '--', alpha = 0.3)
        plt.title(expr + ': diff_wadg_nsg_random')
        # plt.show()
        plt.savefig(save_figs_base_path + f'{expr}_diff_wadg_nsg_random.png')
        plt.clf()

        plt.figure(figsize=(10, 6))
        plt.xlabel('Queries')
        plt.plot(x_coords, diff_wadg_nsg_no_random, label = 'diff_wadg_nsg_no_random')
        plt.legend()
        plt.grid(linestyle = '--', alpha = 0.3)
        plt.title(expr + ': diff_wadg_nsg_no_random')
        # plt.show()
        plt.savefig(save_figs_base_path + f'{expr}_diff_wadg_nsg_no_random.png')
        plt.clf()


if __name__ == "__main__":
    # """
    # query_num:
    # wadg_parameters: hot_point_num - window_size - cluster_num
    # """
    # p = Plot1(10000, "2K-2K-K", 100, './nsg.csv', './wadg.csv')
    # # p.comparison_curve_chart()
    # p.comparison_line_chart()

    exprs = {'start_point_dist': './anals/pure_data_start_points.csv',
             'search_points_count': './anals/pure_data_search_points.csv',
             'max_search_length': './anals/pure_data_max_search_length.csv',
             'search_time_before': './anals/pure_data_search_time_before.csv',
             'search_time_of': './anals/pure_data_search_time_of.csv'
            }
    save_figs_base_path = './anals/debug_figs/'
    exprs = {
        
    }

    num = 1000

    p = Plot2(exprs, save_figs_base_path, num)