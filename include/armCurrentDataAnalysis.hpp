#ifndef ARMCURRENTDATAANALYSIS_HPP
#define ARMCURRENTDATAANALYSIS_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <regex>
#include <map>
#include <Eigen/Dense>
#include "eigen3/Eigen/Core"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace std;
using namespace Eigen;

// 声明变量
// 存放电流数据的三维容器()
extern std::vector<std::vector<std::vector<double>>> CurrentDatas1;
extern std::vector<std::vector<std::vector<double>>> CurrentDatas2;
extern std::vector<std::vector<std::vector<double>>> CurrentDatas3;
// 电流均值容器
extern std::vector<double> currentMean1;
extern std::vector<double> currentMean2;
extern std::vector<double> currentMean3;
// 准备训练数据
extern vector<double> trainingCurrents;
extern vector<double> trainingWeights;
// 滤波容器
extern std::vector<std::vector<std::vector<double>>> filterdata1;
extern std::vector<std::vector<std::vector<double>>> filterdata2;
extern std::vector<std::vector<std::vector<double>>> filterdata3;

class extractData{                              // 从txt文件中提取数据类
    public:
        extractData(){}
        ~extractData(){}
        // 分割字符串
        std::vector<std::string> split(const std::string& s, char delimiter);
        // 提取电流数据(将8个夹臂的电流数据全部提取出来)
        std::vector<std::vector<std::vector<double>>> extractCurrentData(const std::string& filename);
        // 提取梯形面积数据(txt文件中)
        std::vector<double> extractData_trapzArea(const std::string& filename);
        // 提取电流均值数据(txt文件中)
        std::vector<double> extractData_currentMean(const std::string& filename);
        // 提取整个夹臂数据的后10%的数据
        std::vector<std::vector<std::vector<double>>> extractLastTenPercent(const std::vector<std::vector<std::vector<double>>>& CurrentDatas_1);
};

class drawPictures{
    public:
        drawPictures(){}
        ~drawPictures(){}
        // 画图
        void drawPicture(std::vector<double>& x_1, std::vector<double>& data_1, std::vector<double> data_2, std::vector<double>& x_2, 
                        std::vector<double> data_3, std::vector<double>& x_3);
        // 绘制电流图
        void drawCurrentPicture(std::vector<std::vector<double>>& datas);
        // 绘制直方图
        void drawHistogram(const std::map<int, int>& counts);
};

class filterData{
    public:
        filterData(){}
        ~filterData(){}
        // 滤波函数1
        std::vector<std::vector<std::vector<double>>> lvbo(std::vector<std::vector<std::vector<double>>>& datas);
        // 滤波函数2(agv上跑)
        std::vector<std::vector<double>> lvbo2(std::vector<std::vector<double>>& datas);
};

class calculateData{
    public:
        calculateData(){}
        ~calculateData(){}
        // 四舍五入并取整
        int roundAndTruncate(double value);
        // 统计每个整数出现的次数
        std::map<int, int> countOccurrences(const std::vector<double>& data);
        // 计算电流均值
        std::vector<double> calculateCurrentMean(std::vector<std::vector<std::vector<double>>>& tempCurrentDatas);
};

class trainingModel{
    public:
        trainingModel(){}
        ~trainingModel(){}
        // 训练模型
        VectorXd trainModel(const vector<double>& weights, const vector<double>& areas);
        // 准备训练数据
        void prepareTrainingData(const vector<double>& currentMean_1, const vector<double>& currentMean_2,
                                const vector<double>& currentMean_3,
                                vector<double>& training_currents,
                                vector<double>& training_weights);
        // 保存训练结果
        void saveWeights(const VectorXd& w, const string& filename);
};

#endif // !ARMCURRENTDATAANALYSIS_HPP#define 