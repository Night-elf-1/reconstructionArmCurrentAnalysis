#include "armCurrentDataAnalysis.hpp"

int main(int argc, char const *argv[])
{
    // 有新的数据往这里填充 下面声明变量出也要相应添加
    std::string filename_1 = "/home/hamster/mycode/reconstructionArmCurrentAnalysis/data/armsdata_J75_11-07(1.05t).txt";    // 1.05吨的电流数据
    std::string filename_2 = "/home/hamster/mycode/reconstructionArmCurrentAnalysis/data/armsdata_J75_11-07(1.2t).txt";     // 1.2吨的电流数据
    std::string filename_3 = "/home/hamster/mycode/reconstructionArmCurrentAnalysis/data/armsdata_J75_11-07(2.2t).txt";     // 2.2吨的电流数据

    // 声明类变量
    extractData extra;
    drawPictures draw;
    filterData filter;
    calculateData cal;
    trainingModel train;

    // 功能选择
    std::cout << "===================================" << std::endl;
    std::cout << "选择功能：" << std::endl;
    std::cout << "1 : 训练模型" << std::endl;
    std::cout << "2 : 画图" << std::endl;
    int siwtchs;
    std::cin >> siwtchs;
    if(siwtchs == 1){
        CurrentDatas1 = extra.extractCurrentData(filename_1);
        CurrentDatas2 = extra.extractCurrentData(filename_2);
        CurrentDatas3 = extra.extractCurrentData(filename_3);

        currentMean1 = extra.extractData_currentMean(filename_1);
        currentMean2 = extra.extractData_currentMean(filename_2);
        currentMean3 = extra.extractData_currentMean(filename_3);

        // 准备训练数据
        // 参数1 —> n-2 为提取到的电流均值容器 参数n-1 和 n为需填充的训练容器
        train.prepareTrainingData(currentMean1, currentMean2, currentMean3, trainingCurrents, trainingWeights);
        // 训练模型
        VectorXd w = train.trainModel(trainingWeights, trainingCurrents);
        std::cout << "===================================" << std::endl;
        cout << "模型训练完成！" << endl;
        cout << "斜率 (w0): " << w(0) << endl;
        cout << "截距 (w1): " << w(1) << endl;
        // 保存权重
        train.saveWeights(w, "/home/hamster/mycode/reconstructionArmCurrentAnalysis/data/weights.txt");
        // 实时预测循环
        char continue_predict = 'y';
        while (continue_predict == 'y' || continue_predict == 'Y')
        {
            double realTimeCurrent;
            std::cout << "===================================" << std::endl;
            cout << "输入实时电流值: ";
            cin >> realTimeCurrent;
            // 预测重量
            double estimatedWeight = w(0) * realTimeCurrent + w(1);
            cout << "预测重量: " << estimatedWeight << " 吨" << endl;

            cout << "是否继续预测？(y/n): ";
            cin >> continue_predict;
        }
        
    }else if(siwtchs == 2){
        std::cout << "===================================" << std::endl;
        std::cout << "1 : 绘制电流夹臂完整数据" << std::endl;
        std::cout << "2 : 绘制单个电流均值数据" << std::endl;
        std::cout << "3 : 绘制直方图" << std::endl;
        int switchs;
        std::cin >> switchs;
        if(switchs == 1){
            CurrentDatas1 = extra.extractCurrentData(filename_1);
            CurrentDatas2 = extra.extractCurrentData(filename_2);
            CurrentDatas3 = extra.extractCurrentData(filename_3);
            filterdata1 = filter.lvbo(CurrentDatas1);
            filterdata2 = filter.lvbo(CurrentDatas2);
            filterdata3 = filter.lvbo(CurrentDatas3);
            draw.drawCurrentPicture(filterdata1[10]);
            draw.drawCurrentPicture(filterdata2[10]);
            draw.drawCurrentPicture(filterdata3[5]);
        }else if(switchs == 2){
            currentMean1 = extra.extractData_currentMean(filename_1);
            currentMean2 = extra.extractData_currentMean(filename_2);
            currentMean3 = extra.extractData_currentMean(filename_3);
            // 获取x轴的信息 画图用
            std::vector<double> x_1(currentMean1.size());                                      // X轴数据
            for (size_t i = 0; i < currentMean1.size(); ++i) {
                x_1[i] = i + 1;
            }
            std::vector<double> x_2(currentMean2.size());                                      // X轴数据
            for (size_t i = 0; i < currentMean2.size(); ++i) {
                x_2[i] = i + 1;
            }
            std::vector<double> x_3(currentMean3.size());                                      // X轴数据
            for (size_t i = 0; i < currentMean3.size(); ++i) {
                x_3[i] = i + 1;
            }
            draw.drawPicture(x_1, currentMean1, currentMean2, x_2, currentMean3, x_3);
        }else if(switchs == 3){
            currentMean1 = extra.extractData_currentMean(filename_1);
            currentMean2 = extra.extractData_currentMean(filename_2);
            currentMean3 = extra.extractData_currentMean(filename_3);
            auto match1 = cal.countOccurrences(currentMean1);
            auto match2 = cal.countOccurrences(currentMean2);
            auto match3 = cal.countOccurrences(currentMean3);
            draw.drawHistogram(match1);
            draw.drawHistogram(match2);
            draw.drawHistogram(match3);
        }else{
            std::cout << "输入有误" << std::endl;
        }
    }else{
        std::cout << "输入有误" << std::endl;
    }


    return 0;
}
