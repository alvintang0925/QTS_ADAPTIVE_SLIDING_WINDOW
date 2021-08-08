//
//  main.cpp
//  QTS ADAPTIVE SLIDING WINDOW
//
//  Created by 唐健恆 on 2021/3/23.
//  Copyright © 2021 Alvin. All rights reserved.
//
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <string>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cfloat>
#include "date.h"
#include "portfolio.h"


using namespace std;
using namespace filesystem;

#define EXPNUMBER 50
#define ITERNUMBER 10000
#define PORTFOLIONUMBER 10
#define FUNDS 10000000
#define DELTA 0.0004
#define QTSTYPE 2 //QTS 0, GQTS 1, GNQTS 2
#define TRENDLINETYPE 0 //linear 0, quadratic 1
#define MODE 1 //extended 0, moved 1
#define STARTDATE "20100930"
#define TRAINRAGNE 65

string file_dir = "0729_1";

bool readData(string filename, vector<vector<string>> &data_vector, int &size, int &day_number) {
    cout << filename << endl;
    ifstream inFile(filename, ios::in);
    string line;
    vector< vector<string> > temp;

    if (!inFile) {
        cout << "Open file failed!" << endl;
        exit(1);
    }
    while (getline(inFile, line)) {
        istringstream delime(line);
        string s;
        vector<string> line_data;
        while (getline(delime, s, ',')) {
            if (s != "\r") {
                s.erase(remove(s.begin(), s.end(), '\r'), s.end());
                line_data.push_back(s);
            }
        }
        temp.push_back(line_data);
    }
    inFile.close();
    
    size = temp[0].size() - 1;
    day_number = temp.size() - 1;
    data_vector = temp;

    return true;
}

bool readSpeData(string filename, string title, vector<string> &myTrainData_vector, int &myTrainData_size) {
    cout << filename << endl;
    ifstream inFile(filename, ios::in);
    string line;
    vector< vector<string> > data;

    if (!inFile) {
        cout << "Open file failed!" << endl;
        return false;
    }
    
    bool sw = false;
    vector<string> line_data;
    while(getline(inFile, line)){
        istringstream delime(line);
        string s;
        
        while(getline(delime, s, ',')){
            if(sw){
                if(s != "\r"){
                    s.erase(remove(s.begin(), s.end(), '\r'), s.end());
                    line_data.push_back(s);
                }
            }
            if(s == title){
                sw = true;
            }
        }
        sw = false;
    }
    inFile.close();
    
    myTrainData_size = line_data.size();
    myTrainData_vector = line_data;

    return true;
}

string** vectorToArray(vector<vector<string>> data_vector){
    string **d = new string*[data_vector.size()];
    for(int j = 0; j < data_vector.size(); j++){
        d[j] = new string[data_vector[0].size()];
        for(int k = 0; k < data_vector[0].size(); k++){
            d[j][k] = data_vector[j][k];
        }
    }
    return d;
}

string* vectorToArray(vector<string> myTrainData_vector){
    string *d = new string[myTrainData_vector.size()];
    for(int j = 0; j < myTrainData_vector.size(); j++){
        d[j] = myTrainData_vector[j];
    }
    return d;
}

void createStock(Stock* stock_list, int size, int range_day_number, string **data, int start_index, int end_index) {
    for (int j = 0; j < size; j++) {
        stock_list[j].idx = j;
        stock_list[j].init(range_day_number + 1);
        stock_list[j].company_name = data[0][j+1];
        for (int k = start_index - 1; k <= end_index; k++) {
            stock_list[j].price_list[k - start_index + 1] = atof(data[k][j+1].c_str());
        }
    }
}

void createDir(string file_dir){
    create_directory(file_dir);
    create_directory(file_dir + "/train");
    create_directory(file_dir + "/test");
    create_directory(file_dir + "/verify");
}

void copyData(string *data_copy, string **data, int day_number){
    for(int j = 0; j < day_number; j++){
        data_copy[j] = data[j+1][0];
        data_copy[j].resize(6);
    }
}

void setWindow(string target_date, string &start_date, string &end_date, int &start_index, int &end_index, string** data, int day_number, int &range_day_number){

    for(int j = 0; j < day_number; j++){
        if(target_date == data[j+1][0]){
            start_date = data[j+1][0];
            end_date = data[j+TRAINRAGNE][0];
            start_index = j + 1;
            end_index = j + TRAINRAGNE;
        }
    }
    range_day_number = end_index - start_index + 1;
}

void initial(double *b, int size) {
    for (int j = 0; j < size; j++) {
        b[j] = 0.5;
    }
}

void initPortfolio(Portfolio *p, int size, int day_number, Stock *stock_list){
    for(int j = 0; j < PORTFOLIONUMBER; j++){
        p[j].init(size, day_number, FUNDS, stock_list);
    }
}

void initPortfolio(Portfolio *p){
    for(int j = 0; j < PORTFOLIONUMBER; j++){
        p[j].init();
    }
}

void genPortfolio(Portfolio* portfolio_list, Stock* stock_list, int portfolio_number, double *beta_, int n, int i) {
    for (int j = 0; j < portfolio_number; j++) {
        portfolio_list[j].exp = n + 1;
        portfolio_list[j].gen = i + 1;
        portfolio_list[j].stock_number = 0;
        for (int k = 0; k < portfolio_list[j].size; k++) {
            double r = (double)rand() / (double)RAND_MAX;
            if (r > beta_[k]) {
                portfolio_list[j].data[k] = 0;
            }
            else {
                portfolio_list[j].data[k] = 1;
            }
        }
        
        for(int k = 0; k < portfolio_list[j].size; k++){
            if(portfolio_list[j].data[k] == 1){
                portfolio_list[j].stock_id_list[portfolio_list[j].stock_number] = k;
                portfolio_list[j].stock_number++;
            }
        }
    }
}

void gen_testPortfolio(Portfolio* portfolio_list, Stock* stock_list, int portfolio_number, string **data, string *myTrainData, int myTrainData_size) {
    for (int j = 0; j < portfolio_number; j++) {
        for(int k = 0; k < portfolio_list[j].size; k++){
            for(int h = 0; h < myTrainData_size; h++){
                if(data[0][k+1] == myTrainData[h]){
                    portfolio_list[j].data[k] = 1;
                    portfolio_list[j].stock_id_list[portfolio_list[j].stock_number] = k;
                    portfolio_list[j].stock_number++;
                    break;
                }
            }
        }
    }
}

void gen_testPortfolio(Portfolio* portfolio_list, Stock* stock_list, int portfolio_number, string **data, Portfolio &temp_portfolio) {
    for (int j = 0; j < portfolio_number; j++) {
        for(int k = 0; k < portfolio_list[j].size; k++){
            if(k == temp_portfolio.stock_id_list[portfolio_list[j].stock_number] && portfolio_list[j].stock_number < temp_portfolio.stock_number){
                portfolio_list[j].data[k] = 1;
                portfolio_list[j].stock_id_list[portfolio_list[j].stock_number] = k;
                portfolio_list[j].stock_number++;
            }else{
                portfolio_list[j].data[k] = 0;
            }
            
        }
    }
}

void capitalLevel(Portfolio* portfolio_list, int portfolio_number, double funds) {
    for (int j = 0; j < portfolio_number; j++) {
        for (int k = 0; k < portfolio_list[j].stock_number; k++) {
            portfolio_list[j].investment_number[k] = portfolio_list[j].getDMoney() / portfolio_list[j].constituent_stocks[portfolio_list[j].stock_id_list[k]].price_list[0];
            portfolio_list[j].remain_fund[k] = portfolio_list[j].getDMoney() - (portfolio_list[j].investment_number[k] * portfolio_list[j].constituent_stocks[portfolio_list[j].stock_id_list[k]].price_list[0]);
        }
//        portfolio_list[j].total_money[0] = funds;
        for (int k = 0; k < portfolio_list[j].day_number; k++) {
            portfolio_list[j].total_money[k] = portfolio_list[j].getRemainMoney();
            for (int h = 0; h < portfolio_list[j].stock_number; h++) {
                portfolio_list[j].total_money[k] += portfolio_list[j].investment_number[h] * portfolio_list[j].constituent_stocks[portfolio_list[j].stock_id_list[h]].price_list[k + 1];
            }
            if(portfolio_list[j].total_money[k] > portfolio_list[j].capital_highest_point){
                portfolio_list[j].capital_highest_point = portfolio_list[j].total_money[k];
            }
            double DD = (portfolio_list[j].capital_highest_point - portfolio_list[j].total_money[k]) / portfolio_list[j].capital_highest_point;
            if(DD > portfolio_list[j].MDD){
                portfolio_list[j].MDD = DD;
            }
        }
    }
}

void countTrend(Portfolio* portfolio_list, int porfolio_number, double funds) {
    
    for (int j = 0; j < porfolio_number; j++) {
        double sum = 0;
        if (TRENDLINETYPE == 0) {
            //portfolio_list[j].countQuadraticYLine();
            double x = 0;
            double y = 0;
            for (int k = 0; k < portfolio_list[j].day_number; k++) {
                x += (k + 1) * (portfolio_list[j].total_money[k] - funds);
                y += (k + 1) * (k + 1);
            }
            if (portfolio_list[j].stock_number != 0) {
                portfolio_list[j].m = x / y;
            }
            for (int k = 0; k < portfolio_list[j].day_number; k++) {
                double Y;
                Y = portfolio_list[j].getNormalY(k + 1);
                sum += (portfolio_list[j].total_money[k] - Y) * (portfolio_list[j].total_money[k] - Y);
            }
        }
        else if (TRENDLINETYPE == 1) {
            portfolio_list[j].countQuadraticYLine();
            for (int k = 0; k < portfolio_list[j].day_number; k++) {
                double Y;
                Y = portfolio_list[j].getQuadraticY(k + 1);
                sum += (portfolio_list[j].total_money[k] - Y) * (portfolio_list[j].total_money[k] - Y);
            }
            portfolio_list[j].m = (portfolio_list[j].getQuadraticY(portfolio_list[j].day_number) - portfolio_list[j].getQuadraticY(1)) / (portfolio_list[j].day_number - 1);
        }
        
        portfolio_list[j].daily_risk = sqrt(sum / (portfolio_list[j].day_number));

        if (portfolio_list[j].m < 0) {
            portfolio_list[j].trend = portfolio_list[j].m * portfolio_list[j].daily_risk;
        }
        else {
            portfolio_list[j].trend = portfolio_list[j].m / portfolio_list[j].daily_risk;
        }
    }
}

void recordGAnswer(Portfolio* portfolio_list, Portfolio& gBest, Portfolio& gWorst, Portfolio& pBest, Portfolio& pWorst) {
    pBest.copyP(portfolio_list[0]);
    pWorst.copyP(portfolio_list[PORTFOLIONUMBER - 1]);
    for (int j = 0; j < PORTFOLIONUMBER; j++) {
        if (pBest.trend < portfolio_list[j].trend) {
            pBest.copyP(portfolio_list[j]);
        }
        if (pWorst.trend > portfolio_list[j].trend) {
            pWorst.copyP(portfolio_list[j]);
        }
    }
    
    if (gBest.trend < pBest.trend) {
        gBest.copyP(pBest);
    }
    
    if (gWorst.trend > pWorst.trend) {
        gWorst.copyP(pWorst);
    }
}

void adjBeta(Portfolio& best, Portfolio& worst, double *beta_) {
    for (int j = 0; j < best.size; j++) {
        if (QTSTYPE == 2) {
            if (best.data[j] > worst.data[j]) {
                if (beta_[j] < 0.5) {
                    beta_[j] = 1 - beta_[j];
                }
                beta_[j] += DELTA;
            }
            else if (best.data[j] < worst.data[j]) {
                if (beta_[j] > 0.5) {
                    beta_[j] = 1 - beta_[j];
                }
                beta_[j] -= DELTA;
            }
        }
        else if (QTSTYPE == 1) {
            if (best.data[j] > worst.data[j]) {
                beta_[j] += DELTA;
            }
            else if (best.data[j] < worst.data[j]) {
                beta_[j] -= DELTA;
            }
        }
        else {
            if (best.data[j] > worst.data[j]) {
                beta_[j] += DELTA;
            }
            else if (best.data[j] < worst.data[j]) {
                beta_[j] -= DELTA;
            }
        }
    }
}

void recordExpAnswer(Portfolio& expBest, Portfolio& gBest) {
    if (expBest.trend < gBest.trend) {
        expBest.copyP(gBest);
        expBest.answer_counter = 1;
    }
    else if (expBest.trend == gBest.trend) {
        expBest.answer_counter++;
    }
}

string getOutputFilePath(string start_date, string end_date, string file_dir, string type){
    return file_dir + "/" + type + "/" + type + "_" + start_date + "_" + end_date + ".csv";
}

void outputFile(Portfolio& portfolio, string file_name, string **data, int start_index) {
    ofstream outfile;
    outfile.open(file_name, ios::out);
    outfile << setprecision(15);

    outfile << "Iteration," << ITERNUMBER << endl;
    outfile << "Element number," << PORTFOLIONUMBER << endl;
    outfile << "Delta," << DELTA << endl;
    outfile << "Exp times," << EXPNUMBER << endl << endl;
    
    outfile << "Init funds," << portfolio.funds << endl;
    outfile << "Final funds," << portfolio.total_money[portfolio.day_number - 1] << endl;
    outfile << "Real award," << portfolio.total_money[portfolio.day_number - 1] - portfolio.funds << endl << endl;
    
    outfile << "MMD," << portfolio.MDD << endl;
    outfile << "PF," << portfolio.PF << endl << endl;
    
    outfile << "m," << portfolio.m << endl;
    outfile << "Daily_risk," << portfolio.daily_risk << endl;
    outfile << "Trend," << portfolio.trend << endl << endl;
    
    if (TRENDLINETYPE == 0) {
        portfolio.countQuadraticYLine();
        double sum = 0;
        for (int k = 0; k < portfolio.day_number; k++) {
            double Y;
            Y = portfolio.getQuadraticY(k + 1);
            sum += (portfolio.total_money[k] - Y) * (portfolio.total_money[k] - Y);
        }
        double c = (portfolio.getQuadraticY(portfolio.day_number) - portfolio.getQuadraticY(1)) / (portfolio.day_number - 1);
        double d = sqrt(sum / (portfolio.day_number));
        
        outfile << "Quadratic trend line," << portfolio.a << "x^2 + " << portfolio.b << "x + " << FUNDS << endl << endl;
        outfile << "Quadratic m," << c << endl;
        outfile << "Quadratic daily risk," << d << endl;
        if(c < 0){
            outfile << "Quadratic trend," << c * d << endl << endl;
        }else{
            outfile << "Quadratic trend," << c / d << endl << endl;
        }
    }
    else {
        outfile << "Quadratic trend line," << portfolio.a << "x^2 + " << portfolio.b << "x + " << FUNDS << endl << endl;
        double x = 0;
        double y = 1;
        double sum = 0;
        for (int k = 0; k < portfolio.day_number - 1; k++) {
            x += (k + 2) * (portfolio.total_money[k + 1] - portfolio.funds);
            y += (k + 2) * (k + 2);
        }

        double c = x / y;
        for (int k = 0; k < portfolio.day_number; k++) {
            double Y;
            Y = c * (k + 1) + portfolio.funds;
            sum += (portfolio.total_money[k] - Y) * (portfolio.total_money[k] - Y);
        }
        double d = sqrt(sum / (portfolio.day_number));

        outfile << "Linear m," << c << endl;
        outfile << "Linear daily risk," << d << endl;
        if(c < 0){
            outfile << "Linear trend," << c * d << endl << endl;
        }else{
            outfile << "Linear trend," << c / d << endl << endl;
        }
    }

    outfile << "Best generation," << portfolio.gen << endl;
    outfile << "Best experiment," << portfolio.exp << endl;
    outfile << "Best answer times," << portfolio.answer_counter << endl << endl;

    outfile << "Stock number," << portfolio.stock_number << endl;
    outfile << "Stock#,";
    for (int j = 0; j < portfolio.stock_number; j++) {
        outfile << portfolio.constituent_stocks[portfolio.stock_id_list[j]].company_name << ",";
    }
    outfile << endl;

    outfile << "Number,";
    for (int j = 0; j < portfolio.stock_number; j++) {
        outfile << portfolio.investment_number[j] << ",";
    }
    outfile << endl;

    outfile << "Distribue funds,";
    for (int j = 0; j < portfolio.stock_number; j++) {
        outfile << portfolio.getDMoney() << ",";
    }
    outfile << endl;

    outfile << "Remain funds,";
    for (int j = 0; j < portfolio.stock_number; j++) {
        outfile << portfolio.remain_fund[j] << ",";
    }
    outfile << endl;

    for (int j = 0; j < portfolio.day_number; j++) {
        outfile << data[start_index + j][0] << ",";
        for (int k = 0; k < portfolio.stock_number; k++) {
            outfile << (portfolio.constituent_stocks[portfolio.stock_id_list[k]].price_list[j + 1] * portfolio.investment_number[k]) + portfolio.remain_fund[k] << ",";
        }
//        outfile << 2 * portfolio.a * (j+1) + portfolio.b << ",";
        outfile << portfolio.total_money[j] << endl;
    }
    outfile << endl;
    outfile.close();

}

void recordCPUTime(double START, double END){
    double total_time = (END - START) / CLOCKS_PER_SEC;
    ofstream outfile_time;
    string file_name = file_dir + "/" + "time.txt";
    outfile_time.open(file_name, ios::out);
    outfile_time << "total time: " << total_time << " sec" << endl;
}

bool isVerifyFinish(Portfolio* portfolio){
    double slope = 2 * portfolio[0].a * portfolio[0].day_number + portfolio[0].b;
    if(slope < 0){
        return true;
    }else{
        return false;
    }

}

bool isVerifyFinish(Portfolio* portfolio, double limit_funds){
    double upperBound = limit_funds * 1.1;
    double lowerBound = limit_funds * 0.9;
    double myFunds =portfolio[0].total_money[portfolio[0].day_number - 1];
    if( myFunds > upperBound || myFunds < lowerBound){
        return true;
    }else{
        return false;
    }
}

void startTrain(Portfolio &result, Stock *stock_list, int size, int range_day_number){
    
    double *beta_ = new double[size];
    Portfolio expBest(size, range_day_number, FUNDS, stock_list);
    Portfolio gBest(size, range_day_number, FUNDS, stock_list);
    Portfolio gWorst(size, range_day_number, FUNDS, stock_list);
    Portfolio pBest(size, range_day_number, FUNDS, stock_list);
    Portfolio pWorst(size, range_day_number, FUNDS, stock_list);
    Portfolio* portfolio_list = new Portfolio[PORTFOLIONUMBER];
    initPortfolio(portfolio_list, size, range_day_number, stock_list);
    
    for(int n = 0; n < EXPNUMBER; n++){
        cout << "___" << n << "___" << endl;
        gBest.init();
        gWorst.init();
        gBest.trend = 0;
        gWorst.trend = DBL_MAX;
        initial(beta_, size);
        for(int i = 0; i < ITERNUMBER; i++){
            pBest.init();
            pWorst.init();
            initPortfolio(portfolio_list);
            genPortfolio(portfolio_list, stock_list, PORTFOLIONUMBER, beta_, n, i);
            capitalLevel(portfolio_list, PORTFOLIONUMBER, FUNDS);
            countTrend(portfolio_list, PORTFOLIONUMBER, FUNDS);
            recordGAnswer(portfolio_list, gBest, gWorst, pBest, pWorst);
            adjBeta(gBest, pWorst, beta_);
        }
        recordExpAnswer(expBest, gBest);
    }
    
    expBest.print();
    delete[] portfolio_list;
    delete[] beta_;
    result.copyP(expBest);

}

int main(int argc, const char * argv[]) {
    srand(114);
    int size;
    int day_number;
    double START, END;
    string** data;
    vector<vector<string>> data_vector;
    
    START = clock();
    createDir(file_dir);
    readData("2011-2020.csv", data_vector, size, day_number);
    data = vectorToArray(data_vector);
    string target_date = STARTDATE;
    string start_date;
    string end_date;
    int start_index;
    int end_index;
    int range_day_number;
    double current_funds = FUNDS;
    double capital_highest_point = 0;
    double MDD = 0;
    double pos_award = 0;
    double neg_award = 0;
    bool isLastDay = false;
    
    while(!isLastDay){
        setWindow(target_date, start_date, end_date, start_index, end_index,data, day_number, range_day_number);
        cout << MODE << "_" << start_date << " - " << end_date << endl;
        
        Stock* stock_list = new Stock[size];
        createStock(stock_list, size, range_day_number, data, start_index, end_index);
        
        string train_start_date = start_date;
        string train_end_date = end_date;
        int train_start_index = start_index;
        int train_end_index = end_index;
        
        Portfolio result(size, range_day_number, FUNDS, stock_list);
        startTrain(result, stock_list, size, range_day_number);
        if(result.trend != 0){
            outputFile(result, getOutputFilePath(train_start_date, train_end_date, file_dir, "train"), data, start_index);
        }else {
                   target_date = data[end_index - TRAINRAGNE + 2][0];
                   continue;
        }
        delete[] stock_list;
        
        int test_start_index = end_index+1;
        int test_end_index;
        string test_start_date = data[test_start_index][0];
        string test_end_date;
        double limit_funds = result.total_money[result.day_number - 1];
        double verify_funds;
        if(MODE == 1){
             verify_funds = result.total_money[0];
        }else{
            verify_funds = FUNDS;
        }
        
        while(true){
            if(MODE == 1){
                start_index++;
                start_date = data[start_index][0];
            }
            end_index++;
            if(end_index == day_number){
                isLastDay = true;
            }
            range_day_number = end_index - start_index + 1;
            end_date = data[end_index][0];
            stock_list = new Stock[size];
            createStock(stock_list, size, range_day_number, data, start_index, end_index);
            Portfolio *new_portfolio = new Portfolio[1];
            new_portfolio[0].init(size, range_day_number, verify_funds, stock_list);
            gen_testPortfolio(new_portfolio, stock_list, 1, data, result);
            capitalLevel(new_portfolio, 1, FUNDS);
            countTrend(new_portfolio, 1, FUNDS);
            new_portfolio[0].countQuadraticYLine();
            if(MODE == 1){
                verify_funds = new_portfolio[0].total_money[0];
            }
            if(isVerifyFinish(new_portfolio) || isLastDay){
                outputFile(new_portfolio[0], getOutputFilePath(start_date, end_date, file_dir, "verify"), data, start_index);
                test_end_date = end_date;
                test_end_index = end_index;
                delete[] stock_list;
                delete[] new_portfolio;
                break;
            }
            
            delete[] stock_list;
            delete[] new_portfolio;
        }
        
        range_day_number = test_end_index - test_start_index + 1;
        stock_list = new Stock[size];
        createStock(stock_list, size, range_day_number, data, test_start_index, test_end_index);
        Portfolio *new_portfolio = new Portfolio[1];
        new_portfolio[0].init(size, range_day_number, current_funds, stock_list, capital_highest_point, MDD);
        gen_testPortfolio(new_portfolio, stock_list, 1, data, result);
        capitalLevel(new_portfolio, 1, current_funds);
        countTrend(new_portfolio, 1, current_funds);
        current_funds = new_portfolio[0].total_money[new_portfolio[0].day_number - 1];
//        capital_highest_point = new_portfolio[0].capital_highest_point;
//        MDD = new_portfolio[0].MDD;
        if(new_portfolio[0].total_money[new_portfolio[0].day_number - 1] - new_portfolio[0].funds >= 0){
            pos_award += new_portfolio[0].total_money[new_portfolio[0].day_number - 1] - new_portfolio[0].funds;
        }else{
            neg_award += -1 * (new_portfolio[0].total_money[new_portfolio[0].day_number - 1] - new_portfolio[0].funds);
        }
        
        if(neg_award == 0){
            new_portfolio[0].PF = -1;
        }else{
            new_portfolio[0].PF = pos_award / neg_award;
        }
        outputFile(new_portfolio[0], getOutputFilePath(test_start_date, test_end_date, file_dir, "test"), data, test_start_index);
        delete[] new_portfolio;
        delete[] stock_list;
        
        cout << test_start_date << " - " << test_end_date << endl;
        cout << "test days: " << range_day_number << endl;
        cout << endl << endl;
        target_date = data[test_end_index - TRAINRAGNE + 1][0];
    }
    
    END = clock();
    recordCPUTime(START, END);
    
    for(int j = 0; j < data_vector.size(); j++){
        data_vector[j].clear();
        delete[] data[j];
    }
    data_vector.clear();
    delete[] data;
    
    return 0;
}
