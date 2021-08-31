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
#define PARTICLENUMBER 10
#define FUNDS 10000000.0
#define DELTA 0.0004
#define QTSTYPE 2 //QTS 0, GQTS 1, GNQTS 2
#define TRENDLINETYPE 0 //linear 0, quadratic 1
#define MODE 0 //IC_NL 0, IC_LTT 1, DC_NL 2
#define STOPTYPE 0 //loss & profit 0, linear 1, quadratic 2, linear + L&P 3
#define STARTDATE "20100930"
#define TRAINRANGE 65

double LOWER = 0.2;
double UPPER = 0.2;
string FILE_DIR = "0824_IC_NL_GNQ_LN_LN+L&P_test";
string DATA_FILE_NAME = "2011-2020.csv";

class TradePeriod{
public:
    string* date_list = NULL;
    
    int start_index = 0;
    int end_index = 0;
    int train_start_index = 0;
    int train_end_index = 0;
    int train_day_number = 0;
    int verify_start_index = 0;
    int verify_end_index = 0;
    int verify_day_number = 0;
    int test_start_index = 0;
    int test_end_index = 0;
    int test_day_number = 0;
    
    Portfolio train_result;
    Portfolio verify_result;
    Portfolio test_result;
    
    void init(int, string**);
    string getTrainStartDate();
    string getTrainEndDate();
    string getVerifyStartDate();
    string getVerifyEndDate();
    string getTestStartDate();
    string getTestEndDate();
    int getTrainDayNumber();
    int getVerifyDayNumber();
    int getTestDayNumber();
    
    string getDate(int);
    TradePeriod();
    ~TradePeriod();
    
};

TradePeriod::TradePeriod(){}
TradePeriod::~TradePeriod(){
    if(this -> date_list != NULL){
        delete[] date_list;
    }
    date_list = NULL;
}

void TradePeriod::init(int day_number, string** data){
    if(this -> date_list != NULL){
        delete[] date_list;
    }
    this -> date_list = new string[day_number];
    for(int j = 0; j < day_number; j++){
        this -> date_list[j] = data[j+1][0];
    }
}

string TradePeriod::getTrainStartDate(){
    return date_list[train_start_index];
}

string TradePeriod::getTrainEndDate(){
    return date_list[train_end_index];
}

string TradePeriod::getVerifyStartDate(){
    return date_list[verify_start_index];
}

string TradePeriod::getVerifyEndDate(){
    return date_list[verify_end_index];
}

string TradePeriod::getTestStartDate(){
    return date_list[test_start_index];
}

string TradePeriod::getTestEndDate(){
    return date_list[test_end_index];
}

string TradePeriod::getDate(int index){
    return date_list[index];
}

int TradePeriod::getTrainDayNumber(){
    return this -> train_end_index - this -> train_start_index + 1;
}

int TradePeriod::getVerifyDayNumber(){
    return this -> verify_end_index - this -> verify_start_index + 1;
}

int TradePeriod::getTestDayNumber(){
    return this -> test_end_index - this -> test_start_index + 1;
}


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

void createStock(Stock* stock_list, int size, string **data, int day_number, int start_index, int end_index) {
    for (int j = 0; j < size; j++) {
        stock_list[j].idx = j;
        stock_list[j].init(day_number + 1);
        stock_list[j].company_name = data[0][j+1];
        for (int k = start_index - 1; k <= end_index; k++) {
            stock_list[j].price_list[k - start_index + 1] = atof(data[k+1][j+1].c_str());
            stock_list[j].date_list[k - start_index + 1] = data[k+1][0];
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

void setWindow(string target_date, int range, int day_number, string *date_list, int &start_index, int &end_index){
    for(int j = 0; j < day_number; j++){
        if(target_date == date_list[j]){
            start_index = j;
            end_index = j + range - 1;
            break;
        }
    }
}

void initial(double *b, int size) {
    for (int j = 0; j < size; j++) {
        b[j] = 0.5;
    }
}

void initPortfolio(Portfolio *p, int size, int day_number, Stock *stock_list){
    for(int j = 0; j < PARTICLENUMBER; j++){
        p[j].init(size, day_number, FUNDS, stock_list);
    }
}

void initPortfolio(Portfolio *p){
    for(int j = 0; j < PARTICLENUMBER; j++){
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

void genTestPortfolio(Portfolio &portfolio, Stock* stock_list, string **data, string *myTrainData, int myTrainData_size) {
    for(int j = 0; j < portfolio.size; j++){
        for(int k = 0; k < myTrainData_size; k++){
            if(data[0][j+1] == myTrainData[k]){
                portfolio.data[j] = 1;
                portfolio.stock_id_list[portfolio.stock_number] = j;
                portfolio.stock_number++;
                break;
            }
        }
    }
}

void genTestPortfolio(Portfolio &portfolio, Stock* stock_list, Portfolio &temp_portfolio) {
    for(int j = 0; j < portfolio.size; j++){
        if(j == temp_portfolio.stock_id_list[portfolio.stock_number] && portfolio.stock_number < temp_portfolio.stock_number){
            portfolio.data[j] = 1;
            portfolio.stock_id_list[portfolio.stock_number] = j;
            portfolio.stock_number++;
        }else{
            portfolio.data[j] = 0;
        }
    }
}


void capitalLevel(Portfolio* portfolio_list, int portfolio_number) {
    for (int j = 0; j < portfolio_number; j++) {
        for (int k = 0; k < portfolio_list[j].stock_number; k++) {
            portfolio_list[j].investment_number[k] = portfolio_list[j].getDMoney() / portfolio_list[j].constituent_stocks[portfolio_list[j].stock_id_list[k]].price_list[0];
            portfolio_list[j].remain_fund[k] = portfolio_list[j].getDMoney() - (portfolio_list[j].investment_number[k] * portfolio_list[j].constituent_stocks[portfolio_list[j].stock_id_list[k]].price_list[0]);
        }

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

void capitalLevel(Portfolio* portfolio_list, int portfolio_number, Portfolio &result) {
    for (int j = 0; j < portfolio_number; j++) {
        for (int k = 0; k < portfolio_list[0].stock_number; k++) {
            portfolio_list[0].investment_number[k] = result.investment_number[k];
            portfolio_list[0].remain_fund[k] = result.remain_fund[k];
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
    pWorst.copyP(portfolio_list[PARTICLENUMBER - 1]);
    for (int j = 0; j < PARTICLENUMBER; j++) {
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

void outputFile(Portfolio& portfolio, string file_name) {
    ofstream outfile;
    outfile.open(file_name, ios::out);
    outfile << setprecision(15);

    outfile << "Iteration," << ITERNUMBER << endl;
    outfile << "Element number," << PARTICLENUMBER << endl;
    outfile << "Delta," << DELTA << endl;
    outfile << "Exp times," << EXPNUMBER << endl << endl;
    
    outfile << "Init funds," << portfolio.funds << endl;
    outfile << "Final funds," << portfolio.total_money[portfolio.day_number - 1] << endl;
    outfile << "Real award," << portfolio.getProfit() << endl << endl;
    
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
        outfile << portfolio.date_list[j] << ",";
        for (int k = 0; k < portfolio.stock_number; k++) {
            outfile << (portfolio.constituent_stocks[portfolio.stock_id_list[k]].price_list[j + 1] * portfolio.investment_number[k]) + portfolio.remain_fund[k] << ",";
        }
        outfile << portfolio.total_money[j] << endl;
    }
    outfile << endl;
    outfile.close();
}

void recordCPUTime(double START, double END){
    double total_time = (END - START) / CLOCKS_PER_SEC;
    ofstream outfile_time;
    string file_name = FILE_DIR + "/" + "time.txt";
    outfile_time.open(file_name, ios::out);
    outfile_time << "total time: " << total_time << " sec" << endl;
}

void recordTotalTestResult(vector<double> &total_fs, vector<string> &total_date, ofstream &outfile_LNLP){
    Portfolio portfolio(1, total_fs.size(), FUNDS);
    portfolio.stock_number = 1;
    for(int j = 0; j < total_fs.size(); j++){
        portfolio.total_money[j] = total_fs[j];
        portfolio.date_list[j] = total_date[j];
    }
    countTrend(&portfolio, 1, FUNDS);
    ofstream outfile_total_result;
    string file_name = FILE_DIR + "/" + "total_test_result.csv";
    outfile_total_result.open(file_name, ios::out);
    outfile_total_result << setprecision(15);
    
    outfile_total_result << "Test Date," << total_date[0] << "-" << total_date[total_date.size() - 1] << endl;
    outfile_total_result << "Total days," << portfolio.day_number << endl;
    outfile_total_result << "Iteration," << ITERNUMBER << endl;
    outfile_total_result << "Element number," << PARTICLENUMBER << endl;
    outfile_total_result << "Delta," << DELTA << endl;
    outfile_total_result << "Exp times," << EXPNUMBER << endl << endl;
    
    outfile_total_result << "Init funds," << portfolio.funds << endl;
    outfile_total_result << "Final funds," << portfolio.total_money[portfolio.day_number - 1] << endl;
    outfile_total_result << "Real award," << portfolio.getProfit() << endl << endl;
    
    outfile_total_result << "m," << portfolio.m << endl;
    outfile_total_result << "Daily_risk," << portfolio.daily_risk << endl;
    outfile_total_result << "Trend," << portfolio.trend << endl << endl;
    
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
        
        outfile_total_result << "Quadratic trend line," << portfolio.a << "x^2 + " << portfolio.b << "x + " << FUNDS << endl << endl;
        outfile_total_result << "Quadratic m," << c << endl;
        outfile_total_result << "Quadratic daily risk," << d << endl;
        if(c < 0){
            outfile_total_result << "Quadratic trend," << c * d << endl << endl;
        }else{
            outfile_total_result << "Quadratic trend," << c / d << endl << endl;
        }
    }
    else {
        outfile_total_result << "Quadratic trend line," << portfolio.a << "x^2 + " << portfolio.b << "x + " << FUNDS << endl << endl;
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

        outfile_total_result << "Linear m," << c << endl;
        outfile_total_result << "Linear daily risk," << d << endl;
        if(c < 0){
            outfile_total_result << "Linear trend," << c * d << endl << endl;
        }else{
            outfile_total_result << "Linear trend," << c / d << endl << endl;
        }
    }
    
    outfile_total_result << "Date,FS" << endl;
    for(int j = 0; j < portfolio.day_number; j++){
        outfile_total_result << portfolio.date_list[j] << "," << portfolio.total_money[j] << endl;
    }
    outfile_total_result.close();
    
    outfile_LNLP << "Real Reward," << portfolio.getProfit() << "Trend Ratio," <<portfolio.trend << "Risk," << portfolio.daily_risk;
}

bool isVerifyFinish(TradePeriod &trade_period, double standard_funds){
    bool result = false;
    
    if(STOPTYPE == 0){
        double lowerBound = standard_funds * (1 - LOWER);
        double upperBound = standard_funds * (1 + UPPER);
        double myFunds = trade_period.verify_result.funds + trade_period.verify_result.getProfit();
        if( myFunds < lowerBound || myFunds > upperBound){
            result = true;
        }
    }else if(STOPTYPE == 1){
        if(trade_period.verify_result.m < 0){
            result = true;
        }
    }else if(STOPTYPE == 2){
        double slope = 2 * trade_period.verify_result.a * trade_period.verify_result.day_number + trade_period.verify_result.b;
        if(slope < 0){
            result = true;
        }
    }else if(STOPTYPE == 3){
        double lowerBound = standard_funds * (1 - LOWER);
        double upperBound = standard_funds * (1 + UPPER);
        double myFunds = trade_period.verify_result.funds + trade_period.verify_result.getProfit();
        if( myFunds < lowerBound || myFunds > upperBound || trade_period.verify_result.m < 0){
            result = true;
        }
    }
    
    if(MODE == 1){
        if(trade_period.getTestDayNumber() == TRAINRANGE){
            result = true;
        }
    }
    
    return result;
}

void startTrain(Portfolio &result, Stock *stock_list, int size, int range_day_number, double funds, int particle_number, int exp_number, int iter_number){
    
    double *beta_ = new double[size];
    Portfolio expBest(size, range_day_number, funds, stock_list);
    Portfolio gBest(size, range_day_number, funds, stock_list);
    Portfolio gWorst(size, range_day_number, funds, stock_list);
    Portfolio pBest(size, range_day_number, funds, stock_list);
    Portfolio pWorst(size, range_day_number, funds, stock_list);
    Portfolio* portfolio_list = new Portfolio[particle_number];
    initPortfolio(portfolio_list, size, range_day_number, stock_list);
    
    for(int n = 0; n < exp_number; n++){
        cout << "___" << n << "___" << endl;
        gBest.init();
        gWorst.init();
        gBest.trend = 0;
        gWorst.trend = DBL_MAX;
        initial(beta_, size);
        for(int i = 0; i < iter_number; i++){
            pBest.init();
            pWorst.init();
            initPortfolio(portfolio_list);
            genPortfolio(portfolio_list, stock_list, particle_number, beta_, n, i);
            capitalLevel(portfolio_list, particle_number);
            countTrend(portfolio_list, particle_number, funds);
            recordGAnswer(portfolio_list, gBest, gWorst, pBest, pWorst);
            adjBeta(gBest, pWorst, beta_);
        }
        recordExpAnswer(expBest, gBest);
    }
    
    
    expBest.print();
    result.copyP(expBest);
    delete[] portfolio_list;
    delete[] beta_;

}

void startVerify(Portfolio &result, Portfolio &train_result, Stock *stock_list, int size, int range_day_number, double funds){
    
    Portfolio* portfolio = new Portfolio[1];
    portfolio[0].init(size, range_day_number, funds, stock_list);
    genTestPortfolio(portfolio[0], stock_list, train_result);
    
    if(MODE == 2){
        capitalLevel(portfolio, 1, train_result);
    }else{
        capitalLevel(portfolio, 1);
    }
    countTrend(portfolio, 1, funds);
    portfolio[0].countQuadraticYLine();
    result.copyP(portfolio[0]);
    delete[] portfolio;
}

void startTest(Portfolio &result, Portfolio &train_result, Stock *stock_list, int size, int range_day_number, double funds){
    Portfolio* portfolio = new Portfolio[1];
    portfolio[0].init(size, range_day_number, funds, stock_list);
    genTestPortfolio(portfolio[0], stock_list, train_result);
    capitalLevel(portfolio, 1);
    countTrend(portfolio, 1, funds);
    portfolio[0].countQuadraticYLine();
    result.copyP(portfolio[0]);
    delete[] portfolio;
}

int main(int argc, const char * argv[]) {
    
    ofstream outfile_LNLP;
    string total_data_name = "total_LN+L&P_30-30.csv";
    outfile_LNLP.open(total_data_name, ios::out);
    outfile_LNLP << setprecision(15);
    
    for(int lo = 30; lo > 0; lo--){
        for(int up = 30; up > 0; up--){
            outfile_LNLP << "Loss," << lo << "Profit," << up;
            LOWER = double(lo) / 100;
            UPPER = double(up) / 100;
            string temp_file_dir = FILE_DIR;
            FILE_DIR += "_" + to_string(lo) + "%" + to_string(up) + "%";
    double START, END;
    START = clock();
    srand(114);
    cout << setprecision(10);
    int size;
    int day_number;
    string target_date = STARTDATE;
    double funds = FUNDS;
    double test_funds = FUNDS;
    double capital_highest_point = 0;
    double MDD = 0;
    double pos_award = 0;
    double neg_award = 0;
    bool isLastDay = false;
    vector<double> total_fs;
    vector<string> total_date;
    string** data;
    vector<vector<string>> data_vector;
    Portfolio result;
    
    createDir(FILE_DIR);
    readData("2011-2020.csv", data_vector, size, day_number);
    data = vectorToArray(data_vector);
    TradePeriod trade_period;
    trade_period.init(day_number, data);
    
    ofstream outfile_total_data;
    string total_data_name = FILE_DIR + "/" + "total_data.csv";
    outfile_total_data.open(total_data_name, ios::out);
    outfile_total_data << "train date,test date,test day number,MDD,profit,PF" << endl;
    
    while(!isLastDay){
        
        setWindow(target_date, TRAINRANGE, day_number, trade_period.date_list, trade_period.train_start_index, trade_period.train_end_index);
        
        //______Train______
        cout << LOWER << "-" << UPPER << "_" << MODE << "_" << trade_period.getTrainStartDate() << " - " << trade_period.getTrainEndDate() << endl;
        Stock* stock_list = new Stock[size];
        createStock(stock_list, size, data, trade_period.getTrainDayNumber(), trade_period.train_start_index, trade_period.train_end_index);
        
        result.init(size, trade_period.getTrainDayNumber(), FUNDS, stock_list);
        
        startTrain(result, stock_list, size, trade_period.getTrainDayNumber(), FUNDS, PARTICLENUMBER, EXPNUMBER, ITERNUMBER);
        
        trade_period.train_result.copyP(result);
        if(result.trend != 0){
            outputFile(trade_period.train_result, getOutputFilePath(trade_period.getTrainStartDate(), trade_period.getTrainEndDate(), FILE_DIR, "train"));
        }else {
            target_date = trade_period.getDate(trade_period.train_start_index + 1);
            continue;
        }
        
        //______Verify______
        cout << "______Verify______" << endl << endl;
        
        double standard_funds;
        double verify_funds;
        
        standard_funds = result.funds + result.getProfit();
        if(MODE == 2){
            verify_funds = result.total_money[0];
        }else{
            verify_funds = FUNDS;
        }
        trade_period.verify_start_index = trade_period.train_start_index;
        trade_period.verify_end_index = trade_period.train_end_index;
        trade_period.test_start_index = trade_period.train_end_index + 1;
        while(true){
            if(MODE == 2){
                trade_period.verify_start_index++;
            }
            trade_period.verify_end_index++;
            trade_period.test_end_index = trade_period.verify_end_index;
            
            if(trade_period.verify_end_index == day_number-1){
                isLastDay = true;
            }
            
            stock_list = new Stock[size];
            createStock(stock_list, size, data, trade_period.getVerifyDayNumber(), trade_period.verify_start_index, trade_period.verify_end_index);
            
            result.init(size, trade_period.getVerifyDayNumber(), verify_funds, stock_list);
            startVerify(result, trade_period.train_result, stock_list, size, trade_period.getVerifyDayNumber(), verify_funds);
            trade_period.verify_result.copyP(result);
            if(MODE == 2){
                verify_funds = result.total_money[0];
            }
            delete[] stock_list;
            
            if(isVerifyFinish(trade_period, standard_funds) || isLastDay){
                outputFile(trade_period.verify_result, getOutputFilePath(trade_period.getVerifyStartDate(), trade_period.getVerifyEndDate(), FILE_DIR, "verify"));
                break;
            }
            
        }
        
        //______Test______
        
        cout << "______Test______" << endl << endl;
        stock_list = new Stock[size];
        createStock(stock_list, size, data, trade_period.getTestDayNumber(), trade_period.test_start_index, trade_period.test_end_index);
        result.init(size, trade_period.getTestDayNumber(), test_funds, stock_list);
        startTest(result, trade_period.train_result, stock_list, size, trade_period.getTestDayNumber(), test_funds);
        trade_period.test_result.copyP(result);
        
        for(int j = 0; j < trade_period.getTestDayNumber(); j++){
            total_fs.push_back(trade_period.test_result.total_money[j]);
            total_date.push_back(trade_period.test_result.date_list[j]);
        }
        delete[] stock_list;
        test_funds = trade_period.test_result.funds + trade_period.test_result.getProfit();
        
        if(trade_period.test_result.getProfit() >= 0){
            pos_award += trade_period.test_result.getProfit();
        }else{
            neg_award += -1 * trade_period.test_result.getProfit();
        }
        
        if(neg_award == 0){
            trade_period.test_result.PF = -1;
        }else{
            trade_period.test_result.PF = pos_award / neg_award;
        }
        
        outfile_total_data << trade_period.getTrainStartDate() << " - " << trade_period.getTrainEndDate() << ",";
        outfile_total_data << trade_period.getTestStartDate() << " - " << trade_period.getTestStartDate() << ",";
        outfile_total_data << trade_period.getTestDayNumber() << ",";
        outfile_total_data << trade_period.test_result.MDD << ",";
        outfile_total_data << trade_period.test_result.getProfit() << ",";
        outfile_total_data << trade_period.test_result.PF << endl;
        
        outputFile(trade_period.test_result, getOutputFilePath(trade_period.getTestStartDate(), trade_period.getTestEndDate(), FILE_DIR, "test"));
        
        cout << trade_period.getTestStartDate() << " - " << trade_period.getTestEndDate() << endl;
        cout << "test days: " << trade_period.getTestDayNumber() << endl;
        cout << endl << endl;

        target_date = trade_period.getDate(trade_period.test_end_index - TRAINRANGE + 1);
    }
    
    recordTotalTestResult(total_fs, total_date, outfile_LNLP);
    END = clock();
    recordCPUTime(START, END);
    
            
    
    for(int j = 0; j < data_vector.size(); j++){
        data_vector[j].clear();
        delete[] data[j];
    }
    
    data_vector.clear();
    delete[] data;
        total_fs.clear();
        total_date.clear();
            FILE_DIR = temp_file_dir;
        }
    }
    
    return 0;
}
