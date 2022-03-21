
#include <string>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;

class Stock {
public:
    string company_name;
    int idx;
    int day_number = 0;
    double* price_list = NULL;
    string* date_list = NULL;
    void init(int);
    void init(Stock&);
    ~Stock();
};

Stock::~Stock(){
    delete[] price_list;
    delete[] date_list;
    this -> price_list = NULL;
    this -> date_list = NULL;
}

void Stock::init(int day_number) {
    this -> day_number = day_number;
    if(this -> price_list != NULL){
        delete[] this -> price_list;
        delete[] this -> date_list;
    }
    this -> price_list = new double[day_number];
    this -> date_list = new string[day_number];
}

void Stock::init(Stock &s){
    this -> init(s.day_number);
    this -> day_number = s.day_number;
    this -> idx = s.idx;
    this -> company_name = s.company_name;
    for(int j = 0; j < this -> day_number; j++){
        this -> price_list[j] = s.price_list[j];
    }
}

class Portfolio {
public:
    int gen = 0;
    int exp = 0;
    int answer_counter = 0;
    int stock_number = 0;
    int size = 0;
    int day_number = 0;
    
    
    double m = 0;
    double a = -1;
    double b = -1;
    double daily_risk = 0;
    double trend = 0;
    double emotion_index = 0;
    double ei_m = 0;
    double ei_daily_risk = 0;
    double remain_money = 0;
    double funds = 0;
    double capital_highest_point = 0;
    double MDD = 0;
    double PF = 0;
    double DPF = 0;
    
    string date = "";
    
    int* data = NULL;
    int* investment_number = NULL;
    int* stock_id_list = NULL;
    double* total_money = NULL;
    double* remain_fund = NULL;
    Stock* constituent_stocks = NULL;
    string* date_list = NULL;
    
    void constr();
    void constr(int, int, double);
    void constr(int, int, double, Stock*);
    void init();
    void init(int, int, double, Stock*);
    void init(int, int, double, Stock*, double, double);
    int getDMoney();
    double getRemainMoney();
    double getNormalY(int);
    double getProfit();
    void countQuadraticYLine();
    double getQuadraticY(int);
    void copyP(Portfolio&);
    void print();
    Portfolio();
    Portfolio(int, int, double);
    Portfolio(int, int, double, Stock*);
    ~Portfolio();
};

Portfolio::Portfolio(){
    constr();
};

Portfolio::Portfolio(int size, int day_number, double funds){
    constr(size, day_number, funds);
}

Portfolio::Portfolio(int size, int day_number, double funds, Stock *stock_list) {

    constr(size, day_number, funds, stock_list);

}

Portfolio::~Portfolio() {
    if(this -> data != NULL){
        delete[] this -> data;
        delete[] this -> investment_number;
        delete[] this -> constituent_stocks;
        delete[] this -> remain_fund;
        delete[] this -> total_money;
        delete[] this -> stock_id_list;
        delete[] this -> date_list;
    }
    
    this -> data = NULL;
    this -> investment_number = NULL;
    this -> constituent_stocks = NULL;
    this -> remain_fund = NULL;
    this -> total_money = NULL;
    this -> stock_id_list = NULL;
    this -> date_list = NULL;
}

void Portfolio::constr() {
    this->gen = 0;
    this->exp = 0;
    this->answer_counter = 0;
    this->stock_number = 0;
    this->m = 0;
    this->a = -1;
    this->b = -1;
    this->daily_risk = 0;
    this->trend = 0;
    this->stock_number = 0;
    this->capital_highest_point = 0;
    this->MDD = 0;
    this->date = "";
    this->emotion_index = 0;
    this->ei_m = 0;
    this->ei_daily_risk = 0;
}

void Portfolio::constr(int size, int day_number, double funds) {
    
    constr();
    this->remain_money = funds;
    this->funds = funds;
    this->size = size;
    this->day_number = day_number;

    this->data = new int[size];
    this->investment_number = new int[size];
    this->remain_fund = new double[size];
    this->total_money = new double[day_number];
    this->stock_id_list = new int[this->size];
    this->date_list = new string[day_number];

    for (int j = 0; j < size; j++) {
        data[j] = 0;
    }
}

void Portfolio::constr(int size, int day_number, double funds, Stock* stock_list) {

    constr(size, day_number, funds);
    this->constituent_stocks = new Stock[size];
    for (int j = 0; j < size; j++) {
        this->constituent_stocks[j].init(stock_list[j]);
    }

    for (int j = 0; j < day_number; j++) {
        this->date_list[j] = stock_list[0].date_list[j + 1];
    }
}

void Portfolio::init() {

    constr();
    this -> remain_money = this -> funds;
    
    if(this -> data != NULL){
        delete[] this -> data;
        delete[] this -> investment_number;
        delete[] this -> remain_fund;
        delete[] this -> total_money;
        delete[] this -> stock_id_list;
    }
    this -> data = new int[this -> size];
    for (int j = 0; j < size; j++) {
        data[j] = 0;
    }
    this -> investment_number = new int[this ->size];
    this -> remain_fund = new double[this -> size];
    this -> total_money = new double[this -> day_number];
    this -> stock_id_list = new int[this -> size];
}

void Portfolio::init(int size, int day_number, double funds, Stock* stock_list) {
    
    constr();
    this -> remain_money = funds;
    this -> funds = funds;
    this -> size = size;
    this -> day_number = day_number;
    
    if (this->data != NULL) {
        delete[] this->data;
        delete[] this->investment_number;
        delete[] this->constituent_stocks;
        delete[] this->remain_fund;
        delete[] this->total_money;
        delete[] this->stock_id_list;
        delete[] this->date_list;
    }
    
    this->data = new int[size];
    this->investment_number = new int[size];
    this->remain_fund = new double[size];
    this->total_money = new double[day_number];
    this->stock_id_list = new int[size];
    this->constituent_stocks = new Stock[size];
    this->date_list = new string[day_number];

    for (int j = 0; j < size; j++) {
        data[j] = 0;
    }
    
    for(int j = 0; j < size; j++){
        this -> constituent_stocks[j].init(stock_list[j]);
    }
    
    for(int j = 0; j < day_number; j++){
        this -> date_list[j] = stock_list[0].date_list[j+1];
    }
}

void Portfolio::init(int size, int day_number, double funds, Stock* stock_list, double capital_highest_point, double MDD) {
    init(size, day_number, funds, stock_list);
    this -> capital_highest_point = capital_highest_point;
    this -> MDD = MDD;
}

int Portfolio::getDMoney() {
    if (this->stock_number != 0) {
        int temp = this -> funds;
        return temp / this->stock_number;
    }
    else {
        return this -> funds;
    }
}

double Portfolio::getRemainMoney() {
    if (this->stock_number != 0) {
        double sum = 0;
        for (int j = 0; j < stock_number; j++) {
            sum += this -> remain_fund[j];
        }
        int temp = this -> funds;
        return (temp % this->stock_number) + sum + (this -> funds - temp);
    }
    else {
        return this -> funds;
    }
}

double Portfolio::getNormalY(int day) {
    return this->m * day + this -> funds;
}

double Portfolio::getProfit(){
    return this -> total_money[this -> day_number - 1] - this -> funds;
}

void Portfolio::countQuadraticYLine() {

    MatrixXd A(this -> day_number, 2);
    VectorXd Y(this -> day_number, 1);
    for (int j = 0; j < this -> day_number; j++) {
        for (int k = 0; k < 2; k++) {
            A(j, k) = pow(j + 1, 2 - k);
        }
        Y(j, 0) = this -> total_money[j] - this -> funds;
    }
    Vector2d theta = A.colPivHouseholderQr().solve(Y);
    this -> a = theta(0, 0);
    this -> b = theta(1, 0);
}


double Portfolio::getQuadraticY(int day) {
    return this->a * pow(day, 2) + this->b * day + this -> funds;
}

void Portfolio::copyP(Portfolio& a) {
    this -> gen = a.gen;
    this -> exp = a.exp;
    this -> answer_counter = a.answer_counter;
    this -> stock_number = a.stock_number;
    this -> m = a.m;
    this -> a = a.a;
    this -> b = a.b;
    this -> daily_risk = a.daily_risk;
    this -> trend = a.trend;
    this -> stock_number = a.stock_number;
    this -> remain_money = a.remain_money;
    this -> funds = a.funds;
    this -> MDD = a.MDD;
    this -> capital_highest_point = a.capital_highest_point;
    this -> PF = a.PF;
    this -> DPF = a.DPF;
    this -> emotion_index = a.emotion_index;
    this->ei_m = a.ei_m;
    this->ei_daily_risk = a.ei_daily_risk;
    
    if(this -> size != a.size || this -> day_number != a.day_number){
        this -> size = a.size;
        this -> day_number = a.day_number;
        if(this -> data != NULL){
            delete[] this -> data;
            delete[] this -> investment_number;
            delete[] this -> constituent_stocks;
            delete[] this -> remain_fund;
            delete[] this -> total_money;
            delete[] this -> stock_id_list;
            delete[] this -> date_list;
        }
        
        this -> data = new int[size];
        this -> investment_number = new int[size];
        this -> remain_fund = new double[size];
        this -> total_money = new double[day_number];
        this -> stock_id_list = new int[size];
        this -> constituent_stocks = new Stock[size];
        this -> date_list = new string[day_number];
        
        for(int j = 0; j < size; j++){
            this -> constituent_stocks[j].init(a.constituent_stocks[j]);
        }
    }
    
    if(this -> constituent_stocks[0].date_list[0] != a.constituent_stocks[0].date_list[0]){
        for(int j = 0; j < size; j++){
            this -> constituent_stocks[j].init(a.constituent_stocks[j]);
        }
    }
    
    for (int j = 0; j < a.size; j++) {
        this -> data[j] = a.data[j];
    }
    for (int j = 0; j < a.stock_number; j++) {
        this -> investment_number[j] = a.investment_number[j];
        this -> remain_fund[j] = a.remain_fund[j];
        this -> stock_id_list[j] = a.stock_id_list[j];
    }
    for (int j = 0; j < a.day_number; j++) {
        this -> total_money[j] = a.total_money[j];
        this -> date_list[j] = a.date_list[j];
    }
}

void Portfolio::print(){
    cout << "exp: " << this -> exp << endl;
    cout << "gen: " << this -> gen << endl;
    cout << "data: " << endl;
    for(int j = 0; j < this -> size; j++){
        if(j > 0){
            cout << ", ";
        }
        cout << this -> data[j];
    }
    cout << endl;
    cout << "Stock number: " << this -> stock_number << endl;
    cout << "Stock ID: " << endl;
    for(int j = 0; j < this -> stock_number; j++){
        if(j > 0){
            cout << ", ";
        }
        cout << this -> stock_id_list[j];
    }
    cout << endl;
    cout << "Trend: " << this -> trend << endl;
    cout << "Real award: " << this -> getProfit() << endl;
    cout << endl;
}
