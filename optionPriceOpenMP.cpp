#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <omp.h>
#include <ctime>

using namespace std;

// Outputs the norm PDF
double norm_pdf(const double& x)
{
    return (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x);
}

// Recursive norm CDF approximator
double norm_cdf(const double& x)
{
    double k = 1.0/(1.0 + 0.2316419*x);
    double k_sum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));

    if (x >= 0.0)
    {
        return (1.0 - (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x) * k_sum);
    }
    else
    {
        return 1.0 - norm_cdf(-x);
    }
}

// Calculator for d1 and d2 inputs into Norm distribution terms in formula
double d_i(const int& i, const double& S, const double& K, const double& r, const double& v, const double& T)
{
    return (log(S/K) + (r + (pow(-1,i-1))*0.5*v*v)*T)/(v*(pow(T,0.5)));
}


// Calculator for European standard call price given standard formula inputs
double call_price(const double& S, const double& K, const double& r, const double& v, const double& T)
{
    return S * norm_cdf(d_i(1, S, K, r, v, T))-K*exp(-r*T) * norm_cdf(d_i(2, S, K, r, v, T));
}

// Calculator for European standard put price given standard formula inputs
double put_price(const double& S, const double& K, const double& r, const double& v, const double& T)
{
    return -S*norm_cdf(-d_i(1, S, K, r, v, T))+K*exp(-r*T) * norm_cdf(-d_i(2, S, K, r, v, T));
}


vector< vector<double> > read_csv()
{
    ifstream in("variables.csv");

    string line, field;

    vector< vector<double> > vars;  // the 2D array
    vector<double> arr;                // array of values for one line only

    while ( getline(in,line) )    // get next line in file
    {
        arr.clear();
        stringstream ss(line);

        while (getline(ss,field,','))  // break line into comma delimitted fields
        {
            arr.push_back(stod(field));  // add each field to the 1D array
        }

        vars.push_back(arr);  // add the 1D array to the 2D array
    }

    return vars;
}

double find_call_min(vector< vector<double> > vars, int rows)
{
    double minCall = 99999.0;
    double call;
    #pragma omp parallel for reduction(min: minCall)
    for (int i = 0; i < rows; i++)
    {
        call = call_price(vars[i][0], // Spot price
                          vars[i][1], // Strike price
                          vars[i][2], // Risk-free rate
                          vars[i][3], // Volatility
                          vars[i][4]);// Time to maturity


        if (minCall > call)
            minCall = call;

    }
    return minCall;
}

double find_call_max(vector< vector<double> > vars, int rows)
{
    double maxCall = 0.0;
    double call;
    #pragma omp parallel for reduction(max: maxCall)
    for (int i = 0; i < rows; i++)
    {
        call = call_price(vars[i][0], // Spot price
                          vars[i][1], // Strike price
                          vars[i][2], // Risk-free rate
                          vars[i][3], // Volatility
                          vars[i][4]);// Time to maturity
        if (maxCall < call)
            maxCall = call;

    }
    return maxCall;
}

double find_put_min(vector< vector<double> > vars, int rows)
{
    double minPut = 99999.0;
    double put;
    #pragma omp parallel for reduction(min: minPut)
    for (int i = 0; i < rows; i++)
    {
        put = put_price(vars[i][0], // Spot price
                        vars[i][1], // Strike price
                        vars[i][2], // Risk-free rate
                        vars[i][3], // Volatility
                        vars[i][4]);// Time to maturity

        if (minPut > put)
            minPut = put;

    }
    return minPut;
}

double find_put_max(vector< vector<double> > vars, int rows)
{
    double maxPut = 0.0;
    double put;
    #pragma omp parallel for reduction(max: maxPut)
    for (int i = 0; i < rows; i++)
    {
        put = put_price(vars[i][0], // Spot price
                        vars[i][1], // Strike price
                        vars[i][2], // Risk-free rate
                        vars[i][3], // Volatility
                        vars[i][4]);// Time to maturity

        if (maxPut < put)
            maxPut = put;

    }
    return maxPut;
}


int main(int argc, char **argv)
{
    #pragma openmp parallel
    int start_s=clock();

    vector< vector<double> > vars = read_csv(); // read S,K,r,v,T variables from csv

    int rows, cols;
    rows = 10000;

    double minCall, minPut, maxCall, maxPut;

    minCall, maxCall = find_call_min(vars,rows),find_call_max(vars,rows);
    minPut, maxPut = find_put_min(vars,rows), find_put_max(vars,rows);

    cout << "Min call price found:        " << minCall << endl;
    cout << "Max call price found:        " << maxCall << endl;
    cout << "Min put price found:         " << minPut << endl;
    cout << "Max put price found:         " << maxPut << endl;

    int stop_s=clock();
    cout << "Time taken:                  " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " seconds." << endl;

    return 0;
}
