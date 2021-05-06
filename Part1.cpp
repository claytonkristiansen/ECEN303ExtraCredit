#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>      //RAND_MAX
#include <stdlib.h>     //srand, rand
#include <time.h>       //time
#include <math.h>       //exp

using std::cout;
using std::cin;
using std::vector;

double Rand0to1()
{
    return ((double)rand()) / RAND_MAX;
}

class ExpRV
{
    double m_lambda = 1;
public:
    ExpRV(double lambda)
    {
        m_lambda = lambda;
    }

    double Sample()                                 //1
    {
        return -(log(1 - Rand0to1())) / m_lambda;
    } 
};

class PoissonRV
{
    double m_lambda = 1;
public:
    PoissonRV(double lambda)
    {
        m_lambda = lambda;
    }

    int Sample()                                 //1
    {
        int m = 0;
        ExpRV exponential(1);
        double sum = exponential.Sample();
        while(sum <= m_lambda)
        {
            ++m;
            sum += exponential.Sample();
        }
        return m;
    } 
};

class TimeRV
{
public:
    double Sample()
    {
        PoissonRV N(5);
        PoissonRV M(2);
        ExpRV X(0.5);
        double totalTime = 0;
        int numTimesToBank = (int)round(N.Sample());
        for(int i = 0; i < numTimesToBank; ++i)
        {
            int numCustomers = (int)round(M.Sample());
            for(int k = 0; k < numCustomers; ++k)
            {
                totalTime += X.Sample();
            }
        }
        return totalTime;
    }
    double AtLeast5TimesSample()
    {
        PoissonRV N(5);
        PoissonRV M(2);
        ExpRV X(0.5);
        double totalTime = 0;
        int numTimesToBank = (int)round(N.Sample());
        while (numTimesToBank < 5) numTimesToBank = (int)round(N.Sample());
        for(int i = 0; i < numTimesToBank; ++i)
        {
            int numCustomers = (int)round(M.Sample());
            for(int k = 0; k < numCustomers; ++k)
            {
                totalTime += X.Sample();
            }
        }
        return totalTime;
    }
    double AtLeast1CustomerSample()
    {
        PoissonRV N(5);
        PoissonRV M(2);
        ExpRV X(0.5);
        double totalTime = 0;
        int numTimesToBank = (int)round(N.Sample());
        for(int i = 0; i < numTimesToBank; ++i)
        {
            int numCustomers = (int)round(M.Sample());
            while (numCustomers < 1) numCustomers = (int)round(M.Sample());
            for(int k = 0; k < numCustomers; ++k)
            {
                totalTime += X.Sample();
            }
        }
        return totalTime;
    }
    double AtLeast5TimesAnd1CustomerSample()
    {
        PoissonRV N(5);
        PoissonRV M(2);
        ExpRV X(0.5);
        double totalTime = 0;
        int numTimesToBank = (int)round(N.Sample());
        while (numTimesToBank < 5) numTimesToBank = (int)round(N.Sample());
        for(int i = 0; i < numTimesToBank; ++i)
        {
            int numCustomers = (int)round(M.Sample());
            while (numCustomers < 1) numCustomers = (int)round(M.Sample());
            for(int k = 0; k < numCustomers; ++k)
            {
                totalTime += X.Sample();
            }
        }
        return totalTime;
    }

    vector<double> NAtLeast5TimesSample(int n)                //2
    {
        vector<double> samples;
        for(int i = 0; i < n; ++i)
        {
            double sample = AtLeast5TimesSample();
            samples.push_back(sample);
            //cout << sample << "\n";
        }
        return samples;
    }

    vector<double> NAtLeast1CustomerSample(int n)                //2
    {
        vector<double> samples;
        for(int i = 0; i < n; ++i)
        {
            double sample = AtLeast1CustomerSample();
            samples.push_back(sample);
            //cout << sample << "\n";
        }
        return samples;
    }

    vector<double> NAtLeast5TimesAnd1CustomerSample(int n)                //2
    {
        vector<double> samples;
        for(int i = 0; i < n; ++i)
        {
            double sample = AtLeast5TimesAnd1CustomerSample();
            samples.push_back(sample);
            //cout << sample << "\n";
        }
        return samples;
    }
};

template<class T>
vector<double> NSamples(int n, T rv)                //2
{
    vector<double> samples;
    for(int i = 0; i < n; ++i)
    {
        double sample = rv.Sample();
        samples.push_back(sample);
        //cout << sample << "\n";
    }
    return samples;
}

template<class T>
double CDF(int n, double x, T rv)                   //3
{
    int sum = 0;
    for(int i = 0; i < n; ++i)
    {
        double sample = rv.Sample();
        if(sample < x)
        {
            ++sum;
        }
    }
    return ((double)sum) / ((double)n);
}

template<class T>
bool VectorToCSV(vector<T> vec, std::ofstream& oFile)
{
    if(!oFile.is_open())
    {
        return false;
    }
    for(T element : vec)
    {
        oFile << element << "\n";
    }
    return true;
}

template<class T>
double Mean(int n, T rv)                            //5
{
    double sum = 0;
    for(int i = 0; i < n; ++i)
    {
        double sample = rv.Sample();
        sum += sample;
    }
    return sum / (double)n;
}

double Mean(vector<double> vec)                            //5
{
    double sum = 0;
    for(double sample : vec)
    {
        sum += sample;
    }
    return sum / (double)vec.size();
}

template<class T>
double Variance(int n, T rv)                        //5
{
    double mean = Mean(n, rv);
    double sum = 0;
    for(int i = 0; i < n; ++i)
    {
        double sample = rv.Sample();
        sum += ((mean - sample) * (mean - sample));
    }
    return sum / ((double)n - 1);
}

double Variance(vector<double> vec)                        //5
{
    double mean = Mean(vec);
    double sum = 0;
    for(double sample : vec)
    {
        sum += ((mean - sample) * (mean - sample));
    }
    return sum / ((double)vec.size() - 1);
}


//-----------------------------------//Part III//-------------------------------------//

double Q9_1(vector<double> samples)
{
    int count = 0;
    for(double sample : samples)
    {
        if(sample <= 20)
        {
            ++count;
        }
    }
    return (double)count / (double)samples.size();    
}

double Q9_2(vector<double> samples)
{
    int count = 0;
    for(double sample : samples)
    {
        if(sample > 20)
        {
            ++count;
        }
    }
    return (double)count / (double)samples.size();    
}

double Q9_3(vector<double> samples)
{
    int count = 0;
    for(double sample : samples)
    {
        if(sample > 20)
        {
            ++count;
        }
    }
    return (double)count / (double)samples.size();    
}

double Q9_4(vector<double> samples)
{
    int count = 0;
    for(double sample : samples)
    {
        if(sample > 20)
        {
            ++count;
        }
    }
    return (double)count / (double)samples.size();    
}

int main()
{
    srand(time(NULL));
    ExpRV erv(0.5);

    vector<double> xVals;
    vector<double> yVals100;
    vector<double> yVals5000;
    std::ofstream xOut("xOut.csv");
    std::ofstream yVals100Out("yVals100.csv");
    std::ofstream yVals5000Out("yVals5000.csv");
    std::ofstream milltionTimes("MillionTimes.csv");

    for(double i = 0; i <= 15; i = i + 0.01)
    {
        xVals.push_back(i);
        yVals100.push_back(CDF(100, i, erv));
        yVals5000.push_back(CDF(5000, i, erv));
    }

    
    VectorToCSV(xVals, xOut);
    VectorToCSV(yVals100, yVals100Out);
    VectorToCSV(yVals5000, yVals5000Out);


    NSamples(10, erv);
    cout << "5a) \tEstimated Mean of X with n=100: \t" << Mean(100, erv) << "\n";
    cout << "5b) \tEstimated Variance of X with n=100: \t" << Variance(100, erv) << "\n";
    cout << "5c) \tEstimated Mean of X with n=5000: \t" << Mean(5000, erv) << "\n";
    cout << "5d) \tEstimated Variance of X with n=5000: \t" << Variance(5000, erv) << "\n";
    cout << "5e) \tCalculated Mean of X: \t\t\t" << 1/0.5 << "\n";
    cout << "5f) \tCalculated Variance of X: \t\t" << 1/(0.5*0.5) << "\n\n";

    PoissonRV ps(5);
    vector<double> samples = NSamples(5000, ps);
    cout << "7a) \tMean of M: \t" << Mean(samples) << "\n";
    cout << "7b) \tVariance of M: \t" << Variance(samples) << "\n\n";

    TimeRV tr;
    vector<double> samples1 = NSamples(1000000, tr);
    vector<double> samples2 = tr.NAtLeast5TimesSample(1000000);
    vector<double> samples3 = tr.NAtLeast1CustomerSample(1000000);
    vector<double> samples4 = tr.NAtLeast5TimesAnd1CustomerSample(1000000);

    cout << "9-i) \tPr: \t\t\t" << Q9_1(samples1) << "\n";
    cout << "9-ii) \tPr: \t\t\t" << Q9_2(samples2) << "\n";
    cout << "9-iii) \tPr: \t\t\t" << Q9_3(samples3) << "\n";
    cout << "9-iv) \tPr: \t\t\t" << Q9_4(samples4) << "\n";

    cout << "10a) \tMean of T: \t\t" << Mean(samples1) << "\n";
    cout << "10b) \tVariance of T: \t\t" << Variance(samples1) << "\n\n";


    return 0;
}
