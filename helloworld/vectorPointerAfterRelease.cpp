#include <iostream>
#include <vector>
using namespace std;

int main(int argc, char const *argv[])
{
    vector<double> a = {1, 2, 3, 4, 5, 6};
    double *b=new double(1.5);
    *b=1.2;
    cout << *(b + 3) << endl;
    // a = vector<double>();
    // for (int i = 0; i < 6; i++)
    // {
    //     cout << sizeof(b) << endl;
    // }
    delete[] b;
    // cout << b + 3 << endl;
    // cout << *(b + 1) << endl;
    return 0;
}
