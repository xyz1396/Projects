#include <vector>
#include <iostream>
using namespace std;

int main(int argc, char const *argv[])
{
    vector<double> a;
    a = {0, 1.999, 2, .333, 4.33, 5.99};
    cout << "a: ";
    for (double i : a)
    {
        cout << i << " ";
    }
    cout << endl;

    // vector to array

    // shallow copy
    double *b;
    b = &a[0];
    cout << "b: ";
    for (int i = 0; i < 6; i++)
    {
        cout << b[i] << " ";
    }
    cout << endl;

    // deep copy
    b = new double[6];
    copy(a.begin(), a.end(), b);
    cout << "data b: ";
    for (int i = 0; i < 6; i++)
    {
        cout << b[i] << " ";
    }
    cout << endl;

    // copy array
    double *d = new double[3];
    copy(b + 2, b + 4, d);
    cout << "\n"
         << "d: " << endl;
    for (int i = 0; i < 2; i++)
    {
        cout << d[i] << " ";
    }
    cout << endl;

    // array to vector
    vector<double> c = vector<double>(b, b + 6);
    c.push_back(7.77);
    cout << "\n"
         << "c: ";
    for (double i : c)
    {
        cout << i << " ";
    }
    cout << endl;
    b[1] = 2.33;
    for (double i : c)
    {
        cout << i << " ";
    }
    cout << endl;

    // new array
    b = new double[6];
    b[0] = 0.55;
    cout << "new b: ";
    for (int i = 0; i < 6; i++)
    {
        cout << b[i] << " ";
    }
    cout << endl;
    delete[] b;
    b = nullptr;
    delete[] b;
    b = new double[5];
    b[1] = 2.33;
    b = nullptr;
    b = &a[0];
    cout << "new b: ";
    for (int i = 0; i < 6; i++)
    {
        cout << b[i] << " ";
    }
    cout << endl;
    return 0;
}
