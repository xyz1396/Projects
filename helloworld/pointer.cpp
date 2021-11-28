# include <iostream>

using namespace std;

class Employee{
    public:
        int id;
        string name;
        float salary;
        void insert(int i, string n, float s){
            id=i;
            name=n;
            salary=s;
        };
        void display(){
            cout << id << "" << name << ""<< salary << endl;
        }
};

int main(){
    int number=30;
    int *p;
    p=&number;
    cout<<"Adreess of number:"<<&number<<endl;
    cout<<"Adress of pointer p:"<<p<<endl;
    cout<<"Value of pinter p:"<<*p<<endl;
    Employee e1;
    e1.insert(201,"tom",250);
    e1.display();
}