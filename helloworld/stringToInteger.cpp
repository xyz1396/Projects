# include <iostream>

using namespace std;

int StrToInt(char* string){
    int number=0;
    while(*string !=0){
        // ASC码48就是'0'，也就是说'0'的值是48，而后依次是'1'到'9'。 
        // 这样正好是char型减去48就是它对应的int值
        number=number*10 + *string -'0';
        printf("%d \n",*string);
        ++string;
    }
    return number;
}

int main(){
    char number[20];
    cout << "Input number: ";
    cin >> number;
    cout << "input number: " << number << endl;
    int intNumber = StrToInt(&number[0]);
    cout << "int number: " << intNumber << endl;
    cout << "int number: " << "2"-"0" << endl;
}

