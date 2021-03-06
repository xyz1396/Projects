# include <iostream>
using namespace std;

bool duplicate(int numbers[], int length, int* duplication)
{
    if (numbers==nullptr || length<=0)
    {
        return false;
    }

    for (int i=0; i<length;i++)
    {
        if(numbers[i]<0 || numbers[i]>length -1)
            return false;
    }

    for (int i=0;i<length;++i)
    {
        while(numbers[i]!=i)
        {
            if(numbers[i]==numbers[numbers[i]])
            {
                *duplication=numbers[i];
                return true;
            }

            // swap numbers[i] and numbers[numbers[i]]
            int temp=numbers[i];
            numbers[i]=numbers[temp];
            numbers[temp]=temp;
        }
    }
}

int main()
{
    int array[7]={3,1,2,0,2,5,3};
    bool foundDup = duplicate(array,7,&array[0]);
    cout<< foundDup << endl;
}