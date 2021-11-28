# include <iostream>
# include <vector>

using namespace std;

// typedef struct listNode
// {
//     int val;
//     struct listNode* next;
//     listNode(int x):
//         val(x),next(NULL){}    
// };

struct listNode
{
    int val;
    struct listNode* next;
    listNode(int x):
        val(x),next(NULL){}
};

listNode* findKthToTail(listNode* pListHead, unsigned int k)
{
    if(pListHead == nullptr)
        return nullptr;
    listNode* pAhead=pListHead;
    listNode* pBehind=nullptr;
    for(unsigned int i =0; i<k-1; ++i)
    {
        pAhead=pAhead->next;
    }

    pBehind=pListHead;

    while(pAhead->next!=nullptr)
    {
        pAhead=pAhead->next;
        pBehind=pBehind->next;
    }

    return pBehind;
}

int main(){
    int num;
    cin >> num;
    listNode* head = new listNode(num);
    listNode* p=head;

    while (cin >> num)
    {
        listNode* q =new listNode(num);
        p->next=q;
        p=p->next;
    }

    listNode* m =head;
    while (m!= nullptr)
    {
        cout << m->val << endl;
        m=m->next;
    }

    cin >> num;
    cout << "last " << num << endl;
    // flush cin
    cin.clear();
    cin.ignore(INT_MAX,'\n');
    cin >> num;
    cout << endl;
    cout << "last " << num << endl;
    cout << "last k node?" << endl;

    if(cin >> num)
    {
        p=findKthToTail(head,num);
        cout << p->val << endl;
    }
    return(0);
}