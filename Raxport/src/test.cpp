
#include <iostream>
#include <windows.h>
#include <Tchar.h>
using namespace std;

namespace MSFileReaderLib
{
    class IXRawfile5
    {
    public:
        virtual ~IXRawfile5();
        virtual long Open(TCHAR *path) = 0;
    };
}

class XRawfile
{
public:
    virtual long Open(TCHAR *path);
};

// https://blog.csdn.net/blacet/article/details/50696109
void testDLL()
{
    class __declspec(uuid("{1d23188d-53fe-4c25-b032-dc70acdbdc02}")) Raw;
    static const CLSID CLSID_Raw = __uuidof(Raw);
    typedef HRESULT(__stdcall * pfn)(REFCLSID, REFIID, void **);
    pfn fn = NULL;
    HINSTANCE hdllInst = LoadLibrary("XRawfile2_x64.dll");
    fn = (pfn)GetProcAddress(hdllInst, "DllGetClassObject");
    if (fn != 0)
    {
        IClassFactory *pcf = NULL;
        HRESULT hr = (fn)(CLSID_Raw, IID_IClassFactory, (void **)&pcf);
        if (SUCCEEDED(hr) && (pcf != NULL))
        {
            MSFileReaderLib::IXRawfile5Ptr comRawFile = NULL;
            hr = pcf->CreateInstance(NULL, IID_IFoo, (void **)&comRawFile);
            if (SUCCEEDED(hr) && (pFoo != NULL))
            {
                comRawFile.CreateInstance("MSFileReader.XRawfile.1");
                printf("777");
                comRawFile.Release();
            }
            pcf->Release();
        }
    }
    FreeLibrary(hdllInst);
}

int main()
{
    testDLL();
    HINSTANCE hDll = ::LoadLibrary("XRawfile2_x64.dll");
    if (hDll != NULL)
        printf("dll load succeed \n");
    else
    {
        printf("dll load failed \n");
        return 0;
    }
    typedef MSFileReaderLib::IXRawfile5 *(*IXRawfile5PtrInstance)();
    IXRawfile5PtrInstance instance = (IXRawfile5PtrInstance)GetProcAddress(hDll, "IXRawfile");
    cout << "ret="
         << "666" << endl;
    MSFileReaderLib::IXRawfile5 *rawPtr = instance();
    cout << "ret="
         << "ret" << endl;
    long ret = rawPtr->Open(_T("D:\\work\\202202\\Ecoli\\201019_CP_6801_E50.raw"));
    cout << "ret=" << ret << endl;
    FreeLibrary(hDll);
}