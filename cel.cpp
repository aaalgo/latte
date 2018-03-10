#include <iostream>
#include <file/CELFileData.h>

// trying to read CEL file using affy sdk

using namespace std;

int main (int argc, char *argv[]) {
    affxcel::CCELFileData cel;
    cel.ReadEx(argv[1]);
    cout << cel.GetAlgorithmParameters() << endl;
    return 0;
}

