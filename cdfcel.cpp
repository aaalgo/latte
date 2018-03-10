#include <iostream>
#include <file/CDFFileData.h>
#include <file/CELFileData.h>

// trying to read CEL file using affy sdk

using namespace std;

int main (int argc, char *argv[]) {
    if (argc != 2) return 0;
    affxcel::CCELFileData cel;
    cel.ReadEx(argv[1]);
    if (cel.GetCols() != 1164) throw 0;
    if (cel.GetRows() != 1164) throw 0;

    affxcdf::CCDFFileData cdf;
    cdf.SetFileName("HG-U133_Plus_2.cdf");
    cdf.Read();
    affxcdf::CCDFFileHeader &header = cdf.GetHeader();
    if (header.GetCols() != 1164) throw 0;
    if (header.GetRows() != 1164) throw 0;
    if (header.GetNumProbeSets() != 54675) throw 0;
    cout << header.GetNumProbeSets() << endl;
    cout << header.GetNumQCProbeSets() << endl;
    int np = header.GetNumProbeSets();
    // foreach probe set
    for (int i = 0; i < np; ++i) {
        affxcdf::CCDFProbeSetInformation info;
        cdf.GetProbeSetInformation(i, info);
        if (info.GetProbeSetType() != affxcdf::ExpressionProbeSetType) throw 0;
        if (info.GetNumGroups() != 1) throw 0;
        affxcdf::CCDFProbeGroupInformation groupInfo;
        info.GetGroupInformation(0, groupInfo);
        if (info.GetNumLists() != groupInfo.GetNumLists()) throw 0;
        if (info.GetNumCells() != groupInfo.GetNumCells()) throw 0;
        if (groupInfo.GetChannel() != 0) throw 0;
        if (groupInfo.GetDirection() != affxcdf::AntiSenseDirection) throw 0;
        if (groupInfo.GetRepType() != affxcdf::UnknownRepType) throw 0;
        if (groupInfo.GetWobbleSituation() != 0) throw 0;
        if (groupInfo.GetAlleleCode() != 0) throw 0;
        if (groupInfo.GetNumCellsPerList() != 2) throw 0;

        cout << groupInfo.GetName() << endl;
        for (int j = 0; j < groupInfo.GetNumCells(); ++j) {
            affxcdf::CCDFProbeInformation probeInfo;
            groupInfo.GetCell(j, probeInfo);
            if (probeInfo.GetProbeLength() != 0) throw 0;
            if (probeInfo.GetProbeGrouping() != 0) throw 0;
            if (probeInfo.GetListIndex() != probeInfo.GetExpos()) throw 0;
            bool is_pm = (probeInfo.GetPBase() != probeInfo.GetTBase());
            int x = probeInfo.GetX();
            int y = probeInfo.GetY();
            affxcel::CELFileEntryType e;
            cel.GetEntry(x, y, e);
            cout << '\t' << x << ' ' << y << ' ' << (is_pm ? "pm" : "mm") << ' ' << e.Intensity << endl;
        }
    }

    return 0;
}

