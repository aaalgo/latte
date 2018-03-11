#include <map>
#include <iostream>
#include <file/CDFFileData.h>
#include <file/CELFileData.h>
#include <opencv2/opencv.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>




using std::map;
using namespace std;
namespace ba = boost::accumulators;

vector<cv::Vec3b> PALETTE_TABLEAU20{
	{31, 119, 180}, {174, 199, 232}, {255, 127, 14}, {255, 187, 120},
	{44, 160, 44}, {152, 223, 138}, {214, 39, 40}, {255, 152, 150},
	{148, 103, 189}, {197, 176, 213}, {140, 86, 75}, {196, 156, 148},
	{227, 119, 194}, {247, 182, 210}, {127, 127, 127}, {199, 199, 199},
	{188, 189, 34}, {219, 219, 141}, {23, 190, 207}, {158, 218, 229} 
};   

// trying to read CEL file using affy sdk
class CDF: map<string, cv::Rect> {
    void save_and_discard (string const &path, cv::Mat image, float max) {
        //cv::normalize(image, image, 0, 255, cv::NORM_MINMAX);
        image *= 255/max;
        image.convertTo(image, CV_8U);
        cv::applyColorMap(image, image, cv::COLORMAP_JET);
        cv::imwrite(path, image);
    }
    void stat_and_discard (vector<float> &vs, float *mean, float *sigma) {
        sort(vs.begin(), vs.end());
        vs.pop_back();
        ba::accumulator_set<float, ba::stats<ba::tag::variance>> acc;
        for (unsigned j = 1; j < vs.size(); ++j) {
            acc(vs[j]);
        }
        *mean = ba::mean(acc);
        *sigma = sqrt(ba::variance(acc));
    }
public:
    CDF () {
        int H = 1164;
        int W = 1164;
        affxcel::CCELFileData cel;
        cel.ReadEx("cel");
        if (cel.GetCols() != W) throw 0;
        if (cel.GetRows() != H) throw 0;

        affxcdf::CCDFFileData cdf;
        cdf.SetFileName("HG-U133_Plus_2.cdf");
        cdf.Read();
        affxcdf::CCDFFileHeader &header = cdf.GetHeader();
        if (header.GetCols() != W) throw 0;
        if (header.GetRows() != H) throw 0;
        if (header.GetNumProbeSets() != 54675) throw 0;
        cout << header.GetNumProbeSets() << endl;
        cout << header.GetNumQCProbeSets() << endl;
        int np = header.GetNumProbeSets();
        // foreach probe set
		//cv::Mat loc(H, W, CV_8UC3, cv::Scalar(0,0,0));
        int H0 = 1700;
        cv::Mat raw(H, W, CV_32F, cv::Scalar(0));
        cv::Mat reorder(H0, W, CV_32F, cv::Scalar(0));
        float max = 0;

        {
            cv::Mat image(H, W, CV_32F, cv::Scalar(0));
            float *ptr = image.ptr<float>(0);
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    affxcel::CELFileEntryType e;
                    cel.GetEntry(j, i, e);
                    float v = std::log(e.Intensity+1);
                    *ptr = v;
                    max = std::max(max, v);
                    ++ptr;
                }
            }
            save_and_discard("/home/wdong/public_html/ma/raw.png", image, max);
        }

        int y_off = 0;
        int x_off = 0;
        int i = 0;

        cv::Mat outliers(H, W, CV_32F, cv::Scalar(0));
        cv::Mat inliers(H, W, CV_32F, cv::Scalar(0));
        for (; i < np; ++i) {
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

            string name = groupInfo.GetName();
            //cout << groupInfo.GetName() << endl;
            int minx = W, maxx = 0, miny = H, maxy = 0;
            if (x_off + groupInfo.GetNumCells()/2 + 1 > W) {
                y_off += 3;
                x_off = 0;
            }
            if (y_off >= H0) break;
            int x0 = x_off;
            int x1 = x_off;
            vector<float> vs_pm;
            vector<float> vs_mm;
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
                float v = std::log(e.Intensity+1);
                if (is_pm) {
                    reorder.at<float>(y_off, x0++) = v;
                    vs_pm.push_back(v);
                }
                else {
                    reorder.at<float>(y_off+1, x1++) = v;
                    vs_mm.push_back(v);
                }
                //loc.at<cv::Vec3b>(y, x) = PALETTE_TABLEAU20[i];
                //cout << '\t' << x << ' ' << y << ' ' << (is_pm ? "pm" : "mm") << ' ' << e.Intensity << endl;
            }
            float mean_pm, sigma_pm;
            float mean_mm, sigma_mm;
            stat_and_discard(vs_pm, &mean_pm, &sigma_pm);
            stat_and_discard(vs_mm, &mean_mm, &sigma_mm);

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
                float v = std::log(e.Intensity+1);
                float mean = is_pm ? mean_pm : mean_mm;
                float sigma = is_pm ? sigma_pm : sigma_mm;
                if (abs(v - mean) >= 3 * sigma) {
                    outliers.at<float>(y, x) = v;
                }
                else {
                    inliers.at<float>(y, x) = v;
                }
            }
        }
        cerr << i << " of " << np << " drawn at y= " << y_off << "." << endl;
        save_and_discard("/home/wdong/public_html/ma/reorder.png", reorder, max);
        save_and_discard("/home/wdong/public_html/ma/outliers.png", outliers, max);
        save_and_discard("/home/wdong/public_html/ma/inliers.png", inliers, max);
    }
};


int main (int argc, char *argv[]) {
    CDF cdf;


    return 0;
}

