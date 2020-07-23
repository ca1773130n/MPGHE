#include <opencv2/imgproc.hpp>
#include <opencv2/videoio.hpp>
#include <opencv2/highgui.hpp>
#include <map>
#include <bitset>
#include <set>

#include "common.h"
#include "mpghe.h"
#include "debug.h"

using namespace utils;

namespace cv {
	static int calculateG(const Mat& input, const float alpha, const float beta, const LocalMethod& method, const int m, Mat& G, double& Gmin, double& Gmax);
	static int normalizeV(const Mat&V, Mat&v, const float beta, const LocalMethod& method, const int m);
	static float calculateV(const Mat& input, const int x, const int y, const LocalMethod& method, const int m);
	static void findLocalMaximums(const Mat& hist, const std::vector<size_t>& indices, std::vector<size_t>& maximums);
	static void findLocalMinimums(const Mat& hist, std::vector<size_t>& minimums, const float Gmin, const float Gmax, const float Hmin, const float Hmax);
	static void equalizeHistogram(const Mat& I, const Mat& hist, Mat& mapTable, std::vector<size_t>& minimums);

	int multiPeakGHE(const Mat& input, Mat& output, const TestOptions& opt) {
		int error = TEST_OK;
		const float alpha = opt.alpha;
		const float beta = opt.beta;
		const int m = opt.m;
		const LocalMethod method = opt.method;

		Mat u, v, w, p, V, G;
		Mat histoG, mapTable;
		double Gmin, Gmax, Hmin, Hmax;
		std::vector<size_t> localMins;

		int channelIndex[1] = { 0 };
		int histSize[1] = { 256 };
		float histRange[2] = { 0.0, 255.0 };
		const float* histRanges[1] = { histRange };

		// calculate G and its histogram
		error = calculateG(input, alpha, beta, method, m, G, Gmin, Gmax);
		if (error != TEST_OK) {
			return error;
		}
		calcHist(&G, 1, channelIndex, Mat(), histoG, 1, histSize, histRanges);
        
		// find mid nodes (local minimums)
		minMaxLoc(histoG, &Hmin, &Hmax);
		findLocalMinimums(histoG, localMins, Gmin, Gmax, Hmin, Hmax);

		// equalize the histogram piecewise and independently
		equalizeHistogram(G, histoG, mapTable, localMins);

		// output the enhanced image
		for (int y = 0; y < input.rows; ++y)
		for (int x = 0; x < input.cols; ++x) {
			uchar index = input.at<uchar>(y, x);
			output.at<uchar>(y, x) = mapTable.at<uchar>(index, 0);
		}

		return error;
	}

	static int normalizeV(const Mat&V, Mat&v, const float beta, const LocalMethod& method, const int m) {
		int error = TEST_OK;
		double Vmin, Vmax;
        float Vxy, vxy;
		minMaxLoc(V, &Vmin, &Vmax);

		for (int y = 0; y < V.rows; ++y)
		for (int x = 0; x < V.cols; ++x) {
			Vxy = V.at<float>(y, x);

			switch (method.type) {
			case LM_LAPLACIAN:
				if (Vxy < 0)
					vxy = -0.5 * pow(Vxy / Vmin, beta);
				else
					vxy = 0.5 * pow(Vxy / Vmax, beta);
				break;
			case LM_LOCAL_MEAN:
				vxy = pow((Vxy - Vmin) / (Vmax - Vmin), beta);
				break;
			default:
				MLOGE("unsupported local method %d.\n", method.index);
				error = TEST_FAIL;
			}

			v.at<float>(y, x) = vxy;
		}

		return error;
	}

	static float calculateV(const Mat& input, const int x, const int y, const LocalMethod& method, const int m) {
		assert(m >= 0);

		float result = 0;
		float sum = 0;
		int w = 1;
	
		if (method.type == LM_LOCAL_MEAN)
			w = m;

		for (int i = -w; i <= w; ++i)
		for (int j = -w; j <= w; ++j) {
			if (y + j < 0 || y + j >= input.rows || x + i < 0 || x + i >= input.cols) continue;
			sum += input.at<uchar>(y + j, x + i);
		}

		switch (method.type) {
		case LM_LAPLACIAN:
			result = input.at<uchar>(y, x) * 9 - sum;
			break;
		case LM_LOCAL_MEAN:
			result = sum / (m * m);
			break;
		}

		return result;
	}
	
	static void findLocalMaximums(const Mat& hist, const std::vector<size_t>& indices, std::vector<size_t>& maximums) {
		float prev, cur, next;

		size_t numIndices = indices.size();
		for (size_t x = 1; x < numIndices - 1; ++x) {
			prev = hist.at<float>(indices[x - 1], 0);
			cur = hist.at<float>(indices[x], 0);
            next = hist.at<float>(indices[x + 1], 0);

            if (cur > prev && cur > next)
                maximums.push_back(indices[x]);
		}
	}

	// [27] H.D. Cheng, Y. Sun, A hierarchical approach to color image segmentation using homogeneity, IEEE Trans. Image Process. 9 (12) (2000).
	static void findLocalMinimums(const Mat& hist, std::vector<size_t>& minimums, const float Gmin, const float Gmax, const float Hmin, const float Hmax) {
		assert(hist.type() == CV_32F);
		assert(hist.cols == 1);

		std::vector<size_t> histIndices;
		std::vector<size_t> localPeakIndices;
		std::vector<size_t> localSPeakIndices;
		std::vector<size_t> localValidPeakIndices;
		std::set<size_t> validPeakSet;

		// find all significant local peaks
		for (int i = 0; i < hist.rows; ++i)
			histIndices.push_back(i);
		findLocalMaximums(hist, histIndices, localPeakIndices);
		findLocalMaximums(hist, localPeakIndices, localSPeakIndices);

		// remove small peaks
		for (size_t i = 0; i < localSPeakIndices.size(); ++i) {
			size_t pi = localSPeakIndices[i];
			if (hist.at<float>(pi, 0) / Hmax >= 0.05) {
				validPeakSet.insert(pi);
				localValidPeakIndices.push_back(pi);
			}
		}
	
		// remove peaks too close each other
		for (size_t i = 0; i < localValidPeakIndices.size() - 1; ++i) {
			size_t pi1 = localValidPeakIndices[i];
			size_t pi2 = localValidPeakIndices[i + 1];
			float h1 = hist.at<float>(pi1, 0);
			float h2 = hist.at<float>(pi2, 0);

			if (validPeakSet.find(pi2) != validPeakSet.end()) {
				if (pi2 - pi1 <= 15) {
					if (h1 > h2)
						validPeakSet.erase(pi2);
					else
						validPeakSet.erase(pi1);
				}
			}
		}

		std::vector<size_t> validPeaks1(validPeakSet.size());
		std::copy(validPeakSet.begin(), validPeakSet.end(), validPeaks1.begin());

		// remove a peak which forms only a small valley with its neighbor peak
		for (size_t i = 0; i < validPeaks1.size() - 1; ++i) {
			size_t pi1 = validPeaks1[i];
			size_t pi2 = validPeaks1[i + 1];
			float h1 = hist.at<float>(pi1, 0);
			float h2 = hist.at<float>(pi2, 0);
			
			float sum = 0;
			for (size_t j = pi1; j <= pi2; ++j)
				sum += hist.at<float>(j, 0);
			float havg = sum / (pi2 - pi1 + 1);

			if ((havg / (0.5 * (h1 + h2))) > 0.75) {
				if (h1 > h2) validPeakSet.erase(pi2);
				else validPeakSet.erase(pi1);
			}
		}
		
		std::vector<size_t> validPeaks2(validPeakSet.size());
		std::copy(validPeakSet.begin(), validPeakSet.end(), validPeaks2.begin());
		
		// find local minimum between every two peaks
		minimums.push_back(Gmin);
		for (size_t i = 0; i < validPeaks2.size() - 1; ++i) {
			size_t pi1 = validPeaks2[i];
			size_t pi2 = validPeaks2[i + 1];
			float h1 = hist.at<float>(pi1, 0);
			float h2 = hist.at<float>(pi2, 0);
			float min = std::min(h1, h2);
			int minIndex = -1;

			for (size_t j = pi1; j < pi2; ++j) {
				float hj = hist.at<float>(j, 0);
				if (hj < min) {
					min = hj;
					minIndex = j;
				}
			}

			if (minIndex >= 0)
				minimums.push_back(minIndex);
		}
		minimums.push_back(Gmax);
	}

    static void equalizeHistogramSegment(const Mat& I, const Mat& hist, Mat& mapTable, size_t startIndex, size_t endIndex) {
        int histSize = endIndex - startIndex + 1;
        float size = 0;
 
        for (int i = 0; i < histSize; i++)
            size += hist.at<float>(startIndex + i, 0);

        int winSize = 15;
        int winMidSize = winSize / 2;
        
        float cumhistogram[histSize + winSize];
        memset(cumhistogram, 0, sizeof(float) * (histSize + winSize));
        cumhistogram[winMidSize] = hist.at<float>(startIndex, 0);

        for (int i = 1; i < histSize; i++)
            cumhistogram[winMidSize + i] = hist.at<float>(startIndex + i, 0) + cumhistogram[winMidSize + i - 1];

        float alpha = 1.0f / cumhistogram[histSize + winMidSize - 1];
        for (int i = winMidSize; i < histSize - winMidSize; ++i) {
            float mean = 0;
            for (int j = i - winMidSize; j <= (i + winMidSize); ++j)
                mean += cumhistogram[j];
            cumhistogram[i] = mean / winSize;
        }

        for (int i = 0; i < histSize; i++)
            mapTable.at<uchar>(startIndex + i, 0) = startIndex + saturate_cast<uchar>((endIndex - startIndex) * cumhistogram[i + winMidSize] * alpha);
    }

	static void equalizeHistogram(const Mat& I, const Mat& hist, Mat& mapTable, std::vector<size_t>& minimums) {
		size_t sizeMins = minimums.size();
		mapTable.create(hist.rows, 1, CV_8U);
        mapTable.setTo(0);

		for (size_t i = 0; i < sizeMins - 1; ++i) {
			size_t start = minimums[i];
			size_t end = minimums[i + 1];
            equalizeHistogramSegment(I, hist, mapTable, start, end);
		}
	}

	static int calculateG(const Mat& input, const float alpha, const float beta, const LocalMethod& method, const int m, Mat& G, double& Gmin, double& Gmax) {
		int error = TEST_OK;
		double Imin, Imax,  pmin, pmax;
		Size s(input.cols, input.rows);
		Mat u, v, p, V;
		
		u.create(s, CV_32F);
		v.create(s, CV_32F);
		p.create(s, CV_32F);
		V.create(s, CV_32F);
		G.create(s, CV_32F);

		minMaxLoc(input, &Imin, &Imax);
		float wxy = alpha / (Imax - Imin);
		Gmin = 0;
		Gmax = 255;

		// calculate u(x, y) and V(x, y)
		for (int y = 0; y < input.rows; ++y)
		for (int x = 0; x < input.cols; ++x) {
			u.at<float>(y, x) = (input.at<uchar>(y, x) - Imin) / (Imax - Imin);
			V.at<float>(y, x) = calculateV(input, x, y, method, m);
		}

		// calculate v(x, y) by normalizing the V(x, y)
		error = normalizeV(V, v, beta, method, m);
		if (error != TEST_OK) {
			return error;
		}

		// calculate p(x, y)
		for (int y = 0; y < input.rows; ++y)
		for (int x = 0; x < input.cols; ++x) {
			float vxy = v.at<float>(y, x);
			float uxy = u.at<float>(y, x);
			p.at<float>(y, x) = uxy + wxy * vxy;
		}

		// calculate G(x, y)
		minMaxLoc(p, &pmin, &pmax);
		for (int y = 0; y < input.rows; ++y)
		for (int x = 0; x < input.cols; ++x) {
			float pxy = p.at<float>(y, x);
			G.at<float>(y, x) = Gmin + (pxy - pmin) * ((Gmax - Gmin) / (pmax - pmin));
		}

		return error;
	}
}
