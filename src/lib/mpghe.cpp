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
	static int calculateG(const Mat& input, const double alpha, const double beta, const LocalMethod& method, const int m, Mat& G, double& Gmin, double& Gmax);
	static int normalizeV(const Mat&V, Mat&v, const double beta, const LocalMethod& method, const int m);
	static double calculateV(const Mat& input, const int x, const int y, const int m);
	static void findLocalMaximums(const Mat& hist, const std::vector<size_t>& indices, std::vector<size_t>& maximums);
	static void findLocalMinimums(const Mat& hist, std::vector<size_t>& minimums, const double Gmin, const double Gmax, const double Hmin, const double Hmax);
	static void equalizeHistogram(const Mat& I, const Mat& hist, Mat& mapTable, std::vector<size_t>& minimums);

	int multiPeakGHE(const Mat& input, Mat& output, const TestOptions& opt) {
		int error = TEST_OK;
		const double alpha = opt.alpha;
		const double beta = opt.beta;
		const int m = opt.m;
		const LocalMethod method = opt.method;

		Mat u, v, w, p, V, G, Gf;
		Mat histoG, mapTable;
		double Gmin, Gmax, Hmin, Hmax;
		std::vector<size_t> localMins;

		int channelIndex[1] = { 0 };
		int histSize[1] = { 256 };
		float histRange[2] = { 0.0, 255.0 };
		const float* histRanges[1] = { histRange };

		output.create(input.rows, input.cols, CV_8U);

		// calculate G and its histogram
		error = calculateG(input, alpha, beta, method, m, G, Gmin, Gmax);
		if (error != TEST_OK) {
			return error;
		}
		G.convertTo(Gf, CV_32F);
		calcHist(&Gf, 1, channelIndex, Mat(), histoG, 1, histSize, histRanges);

		// find mid nodes (local minimums)
		minMaxLoc(histoG, &Hmin, &Hmax);
		findLocalMinimums(histoG, localMins, Gmin, Gmax, Hmin, Hmax);

		// equalize the histogram piecewise and independently
		equalizeHistogram(G, histoG, mapTable, localMins);

		// output the enhanced image
		for (int y = 0; y < input.rows; ++y)
		for (int x = 0; x < input.cols; ++x) {
			int index = G.at<double>(y, x);
			output.at<uchar>(y, x) = mapTable.at<float>(index, 0);
		}

		return error;
	}

	static int normalizeV(const Mat&V, Mat&v, const double beta, const LocalMethod& method, const int m) {
		int error = TEST_OK;
		double Vmin, Vmax, Vxy, vxy;
		minMaxLoc(V, &Vmin, &Vmax);

		for (int y = 0; y < V.rows; ++y)
		for (int x = 0; x < V.cols; ++x) {
			Vxy = V.at<double>(y, x);

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

			v.at<double>(y, x) = vxy;
		}

		return error;
	}

	static double calculateV(const Mat& input, const int x, const int y, const LocalMethod& method, const int m) {
		assert(m >= 0);

		double result = 0;
		double sum = 0;
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
		float prev, cur;
		float min = 987654321;
		float max = 0;
		bool goingUP = true;

		size_t numIndices = indices.size();
		for (size_t x = 1; x < numIndices; ++x) {
			prev = hist.at<float>(indices[x - 1], 0);
			cur = hist.at<float>(indices[x], 0);

			if (prev > max) {
				max = prev;
			}
			else if (prev < min) {
				min = prev;
			}

			if (cur > prev) {
				if (!goingUP) goingUP = true;
			}
			else if (cur < prev) {
				if (goingUP) {
					// found a local peak
					maximums.push_back(indices[x - 1]);
					goingUP = false;
				}
			}
		}
	}

	// [27] H.D. Cheng, Y. Sun, A hierarchical approach to color image segmentation using homogeneity, IEEE Trans. Image Process. 9 (12) (2000).
	static void findLocalMinimums(const Mat& hist, std::vector<size_t>& minimums, const double Gmin, const double Gmax, const double Hmin, const double Hmax) {
		assert(hist.type() == CV_32F);
		assert(hist.cols == 1);

		std::vector<size_t> histIndices;
		std::vector<size_t> localPeakIndices;
		std::vector<size_t> localSPeakIndices;
		std::vector<size_t> localValidPeakIndices;
		std::set<size_t> validPeakSet;

		// find all significant local peaks
		for (size_t i = 0; i < hist.rows; ++i)
			histIndices.push_back(i);
		findLocalMaximums(hist, histIndices, localPeakIndices);
		findLocalMaximums(hist, localPeakIndices, localSPeakIndices);

		// remove small peaks
		for (size_t i = 0; i < localSPeakIndices.size(); ++i) {
			size_t pi = localSPeakIndices[i];
			if (hist.at<float>(pi, 0) >= 0.05 *  Hmax) {
				validPeakSet.insert(pi);
				localValidPeakIndices.push_back(pi);
			}
		}
		
		// remove peaks too close each other
		size_t lastIndex;
		for (size_t i = 0; i < localValidPeakIndices.size() - 1; ++i) {
			size_t pi1 = localValidPeakIndices[i];
			size_t pi2 = localValidPeakIndices[i + 1];
			double h1 = hist.at<float>(pi1, 0);
			double h2 = hist.at<float>(pi2, 0);
			lastIndex = pi1;

			if (validPeakSet.find(pi2) != validPeakSet.end()) {
				if (pi2 - lastIndex < 15) {
					if (h1 > h2) {
						validPeakSet.erase(pi2);
					} else {
						validPeakSet.erase(lastIndex);
						lastIndex = pi2;
					}
				}
				else lastIndex = pi2;
			}
		}

		std::vector<size_t> validPeaks1(validPeakSet.size());
		std::copy(validPeakSet.begin(), validPeakSet.end(), validPeaks1.begin());

		// remove a peak which forms only a small valley with its neighbor peak
		for (size_t i = 0; i < validPeaks1.size() - 1; ++i) {
			size_t pi1 = validPeaks1[i];
			size_t pi2 = validPeaks1[i + 1];
			double h1 = hist.at<float>(pi1, 0);
			double h2 = hist.at<float>(pi2, 0);
			
			double sum = 0;
			for (size_t j = pi1 + 1; j < pi2; ++j)
				sum += hist.at<float>(j, 0);
			double havg = sum / (pi2 - pi1 - 1);

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
			double h1 = hist.at<float>(pi1, 0);
			double h2 = hist.at<float>(pi2, 0);
			double min = std::min(h1, h2);
			int minIndex = -1;

			for (size_t j = pi1; j < pi2; ++j) {
				double hj = hist.at<float>(j, 0);
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

	static void equalizeHistogram(const Mat& I, const Mat& hist, Mat& mapTable, std::vector<size_t>& minimums) {
		Mat cdf;
		float accum = hist.at<float>(0, 0);
		size_t sizeMins = minimums.size();
		cdf.create(hist.rows, 1, CV_32F);
		mapTable.create(hist.rows, 1, CV_32F);

		for (size_t i = 0; i < sizeMins - 1; ++i) {
			size_t start = minimums[i];
			size_t end = minimums[i + 1];
			cdf.setTo(0);
			accum = 0;

			for (size_t j = 0; j <= end - start; ++j) {
				accum += hist.at<float>(start + j, 0);
				cdf.at<float>(j, 0) = accum;
			}

			normalize(cdf.clone(), cdf, start, end, NORM_MINMAX);
			for (size_t j = 0; j <= end - start; ++j)
				mapTable.at<float>(start + j, 0) = cdf.at<float>(j, 0);
		}
	}

	static int calculateG(const Mat& input, const double alpha, const double beta, const LocalMethod& method, const int m, Mat& G, double& Gmin, double& Gmax) {
		int error = TEST_OK;
		double Imin, Imax,  pmin, pmax;
		Size s(input.cols, input.rows);
		Mat u, v, p, V;
		
		u.create(s, CV_64F);
		v.create(s, CV_64F);
		p.create(s, CV_64F);
		V.create(s, CV_64F);
		G.create(s, CV_64F);

		minMaxLoc(input, &Imin, &Imax);
		double wxy = alpha / (Imax - Imin);
		Gmin = 0;
		Gmax = 255;

		// calculate u(x, y) and V(x, y)
		for (int y = 0; y < input.rows; ++y)
		for (int x = 0; x < input.cols; ++x) {
			u.at<double>(y, x) = (input.at<uchar>(y, x) - Imin) / (Imax - Imin);
			V.at<double>(y, x) = calculateV(input, x, y, method, m);
		}

		// calculate v(x, y) by normalizing the V(x, y)
		error = normalizeV(V, v, beta, method, m);
		if (error != TEST_OK) {
			return error;
		}

		// calculate p(x, y)
		for (int y = 0; y < input.rows; ++y)
		for (int x = 0; x < input.cols; ++x) {
			double vxy = v.at<double>(y, x);
			double uxy = u.at<double>(y, x);
			p.at<double>(y, x) = uxy + wxy * vxy;
		}

		// calculate G(x, y)
		minMaxLoc(p, &pmin, &pmax);
		for (int y = 0; y < input.rows; ++y)
		for (int x = 0; x < input.cols; ++x) {
			double pxy = p.at<double>(y, x);
			double Gxy = Gmin + (pxy - pmin) * ((Gmax - Gmin) / (pmax - pmin));
			G.at<double>(y, x) = Gxy;
		}

		return error;
	}
}
