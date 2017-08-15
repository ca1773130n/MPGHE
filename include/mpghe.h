#pragma once

#include <opencv2/core.hpp>

namespace cv {
	int multiPeakGHE(const Mat& input, Mat& output, const TestOptions& opt);
}