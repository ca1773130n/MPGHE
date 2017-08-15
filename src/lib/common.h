#pragma once

#define TEST_OK 0
#define TEST_FAIL -1
#define NUM_CHANNELS 3

// default options
#define DEFAULT_CAM_W		640
#define DEFAULT_CAM_H		360
#define DEFAULT_ALPHA		50
#define DEFAULT_ALPHA_STEP	5
#define DEFAULT_BETA		0.01
#define DEFAULT_BETA_STEP	0.005
#define DEFAULT_M			1
#define DEFAULT_METHOD	LM_LAPLACIAN
#define DEFAULT_SOURCE		IS_WEBCAM

// internal options
#define HISTOGRAM_RANGE	256

enum LMType {
	LM_LAPLACIAN,
	LM_LOCAL_MEAN,
	NUM_LM_TYPES,
	LM_INVALID
};

union LocalMethod {
	enum LMType type;
	int index;
};

enum InputSource {
	IS_IMAGE,
	IS_WEBCAM,
	IS_INVALID
};

struct TestOptions {
	double alpha;
	double beta;
	int m;
	char *filepath;
	size_t screenW;
	size_t screenH;
	size_t histWinW;
	size_t histWinH;
	LocalMethod method;
	InputSource source;
};
