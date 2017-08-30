#include <opencv2/imgproc.hpp>
#include <opencv2/videoio.hpp>
#include <opencv2/highgui.hpp>
#include <cassert>
#include <iostream>
#include <inttypes.h>

#include "common.h"
#include "platform.h"
#include "mpghe.h"
#include "debug.h"
#include "profile.h"

#if defined(_WIN32)
#include "getopt_win.h"
#else
#include <getopt.h>
#endif

#define WINNAME_INPUT "input"
#define WINNAME_OUTPUT_HE "output (HE)"
#define WINNAME_OUTPUT_MPGHE "output (MPGHE)"
#define WINNAME_HISTOGRAM_INPUT "histogram (input)"
#define WINNAME_HISTOGRAM_HE "histogram (HE)"
#define WINNAME_HISTOGRAM_MPGHE "histogram (MPGHE)"

using namespace utils;

static int SetupOptions(int argc, char *argv[], TestOptions& options)
{
	assert(argc >= 1);
	assert(argv != NULL);

	int error = TEST_OK;

	// set to default values
	options.alpha = DEFAULT_ALPHA;
	options.beta = DEFAULT_BETA;
	options.m = DEFAULT_M;
	options.method.type = DEFAULT_METHOD;
	options.source = DEFAULT_SOURCE;

	if (argc == 1) {
		printf("Usage: testMPGHE {options} {filepath}\n");
		printf("\n");
		printf("Options\n");
		printf("-a : alpha value (double)\n");
		printf("-b : beta value (double)\n");
		printf("-f : input filebeta value (double)\n");
		printf("-i : set input source (string, [image|video|webcam]) \n");
		printf("-m : local window size (integer)\n");
		printf("-v : local method (string, [laplacian|local_mean])\n");
		printf("\n");
		printf("Shortcuts to modify options\n");
		printf("[A|a] : increase / decrease the alpha\n");
		printf("[B|b] : increase / decrease the beta\n");
		printf("[M|m] : increase / decrease the m\n");
		printf("[V|v] : change the local method by increasing / decreasing the method index\n");
		error = TEST_FAIL;
	} else {
		int opt = -1;
		char *endPtr;

		optarg = NULL;
		optind = 1;
		opterr = 0;	
		optopt = 0;

		while ((opt = getopt(argc, argv, ":a:b:i:m:v:")) != -1) {
			switch (opt) {
			case 'a':
				endPtr = NULL;
				options.alpha = strtoimax(optarg, &endPtr, 10);
				break;
			case 'b':
				endPtr = NULL;
				options.beta = atof(optarg);
				break;
			case 'i':
				if (!strncmp("image", optarg, 5)) {
					options.source = IS_IMAGE;
				} else if (!strncmp("webcam", optarg, 6)) {
					options.source = IS_WEBCAM;
				} else if (!strncmp("video", optarg, 5)) {
					options.source = IS_VIDEO;
				} else {
					printf("input source %s is not supported.\n");
					options.source = IS_INVALID;
					error = TEST_FAIL;
				}
				break;
			case 'm':
				endPtr = NULL;
				options.m = strtoimax(optarg, &endPtr, 10);
				break;
			case 'v':
				if (!strncmp("laplacian", optarg, 9)) {
					options.method.type = LM_LAPLACIAN;
				} else if (!strncmp("local_mean", optarg, 10)) {
					options.method.type = LM_LOCAL_MEAN;
				} else {
					printf("local method %s is not supported.\n");
					options.method.type = LM_INVALID;
					error = TEST_FAIL;
				}
				break;
			case '?':
				printf("ignored invalid option : -%c.\n", optopt);
				break;
			}
		}

		if (optind < argc) {
			char *filePath = argv[optind];
			if (FileExists(filePath))
				options.filepath = argv[optind];
			else {
				printf("input file %s does not exists.\n", filePath);
				error = TEST_FAIL;
			}
		}

		GetScreenResolution(options.screenW, options.screenH);
	}

	return error;
}

static void Preprocess(const cv::Mat& frame, cv::Mat& resized, cv::Mat *channels, TestOptions& opt) {
	double frameRatio = (double)frame.rows / frame.cols;
	cv::Size s;
	
	opt.histWinW = opt.screenW / 3;
	opt.histWinH = opt.screenH / 4;
	s.width = opt.screenW / 3;
	s.height = s.width * frameRatio;

	while (opt.screenH - s.height < opt.histWinH) {
		s.width *= 0.8;
		s.height = s.width * frameRatio;
	}

	cv::resize(frame, resized, s);
	cv::split(resized, channels);

	// set window size considering title bar and task bar
	opt.histWinH = opt.screenH - s.height - 32 - 100;
}

static void ProcessKey(const int key, TestOptions& options) {
	switch (key) {
	case 'A': options.alpha+= DEFAULT_ALPHA_STEP; break;
	case 'a': if ((options.alpha -= DEFAULT_ALPHA_STEP) < 0) options.alpha = 0; break;
	case 'B': options.beta += DEFAULT_BETA_STEP; break;
	case 'b': if ((options.beta -= DEFAULT_BETA_STEP) < 0) options.beta = 0; break;
	case 'M': options.m++; break;
	case 'm': if (--options.m < 1) options.m = 1; break;
	case 'V': ++options.method.index; options.method.index %= NUM_LM_TYPES;  break;
	case 'v': if (--options.method.index < 0) options.method.index = NUM_LM_TYPES - 1; break;
	}
}

static void PrintResult(const TestOptions& opt, const float timeMs) {
	std::cout << "\r";
	std::cout << "Settings ";
	std::cout << "[alpha: " << opt.alpha << "]";
	std::cout << "[beta: " << opt.beta << "]";
	std::cout << "[m: " << opt.m << "]";
	std::cout << "[method: " << opt.method.index << "]";
	std::cout << " / Time [" << timeMs << " ms]";
	std::cout << std::flush;
}

int main(int argc, char *argv[])
{
	TestOptions options;
	Logger *logger = Logger::getInstance();
	logger->init(LOG_LEVEL_INFO, LOG_LEVEL_DEBUG, "log.txt");

	if (SetupOptions(argc, argv, options) == TEST_OK) {
		int error = TEST_OK;
		cv::VideoCapture *vidCap = NULL;
		cv::Mat frame, input;
		cv::Mat outputHE, outputMPGHE;
		cv::Mat histoInput[NUM_CHANNELS];
		cv::Mat histoHE[NUM_CHANNELS];
		cv::Mat histoMPGHE[NUM_CHANNELS];
		cv::Mat histoDrawInput, histoDrawHE, histoDrawMPGHE;
		cv::Mat srcChannels[NUM_CHANNELS];
		cv::Mat outChannelsHE[NUM_CHANNELS];
		cv::Mat outChannelsMPGHE[NUM_CHANNELS];

		cv::namedWindow(WINNAME_INPUT);
		cv::namedWindow(WINNAME_OUTPUT_HE);
		cv::namedWindow(WINNAME_OUTPUT_MPGHE);
		cv::namedWindow(WINNAME_HISTOGRAM_INPUT);
		cv::namedWindow(WINNAME_HISTOGRAM_HE);
		cv::namedWindow(WINNAME_HISTOGRAM_MPGHE);

        const char *videoInput = nullptr;
		while (error == TEST_OK) {
			Timer timer("processing time of a frame");
			
			switch (options.source) {
			case IS_IMAGE:
				if (!frame.cols) {
					frame = cv::imread(options.filepath);
					Preprocess(frame, input, srcChannels, options);
				}
				break;
            case IS_VIDEO:
                videoInput = options.filepath;
			case IS_WEBCAM:
				if (!vidCap)
                    if (videoInput)
    					vidCap = new cv::VideoCapture(videoInput);
    				else
                        vidCap = new cv::VideoCapture(0);

				if (!vidCap->isOpened()) {
					printf("unable to open the webcam.\n");
					error = TEST_FAIL;
				}

				if (!vidCap->read(frame)) {
					printf("webcam disconnected.\n");
					error = TEST_FAIL;
				} else Preprocess(frame, input, srcChannels, options);
				break;
			default:
				// never be here
				assert(false);
			}

			if (!histoDrawInput.cols) {
				cv::Size hWinSize(options.histWinW, options.histWinH);
				histoDrawInput.create(hWinSize, CV_8UC3);
				histoDrawHE.create(hWinSize, CV_8UC3);
				histoDrawMPGHE.create(hWinSize, CV_8UC3);
			}

			// calculate MPGHE
			for (size_t c = 0; c < NUM_CHANNELS; ++c) {
				error = cv::multiPeakGHE(srcChannels[c], outChannelsMPGHE[c], options);
				if (error != TEST_OK) {
					MLOGE("multiPeakGHE() for channel %d failed.\n", c);
				}
			}
			cv::merge(outChannelsMPGHE, NUM_CHANNELS, outputMPGHE);

			double timeMs = timer.getElapsed();
			ProcessKey(cv::waitKey(1), options);
			PrintResult(options, timeMs);

			// calculate the reference Histogram Equalization result
			for (size_t c = 0; c < NUM_CHANNELS; ++c)
				cv::equalizeHist(srcChannels[c], outChannelsHE[c]);
			cv::merge(outChannelsHE, NUM_CHANNELS, outputHE);

			// calculate and draw histograms both HE and MPGHE
			int channelIndex[1];
			int histSize[1] = { HISTOGRAM_RANGE };
			float histRange[2] = { 0.0, HISTOGRAM_RANGE - 1 };
			const float* histRanges[1] = { histRange };
			size_t histWinH = options.histWinH * 0.8;

			for (size_t c = 0; c < NUM_CHANNELS; ++c) {
				channelIndex[0] = c;
				cv::calcHist(&input, 1, channelIndex, cv::Mat(), histoInput[c], 1, histSize, histRanges);
				cv::calcHist(&outputHE, 1, channelIndex, cv::Mat(), histoHE[c], 1, histSize, histRanges);
				cv::calcHist(&outputMPGHE, 1, channelIndex, cv::Mat(), histoMPGHE[c], 1, histSize, histRanges);
				cv::normalize(histoInput[c].clone(), histoInput[c], 0, histWinH, cv::NORM_MINMAX);
				cv::normalize(histoHE[c].clone(), histoHE[c], 0, histWinH, cv::NORM_MINMAX);
				cv::normalize(histoMPGHE[c].clone(), histoMPGHE[c], 0, histWinH, cv::NORM_MINMAX);
			}

			// draw histograms
			int valCurr, valPrev;
			cv::Scalar red(0, 0, 255);
			cv::Scalar green(0, 255, 0);
			cv::Scalar blue(255, 0, 0);
			cv::Scalar black(0, 0, 0);
			cv::Scalar *colors[NUM_CHANNELS] = { &blue, &green, &red };

			double pointRatio = (double)options.histWinW / HISTOGRAM_RANGE;
			histoDrawInput.setTo(black);
			histoDrawHE.setTo(black);
			histoDrawMPGHE.setTo(black);

			for (size_t i = 1; i < HISTOGRAM_RANGE; ++i) {
				for (size_t c = 0; c < NUM_CHANNELS; ++c) {
					valPrev = histoInput[c].at<float>(i - 1, 0);
					valCurr = histoInput[c].at<float>(i, 0);
					cv::line(histoDrawInput, cv::Point((i - 1) * pointRatio, options.histWinH - valPrev), cv::Point(i * pointRatio, options.histWinH - valCurr), *colors[c]);
					valPrev = histoHE[c].at<float>(i - 1, 0);
					valCurr = histoHE[c].at<float>(i, 0);
					cv::line(histoDrawHE, cv::Point((i - 1) * pointRatio, options.histWinH - valPrev), cv::Point(i * pointRatio, options.histWinH - valCurr), *colors[c]);
					valPrev = histoMPGHE[c].at<float>(i - 1, 0);
					valCurr = histoMPGHE[c].at<float>(i, 0);
					cv::line(histoDrawMPGHE, cv::Point((i - 1) * pointRatio, options.histWinH - valPrev), cv::Point(i * pointRatio, options.histWinH - valCurr), *colors[c]);
				}
			}

			// show the result images
			cv::Size fWinSize(options.histWinW, options.screenH - input.rows);
			cv::moveWindow(WINNAME_INPUT, 0, 0);
			cv::moveWindow(WINNAME_OUTPUT_HE, fWinSize.width, 0);
			cv::moveWindow(WINNAME_OUTPUT_MPGHE, fWinSize.width * 2, 0);
			cv::moveWindow(WINNAME_HISTOGRAM_INPUT, 0, input.rows + 32);
			cv::moveWindow(WINNAME_HISTOGRAM_HE, fWinSize.width, input.rows + 32);
			cv::moveWindow(WINNAME_HISTOGRAM_MPGHE, fWinSize.width * 2, input.rows + 32);

			cv::imshow(WINNAME_INPUT, input);
			cv::imshow(WINNAME_OUTPUT_HE, outputHE);
			cv::imshow(WINNAME_OUTPUT_MPGHE, outputMPGHE);
			cv::imshow(WINNAME_HISTOGRAM_INPUT, histoDrawInput);
			cv::imshow(WINNAME_HISTOGRAM_HE, histoDrawHE);
			cv::imshow(WINNAME_HISTOGRAM_MPGHE, histoDrawMPGHE);
		}
	}

	return 0;
}
