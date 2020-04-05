#pragma once

#include <map>
#include "debug.h"
#include "pattern.h"

#if defined(WIN32) || defined(_WIN32) || defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#elif defined(__APPLE__) || defined(linux) || defined(__linux__) || defined(__unix__)
#include <sys/time.h>
#endif

#define PROFILE(name, code) { AutoTimer _UniqueTimer(name); code }
#define PROFILE_START(name) StopWatch::getInstance()->startTimer(name);
#define PROFILE_END(name) StopWatch::getInstance()->stopTimer(name);
#define PROFILE_PAUSE(name) StopWatch::getInstance()->pauseTimer(name);
#define PROFILE_STAT() StopWatch::getInstance()->printResults();

namespace utils {

#if defined(WIN32) || defined(_WIN32) || defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
	struct timezone {
		int tz_minuteswest;
		int tz_dsttime;
	};

	struct timeval {
		long tv_sec;
		long tv_usec;
	};

	int gettimeofday(struct timeval *tv, struct timezone *tz);
#endif

	class Timer {
		public:
			Timer(std::string name) {
				mName.assign(name);
				reset();
			}

			void reset(void) {
				gettimeofday(&mStartTime, NULL);
			}

			float getElapsed(bool print = false) {
				long mSec = 0;
				gettimeofday(&mCurTime, NULL);
				mSec = (mCurTime.tv_sec - mStartTime.tv_sec) * 1000000 + (mCurTime.tv_usec - mStartTime.tv_usec);
				if (print) {
					MLOGD("Time elapsed for [%s] : %.8f msec\n", mName.c_str(), mSec / 1000.f);
				}
				return mSec / 1000.f;
			}

		protected:
			struct timeval mStartTime;
			struct timeval mCurTime;
			size_t mElapsedMsec;
			std::string mName;
	};

	class StopWatch;
	class StopWatch : public ISingleton<StopWatch> {
		public:
			StopWatch() {}
			~StopWatch() {
				for (auto t : mTimers) {
					delete t.second;
				}
			}

			void startTimer(const char *name);
			void pauseTimer(const char *name);
			void stopTimer(const char *name);
			void printResults(bool flush = true);

		protected:
			std::map<std::string, Timer *> mTimers;
			std::map<std::string, float> mResults;
	};

	class AutoTimer : public Timer {
		public:
			AutoTimer(char *name) : Timer(name) {}
			~AutoTimer() {
				getElapsed();
			}
	};
}
