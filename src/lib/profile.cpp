#if defined(_WIN32) || defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#include <Windows.h>
#include <time.h>
#endif

#include "profile.h"

namespace utils {
#if defined(_WIN32) || defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
	int gettimeofday(struct timeval *tv, struct timezone *tz)
	{
		FILETIME ft;
		unsigned __int64 tmpres = 0;
		static int tzflag;

		if (NULL != tv)
		{
			GetSystemTimeAsFileTime(&ft);

			tmpres |= ft.dwHighDateTime;
			tmpres <<= 32;
			tmpres |= ft.dwLowDateTime;

			tmpres -= DELTA_EPOCH_IN_MICROSECS;
			tmpres /= 10;

			tv->tv_sec = (tmpres / 1000000UL);
			tv->tv_usec = (tmpres % 1000000UL);
		}

		if (NULL != tz)
		{
			if (!tzflag)
			{
				_tzset();
				tzflag++;
			}
			tz->tz_minuteswest = _timezone / 60;
			tz->tz_dsttime = _daylight;
		}

		return 0;
	}
#endif

	void StopWatch::startTimer(const char *name) {
		std::string strName(name);
		if (mTimers[strName] != NULL) {
			mTimers[strName]->reset();
		}
		else {
			Timer *newTimer = new Timer(strName);
			mTimers[strName] = newTimer;
		}
	}

	void StopWatch::pauseTimer(const char *name) {
		std::string strName(name);
		mResults[strName] += mTimers[strName]->getElapsed();
	}

	void StopWatch::stopTimer(const char *name) {
		std::string strName(name);
		mResults[strName] = mTimers[strName]->getElapsed();
	}

	void StopWatch::printResults(bool flush) {
		for (auto t : mResults)
			MLOGD("%s : [%.3f] ms\n", t.first.c_str(), t.second);
	}
}
