#pragma once
#include <stdarg.h>
#include "pattern.h"

#define LOG_FILENAME			"log.txt"
#define LOG_LEVEL_VERBOSE		0
#define LOG_LEVEL_DEBUG			1
#define LOG_LEVEL_INFO			2
#define LOG_LEVEL_WARNING		3
#define LOG_LEVEL_ERROR			4
#define NUM_LOG_LEVELS			5
#define LOG_BUFSIZE				4096

namespace utils {
	enum LogDevice {
		LOG_DEVICE_CONSOLE = 0,
		LOG_DEVICE_FILE,
		NUM_LOG_DEVICES,
		LOG_DEVICE_INVALID
	};

	void _LOG(bool showTag, int level, const char *fmt, ...);

#define MLOGV(fmt, ...) _LOG(true, LOG_LEVEL_VERBOSE, fmt, __VA_ARGS__)
#define MLOGD(fmt, ...) _LOG(true, LOG_LEVEL_DEBUG, fmt, __VA_ARGS__)
#define MLOGI(fmt, ...) _LOG(true, LOG_LEVEL_INFO, fmt, __VA_ARGS__)
#define MLOGW(fmt, ...) _LOG(true, LOG_LEVEL_WARNING, fmt, __VA_ARGS__)
#define MLOGE(fmt, ...) _LOG(true, LOG_LEVEL_ERROR, fmt, __VA_ARGS__)

	class Logger;
	class Logger : public ISingleton<Logger> {
		public:
			Logger() {
			}
			~Logger() {
			}

			bool init(int cLevel, int fLevel, const char *logFilePath);

			inline int getLevel(int deviceID) {
				return mLevels[deviceID];
			}

			inline FILE *getFP(void) {
				return mFP;
			}

			inline const char *getTag(int level) {
				return mTags[level];
			}

			inline void sync(void) {
				std::fflush(mFP);
			}

		protected:
			FILE *mFP;
			int mLevels[NUM_LOG_DEVICES];
			char *mTags[NUM_LOG_LEVELS];
	};

}
