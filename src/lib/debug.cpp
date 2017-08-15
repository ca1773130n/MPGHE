#pragma once
#include <debug.h>
#include <string.h>

namespace utils {
	bool Logger::init(int cLevel, int fLevel, const char *logFilePath) {
		bool result = true;
		mLevels[LOG_DEVICE_CONSOLE] = cLevel;
		mLevels[LOG_DEVICE_FILE] = fLevel;

		mFP = fopen(logFilePath, "w");
		if (!mFP) {
			printf("unable to open a file to log.\n");
			result = false;
		}

		mTags[0] = "[VERBOSE]";
		mTags[1] = "[DEBUG]   ";
		mTags[2] = "[INFO]    ";
		mTags[3] = "[WARNING] ";
		mTags[4] = "[ERROR]   ";

		return result;
	}

	void _LOG(bool showTag, int level, const char *fmt, ...) {
		char buffer[LOG_BUFSIZE] = { 0, };
		bool bufferCopied = false;
		Logger *logger = Logger::getInstance();
		int logLevel = logger->getLevel(LOG_DEVICE_CONSOLE);
		int fileLogLevel = logger->getLevel(LOG_DEVICE_FILE);
		FILE *fp = logger->getFP();

		va_list lpStart;
		va_start(lpStart, fmt);

		if (showTag)
			strcpy(buffer, logger->getTag(level));

		if (logLevel <= level) {
			vsprintf(buffer + strlen(buffer), fmt, lpStart);
			printf(buffer);
			bufferCopied = true;
		}

		if (fileLogLevel <= level) {
			if (!bufferCopied)
				vsprintf(buffer + strlen(buffer), fmt, lpStart);
			fprintf(fp, buffer);
			logger->sync();
		}
	}
}
