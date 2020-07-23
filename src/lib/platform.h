#pragma once
#if defined(WIN32) || defined(_WIN32) || defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#include <Windows.h>
#elif defined(linux) || defined(__linux__) || defined(__unix__)
#include <sys/stat.h>
#include <X11/Xlib.h>
#endif

bool FileExists(const char *filePath);
void GetScreenResolution(int& width, int& height);
