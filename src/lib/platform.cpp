#include "platform.h"

#if defined(_WIN32)
bool FileExists(const char *filePath) {
	DWORD dwAttrib = GetFileAttributes(filePath);
	return (dwAttrib != INVALID_FILE_ATTRIBUTES && !(dwAttrib & FILE_ATTRIBUTE_DIRECTORY));
}

void GetScreenResolution(int& width, int& height) {
	RECT desktop;
	const HWND hDesktop = GetDesktopWindow();
	GetWindowRect(hDesktop, &desktop);
	width = desktop.right;
	height = desktop.bottom;
}
#elif defined(__APPLE__) || defined(linux) || defined(__linux__) || defined(__unix__)
#include <fcntl.h>
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xresource.h>

bool FileExists(const char *filePath) {
	return (access(filePath, F_OK ) != -1);
}

void GetScreenResolution(int& width, int& height) {
	Display* d = XOpenDisplay(NULL);
	Screen* s = DefaultScreenOfDisplay(d);
	width = s->width;
	height = s->height;
}
#endif
