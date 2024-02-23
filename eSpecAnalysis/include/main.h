#pragma once

#ifndef __main_h__
#define __main_h__

#include <string>

// extern "C" is necessary so that the analyze_espec_images function name 
// does not get mangled in the shared library
#ifdef __cplusplus
extern "C" {
#endif

int analyze_espec_images(const char* filepath_eScreenA, const char* filepath_eScreenB, const char* filepath_ePointing, const char* filepath_spectrum); 

#ifdef __cplusplus
}
#endif

#endif
