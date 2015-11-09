#ifdef WIN32
    #include <direct.h>
    #define getCurrentWorkingDir _getcwd
#else
    #include <unistd.h>
    #define getCurrentWorkingDir getcwd
#endif

#define TCPATH "TCPATH"
#define TC3_HOME "TC3_HOME"

#ifdef WIN32
#define TEMP "TEMP"
#else
#define TEMP "TMPDIR"
#endif

#ifdef WIN32
#define SLASH "\\"
#else
#define SLASH "/"
#endif


/************************************************
*   Get path from the values of
*   Thermo-Calc environment variables
*   TC3_HOME (Thermo-Calc 3.0)
*   TCPATH Older versions and fallback
************************************************/
void getThermoCalcEnvironmentPath(char* pathBuffer);
//------------------------------------------------------------------------------

/************************************************
*   Get path to temp directory
*   If it can't find it - default
*   to current working directory
************************************************/
void getTempEnvironmentPath(char* pathBuffer);
//------------------------------------------------------------------------------
