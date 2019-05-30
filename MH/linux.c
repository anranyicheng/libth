#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int
parseLine (char* line) {
  // This assumes that a digit will be found and the line ends in " Kb".
  int i = strlen(line);
  const char* p = line;
  while (*p <'0' || *p > '9') p++;
  line[i-3] = '\0';
  i = atoi(p);
  return i;
}

long
get_current_memory_usage () {
  FILE* file = fopen("/proc/self/status", "r");
  long result = -1;
  char line[128];

  while (fgets(line, 128, file) != NULL) {
    // if (strncmp(line, "VmSize:", 7) == 0) {
    if (strncmp(line, "VmRSS:", 6) == 0) {
      result = parseLine(line);
      break;
    }
  }
  fclose(file);
  if (result > 0)
    return result * 1024;
  else
    return result;
}
