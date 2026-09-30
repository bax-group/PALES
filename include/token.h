#ifndef PALES_TOKEN_H
#define PALES_TOKEN_H

int cntToken(const char *text);
int cntTokenC(const char *text, const char *separators);
char *getToken(int index, const char *text);
char *getTokenC(int index, const char *text, const char *separators);
int getTokenI(int index, const char *text);
float getTokenR(int index, const char *text);

#endif
