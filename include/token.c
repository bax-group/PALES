#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "token.h"

#define TOKEN_BUFFER_SIZE 4096

static char tokenBuffer[TOKEN_BUFFER_SIZE];

static int isSeparator(int c, const char *separators)
{
    if (separators)
        return strchr(separators, c) != NULL;
    return isspace((unsigned char)c) || c == ',';
}

int cntTokenC(const char *text, const char *separators)
{
    int count;

    count = 0;
    if (!text)
        return 0;

    while (*text) {
        while (*text && isSeparator((unsigned char)*text, separators))
            ++text;
        if (!*text)
            break;
        ++count;
        while (*text && !isSeparator((unsigned char)*text, separators))
            ++text;
    }
    return count;
}

int cntToken(const char *text)
{
    return cntTokenC(text, NULL);
}

char *getTokenC(int index, const char *text, const char *separators)
{
    int current;
    size_t length;
    const char *start;

    tokenBuffer[0] = '\0';
    if (!text || index < 1)
        return tokenBuffer;

    current = 0;
    while (*text) {
        while (*text && isSeparator((unsigned char)*text, separators))
            ++text;
        if (!*text)
            break;
        start = text;
        while (*text && !isSeparator((unsigned char)*text, separators))
            ++text;
        if (++current == index) {
            length = (size_t)(text - start);
            if (length >= sizeof(tokenBuffer))
                length = sizeof(tokenBuffer) - 1;
            memcpy(tokenBuffer, start, length);
            tokenBuffer[length] = '\0';
            return tokenBuffer;
        }
    }
    return tokenBuffer;
}

char *getToken(int index, const char *text)
{
    return getTokenC(index, text, NULL);
}

int getTokenI(int index, const char *text)
{
    return (int)strtol(getToken(index, text), NULL, 10);
}

float getTokenR(int index, const char *text)
{
    return (float)strtod(getToken(index, text), NULL);
}
