//---------------------------------------------------------------------------

#include <stdio.h>

#pragma hdrstop

#include "my_bcb.h"


//---------------------------------------------------------------------------

#pragma package(smart_init)

// **** MyString class mimicking limited AnsiString functionality **************

// initialise empty string
MyString::MyString()
{
for (int i = 0; i < S_SIZE; i++) sbuf[i] = 0;
length = 0;
}

// initialise string by copying provided char buf value
MyString::MyString(char* ini_string)
{
for (int i = 0; i < S_SIZE; i++) sbuf[i] = 0;
length = 0;
for (int i = 0; i < S_SIZE-1; i++) // do not exceed buffer size
    {
    sbuf[i] = ini_string[i];
    if (sbuf[i] == 0) break;
    length++;
    }
}

// initialise string by copying provided MyString
MyString::MyString(MyString* ini_string)
{
for (int i = 0; i < S_SIZE; i++) sbuf[i] = 0;
length = ini_string->length;
for (int i = 0; i < length; i++) sbuf[i] = ini_string->sbuf[i];
}

// initialise string by provided integer
MyString::MyString(int ini_n)
{
for (int i = 0; i < S_SIZE; i++) sbuf[i] = 0;
length = 0;
sprintf(sbuf, "%d", ini_n);
while (sbuf[length] != 0) length++;
}

// similar copy constructor
MyString::MyString(const MyString &ini_string)
{
for (int i = 0; i < S_SIZE; i++) sbuf[i] = 0;
length = ini_string.length;
for (int i = 0; i < length; i++) sbuf[i] = ini_string.sbuf[i];
}

// equivalence operator
bool MyString::operator==(const MyString &rhs)
{
bool retval = true;
if (length != rhs.length) retval = false;
else  for (int i = 0; i < length; i++) if (sbuf[i] != rhs.sbuf[i]) retval = false;
return retval;
}

// + operator

MyString& MyString::operator+(const MyString &rhs)
{
MyString* new_str = new MyString();
new_str->length = length+rhs.length;
if (new_str->length >= S_SIZE-1) new_str->length = S_SIZE-1;
for (int i = 0; i < length; i++) new_str->sbuf[i] = sbuf[i];
for (int i = length; i < new_str->length; i++) new_str->sbuf[i] = rhs.sbuf[i-length];
return *new_str;
}

// **** MyList classes mimicking limited TList functionality *******************
// **** currently all functions defined in .h file *****************************

