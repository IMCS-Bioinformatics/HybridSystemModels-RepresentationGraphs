//---------------------------------------------------------------------------

#ifndef my_bcbH
#define my_bcbH
//---------------------------------------------------------------------------

// **** MyString class mimicking limited AnsiString functionality **************

class MyString;

// maximum size of strings
#define S_SIZE 4096

class MyString
{
public:
    // constructors
    MyString(); // create empty string
    MyString(char* ini_string); // copy string provided as argument
    MyString(MyString* ini_string); // copy MyString provided as argument
    MyString(int ini_n); // integer provides as argument
    // copy constructor
    MyString(const MyString &ini_string);
    // equivalence operator
    bool operator==(const MyString &rhs);
    // + operator
    MyString& operator+(const MyString &rhs);
    // values
    char* c_str() {return sbuf;} // return char pointer to string buffer
    char sbuf[S_SIZE];
    int length;
};

// **** MyList classes mimicking limited TList functionality *******************
// **** no control for exceeding array size ************************************

// maximum size of TLists
# define TL_LARGE 1048576 // value 2^20
# define TL_SMALL 128 // value 2^7

class MyListS
{
public:
    // constructor
    MyListS() {Count=0;}
    // function
    void Clear() {Count=0;}
    void Add(void* item) {Items[Count]=item; Count++;}
    // values
    void* Items[TL_SMALL];
    int Count;
};

class MyListL
{
public:
    // constructor
    MyListL() {Count=0;}
    // function
    void Clear() {Count=0;}
    void Add(void* item) {Items[Count]=item; Count++;}
    // values
    void* Items[TL_LARGE];
    int Count;
};

#endif
