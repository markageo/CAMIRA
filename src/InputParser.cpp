#include "InputParser.h"

#define NDEBUG

#include "boost/property_tree/ptree.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <map>
#include <cstdlib>
#include <utility>

// Characters used to define structure of input files
#define COMMENT_CHAR           '#'
#define BLOCK_OPEN_CHAR        '{'
#define BLOCK_CLOSE_CHAR       '}'
#define ASSIGNMENT_CHAR        '='

namespace pt = boost::property_tree;

/*-------------------------------------------------------------------------------------*\
                                     Helper functions
\*-------------------------------------------------------------------------------------*/

namespace 
{

// Error types
enum ErrorType {
    fileReadError,
    unmatchedBlockError,
    unbalancedBraceError,
    noKeyGivenError,
    noDataGivenError,
    noClosingQuotesError
};

// Parser states
enum ParserState {
    expectingKey,     // Expects key
    expectingData     // Expects data
};

//  Declarations
class ReaderStream;
void ParseStream(ReaderStream &, pt::ptree &);
void ProcessExpectingKeyState(ReaderStream &, std::stack<pt::ptree *> &, pt::ptree *&, ParserState &);
void ProcessExpectingDataState(ReaderStream &, std::stack<pt::ptree *> &, pt::ptree *&, ParserState &);
void DisplayErrorCode(ErrorType, ReaderStream &);


class ReaderStream
{
    public:

        ReaderStream(const std::string &inputFileName) : 
            m_inputFileName(inputFileName), m_inputFileStream(inputFileName), m_lineNumber(0) {};

        // Wrapper for getline function which reads new line and updates state of reader
        ReaderStream &ReadInputLine() 
        {
            std::getline(m_inputFileStream, m_inputLine);
            m_lineNumber++;
            m_linePos = m_inputLine.begin();
            SkipWhitespace();
            return *this;
        }

        // Return current character
        char CurrentCharacter() const 
        { return *m_linePos; }

        // Check if at end of current line
        bool AtEndOfLine() const 
        { return m_linePos == m_inputLine.end(); }

        // Go to next character in line (skipping whitespace)
        void GoToNextCharacter() { 
            SkipWhitespace();
            ++m_linePos; 
        }

        // Return line number
        unsigned LineNumber() const 
        { return m_lineNumber; }

        // Return input file name
        std::string InputFileName() const 
        { return m_inputFileName; }


        // Read key, which cannot have whitespace
        std::string ReadKey() {
            std::string key;
            SkipWhitespace();
            while (!SeperatorChar(*m_linePos) && m_linePos != m_inputLine.end() && !std::isspace(*m_linePos)) {
                key += *m_linePos;
                ++m_linePos;
            }
            SkipWhitespace();
            return key;
        }

        // Read data. Whitespace is removed unless it is inside quotes
        std::string ReadData() {
            std::string data;
            char quoteChar;
            SkipWhitespace();
            if        (*m_linePos == '"') {  // Double quotes
                quoteChar = '\"';
                ++m_linePos;
            } else if (*m_linePos == '\'') {  // Single quotes
                quoteChar = '\'';
                ++m_linePos;
            } else {                          // No quotes, ignore whitespace
                quoteChar = '\0';
            }

            while (1) {
                if (quoteChar == '\0') {
                    SkipWhitespace();
                    if (m_linePos == m_inputLine.end()) { 
                        break; 
                    }
                } else {
                    if (*m_linePos == quoteChar) {
                        ++m_linePos;
                        break;
                    }
                    if (m_linePos == m_inputLine.end()) {
                        throw noClosingQuotesError;
                        break;
                    }
                }
                data += *m_linePos;
                ++m_linePos;
            }
            SkipWhitespace();
            return data;
        }
        

        // Status checks
        explicit operator bool() const 
        { return !m_inputFileStream.fail(); }

        bool operator!() const 
        { return m_inputFileStream.fail(); }

    private:

         // State of the reader
        std::string m_inputFileName;
        std::ifstream m_inputFileStream;
        std::string m_inputLine;
        std::string::iterator m_linePos;
        unsigned m_lineNumber;

        // Advance line iterator until there is no whitespace
        void SkipWhitespace() {
            while(std::isspace(*m_linePos)) {
                ++m_linePos;
            }
        }

        // Check if character is a seperator charactor, i.e. block opening, comment, 
        // or equals sign
        bool SeperatorChar(const char ch) {
            if (ch == BLOCK_OPEN_CHAR       ||
                ch == BLOCK_CLOSE_CHAR      ||
                ch == ASSIGNMENT_CHAR ||
                ch == COMMENT_CHAR          ) { 
                    return true;
                } 
                else { 
                    return false;
                }
        }
};


// Outer parser loop function
void ParseStream(ReaderStream &readerStream, pt::ptree &tree) {
    
    // Initialise parser state to expect a key
    ParserState parserState = expectingKey;

    // Pointer to last created ptree
    pt::ptree *ptLast = NULL;

    // ptree stack to handle nesting
    std::stack<pt::ptree *> ptStack;
    ptStack.push(&tree);

    // Iterate each line
    while(readerStream.ReadInputLine()) {

        // Read through the line
        while (!readerStream.AtEndOfLine()) {
            
            // Ignore rest of line if there is comment
            if (readerStream.CurrentCharacter() == COMMENT_CHAR) { break; }

            switch (parserState)
            {
                case expectingKey:
                    ProcessExpectingKeyState(readerStream, ptStack, ptLast, parserState);
                    break;

                case expectingData:
                    ProcessExpectingDataState(readerStream, ptStack, ptLast, parserState);
                    break;
            }

        }

        // Key data pairs can only be given on a single line
        parserState = expectingKey;

    }
}


// Prcess data when parser is expecting key
void ProcessExpectingKeyState(ReaderStream &readerStream, std::stack<pt::ptree *> &ptStack, pt::ptree *&ptLast, ParserState &parserState) 
{       
    std::string key;
    switch (readerStream.CurrentCharacter()) 
    {
        case BLOCK_OPEN_CHAR:
            if (!ptLast) {
                throw unmatchedBlockError;
            }
            ptStack.push(ptLast);
            ptLast = NULL;
            readerStream.GoToNextCharacter();
            break;

        case BLOCK_CLOSE_CHAR:
            if (ptStack.size() <= 1) {
                throw unbalancedBraceError;
            }
            ptStack.pop();
            ptLast = NULL;
            readerStream.GoToNextCharacter();
            break;

        case ASSIGNMENT_CHAR:
            throw noKeyGivenError;
            break;

        default:
            key = readerStream.ReadKey();
            ptLast = &ptStack.top()->push_back( std::make_pair(key, pt::ptree()) )->second; // Add pointer to child of tree at top of the stack
            parserState = expectingData;
    }
}


// Prcess data when parser is expecting data
void ProcessExpectingDataState(ReaderStream &readerStream, std::stack<pt::ptree *> &ptStack, pt::ptree *&ptLast, ParserState &parserState) 
{       
    std::string data;
    switch (readerStream.CurrentCharacter()) 
    {
        case BLOCK_OPEN_CHAR:
            ptStack.push(ptLast);
            ptLast = NULL;
            readerStream.GoToNextCharacter();
            parserState = expectingKey;
            break;

        case BLOCK_CLOSE_CHAR:
            if (ptStack.size() <= 1) {
                throw unbalancedBraceError;
            }
            ptStack.pop();
            ptLast = NULL;
            readerStream.GoToNextCharacter();
            parserState = expectingKey;
            break;

        case ASSIGNMENT_CHAR:
            readerStream.GoToNextCharacter();
            data = readerStream.ReadData();
            if (data.empty()) {
                throw noDataGivenError;
            }
            ptLast->data() = data;
            parserState = expectingKey;
            break;

        // default:
     
    }
}


// Display error message corresponding to error code
void DisplayErrorCode(ErrorType error, ReaderStream &readerStream) 
{
    std::cout << "\n" << "FILE PARSING ERROR" << "\n";
    switch (error)
    {
        case fileReadError:
            std::cout <<  "Unable to open file '" + readerStream.InputFileName() << "." << "\n";
            break;

        case unmatchedBlockError:
            std::cout << "Line " << readerStream.LineNumber() << ": "
                      << "Block opened with '{'  without matching block name." << "\n";
            break;

        case unbalancedBraceError:
            std::cout << "Line " << readerStream.LineNumber() << ": "
                      << "Unbalanced block braces '{' and '}'." << "\n";
            break;

        case noKeyGivenError:
            std::cout << "Line " << readerStream.LineNumber() << ": "
                      << "Assignment operator '=' without corresponding key." << "\n";
            break;

        case noDataGivenError:
            std::cout << "Line " << readerStream.LineNumber() << ": "
                      << "Key assignment '=' specified without data." << "\n";
            break;

        case noClosingQuotesError:
            std::cout << "Line " << readerStream.LineNumber() << ": "
                      << "Unclosed string." << "\n";
            break;
        
    }
    std::cout << std::endl;
}


}


/*-------------------------------------------------------------------------------------*\
                                  User Functions
\*-------------------------------------------------------------------------------------*/

// Read input from file and store in SimData object
// return 0 for success
// return -1 for failure
int ReadInput(pt::ptree &pt, const std::string &inputFileName) 
{

    ReaderStream readerStream(inputFileName);

    try {

        if (!readerStream) {
            throw fileReadError;
        }
        
        // Start reading from the root space of the input file
        ParseStream(readerStream, pt);

        // Success
        return 0;

    } catch (ErrorType err) {

        // Display appropriate error message
        DisplayErrorCode(err, readerStream);
        return -1;

    }
    
}