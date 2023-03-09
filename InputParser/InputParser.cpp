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
#define COMMENT_CHAR_1         '/'
#define COMMENT_CHAR_2         '/' // Can make null character (\0) if comment symbol is only one character
#define BLOCK_OPEN_CHAR        '{'
#define BLOCK_CLOSE_CHAR       '}'
#define ASSIGNMENT_CHAR        '='
#define DIRECTIVE_CHAR         '#'
#define INCLUDE_DIRECTIVE      "include"

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
    noClosingQuotesError,
    invalidDirective,
    expectEndOfLineDirective,
    expectingString
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
void ProcessDirective(ReaderStream &, std::stack<pt::ptree *> &);
void DisplayErrorCode(ErrorType, ReaderStream &);


class ReaderStream
{
    public:

        ReaderStream(const std::string &inputFileName) : 
            m_filename(inputFileName), m_inputFileStream(inputFileName), m_lineNumber(0) {};

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

        // Check if at comment symbol
        bool AtCommentSymbol() const 
        { 
            if (*m_linePos == COMMENT_CHAR_1) {

                if (COMMENT_CHAR_2 == '\0')
                    return true; 

                if ( *(m_linePos+1) == COMMENT_CHAR_2 ) 
                    return true;
            } 
            return false;
        }

        // Go to next character in line (skipping whitespace)
        void GoToNextCharacter() { 
            SkipWhitespace();
            ++m_linePos; 
        }

        // Return line number
        unsigned LineNumber() const 
        { return m_lineNumber; }

        // Return input file name
        std::string Filename() const 
        { return m_filename; }


        // Read key, whitespace is removed
        std::string ReadKey() {
            std::string key;
            SkipWhitespace();
            while (1) {
                SkipWhitespace();
                if (SeperatorChar(*m_linePos))
                    break;

                if (AtCommentSymbol()) 
                    break;

                if (m_linePos == m_inputLine.end())
                    break;

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
                    if (m_linePos == m_inputLine.end()) 
                        break; 

                    if (SeperatorChar(*m_linePos))
                        break;
                    
                    if (AtCommentSymbol()) 
                        break;

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


        // Read type of directive, stops at whitespace
        std::string ReadDirectiveType() {
            std::string directive;
            if (*m_linePos == DIRECTIVE_CHAR)
                m_linePos++;
            while (1) {
                if (m_linePos == m_inputLine.end())
                    break;
                if (std::isspace(*m_linePos))
                    break;
                directive += *m_linePos;
                ++m_linePos;
            }
            SkipWhitespace();
            return directive;
        }
        

        // Status checks
        explicit operator bool() const 
        { return !m_inputFileStream.fail(); }

        bool operator!() const 
        { return m_inputFileStream.fail(); }

        // Advance line iterator until there is no whitespace
        void SkipWhitespace() {
            while(std::isspace(*m_linePos)) {
                ++m_linePos;
            }
        }

    private:

         // State of the reader
        std::string m_filename;
        std::ifstream m_inputFileStream;
        std::string m_inputLine;
        std::string::iterator m_linePos;
        unsigned m_lineNumber;


        // Check if character is a seperator charactor, i.e. block opening/closing 
        // or assignment
        bool SeperatorChar(const char ch) {
            if (ch == BLOCK_OPEN_CHAR    ||
                ch == BLOCK_CLOSE_CHAR   ||
                ch == ASSIGNMENT_CHAR   ) { 
                return true;
            } else { 
                return false;
            }
        }
};


// Outer parser loop function
void ParseStream(ReaderStream &readerStream, pt::ptree &tree) {
    
    // Initialise parser state to expect a key
    ParserState parserState = expectingKey;

    // Pointer to last created ptree
    pt::ptree *ptLast = nullptr;

    // ptree stack to handle nesting
    std::stack<pt::ptree *> ptStack;
    ptStack.push(&tree);

    // Iterate each line
    while(readerStream.ReadInputLine()) {

        // Skip whitespace so comment and directive chacaters can be seen
        readerStream.SkipWhitespace();

        // Ignore entire line if it starts with a comment
        if (readerStream.AtCommentSymbol()) 
            continue;

        // Directive
        if (readerStream.CurrentCharacter() == DIRECTIVE_CHAR) 
        { 
            ProcessDirective(readerStream, ptStack);  
            continue;
        }

        // Read through the line
        while (!readerStream.AtEndOfLine()) {
            
            // Ignore rest of line if there is comment
            if (readerStream.AtCommentSymbol()) 
                break;

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


// Process directives
void ProcessDirective(ReaderStream &readerStream, std::stack<pt::ptree *> &ptStack) 
{       

    // Directive type
    std::string directive;
    directive = readerStream.ReadDirectiveType();
    if (directive == INCLUDE_DIRECTIVE) {
        
        readerStream.SkipWhitespace();
        if (readerStream.CurrentCharacter() != '"' && readerStream.CurrentCharacter() != '\'')
            throw expectingString;

        ReaderStream readerStreamInclude( readerStream.ReadData() );
        if (!readerStream) 
            throw fileReadError;
        
        // Recursive call
        ParseStream(readerStreamInclude, *ptStack.top());

    } else {
        throw invalidDirective;
    }

    // There cannot be anything after the directive
    if ( !readerStream.AtEndOfLine() ) 
        throw expectEndOfLineDirective;


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
            ptLast = nullptr;
            readerStream.GoToNextCharacter();
            break;

        case BLOCK_CLOSE_CHAR:
            if (ptStack.size() <= 1) {
                throw unbalancedBraceError;
            }
            ptStack.pop();
            ptLast = nullptr;
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
            ptLast = nullptr;
            readerStream.GoToNextCharacter();
            parserState = expectingKey;
            break;

        case BLOCK_CLOSE_CHAR:
            if (ptStack.size() <= 1) {
                throw unbalancedBraceError;
            }
            ptStack.pop();
            ptLast = nullptr;
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
    std::cout << "\n" << "FILE PARSING ERROR: " 
              << "'" << readerStream.Filename() << "'" 
              << "\n";
    switch (error)
    {
        case fileReadError:
            std::cout <<  "Unable to open file '" + readerStream.Filename() << "." << "\n";
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

        case invalidDirective:
            std::cout << "Line " << readerStream.LineNumber() << ": "
                      << "Invalid directive." << "\n";
            break;

        case expectEndOfLineDirective:
            std::cout << "Line " << readerStream.LineNumber() << ": "
                      << "Expected end of line after directive." << "\n";
            break;

        case expectingString:
            std::cout << "Line " << readerStream.LineNumber() << ": "
                      << "Expected string." << "\n";
            break;
        
    }
    std::cout << std::endl;
}


}


/*-------------------------------------------------------------------------------------*\
                                  User Functions
\*-------------------------------------------------------------------------------------*/

// Read input from file and store in returned property tree
std::optional<pt::ptree> ParseFile(const std::string &inputFileName) 
{

    ReaderStream readerStream(inputFileName);
    pt::ptree pt;

    try {

        if (!readerStream)
            throw fileReadError;
        
        // Start reading from the root space of the input file
        ParseStream(readerStream, pt);

        // Success
        return pt;

    } catch (ErrorType err) {

        // Display appropriate error message
        DisplayErrorCode(err, readerStream);
        return {};
    }
    
}