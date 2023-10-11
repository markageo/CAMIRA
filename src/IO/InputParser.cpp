#include "InputParser.h"

#include "boost/property_tree/ptree.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <map>
#include <cstdlib>
#include <utility>
#include <filesystem>


namespace INP
{

// Characters used to define structure of input files
#define COMMENT_CHAR_1         '/'
#define COMMENT_CHAR_2         '/' // Can make null character (\0) if comment symbol is only one character
#define BLOCK_OPEN_CHAR        '{'
#define BLOCK_CLOSE_CHAR       '}'
#define ASSIGNMENT_CHAR        '='
#define DIRECTIVE_CHAR         '#'

// Directives
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
std::string ErrorString(const ErrorType &, const ReaderStream &);


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

        // Return relative path to input file (LINUX FILESYSTEMS ONLY)
        std::string RelativePath() const
        { 
            std::string::const_iterator stringIterator = m_filename.end();
            for ( /* NULL */ ; *stringIterator != '/' && stringIterator != m_filename.begin(); stringIterator--) {};
            return std::string( m_filename.begin(), stringIterator );
        }


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

            enum DataState {
                Quoted, Unquoted
            };

            std::string data;
            DataState dataState = Unquoted;
            bool keepReadingLine = true;

            while ( m_linePos != m_inputLine.end() && keepReadingLine ) {
                
                switch ( dataState ) {
                    case Quoted:

                        if ( *m_linePos == '"' ) {
                            dataState = Unquoted;
                            break;
                        }

                        if ( m_linePos == ( m_inputLine.end()-1 )  &&  *m_linePos != '"' ) {
                            throw std::runtime_error( ErrorString(noClosingQuotesError, *this) );
                        }

                        data += *m_linePos;
                        break;

                    case Unquoted:
                        
                        SkipWhitespace();

                        if ( AtCommentSymbol() || SeperatorChar(*m_linePos) ) {
                            keepReadingLine = false;
                            break;
                        }
                            

                        if ( *m_linePos == '"' ) {
                            dataState = Quoted;
                            break;
                        }

                        data += *m_linePos;
                        break;
                }

                m_linePos++;
                if ( dataState == Unquoted )
                    SkipWhitespace();

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
            throw std::runtime_error( ErrorString(expectingString, readerStream) );

        std::string includeFilename = readerStream.ReadData();
        std::string processPath = std::filesystem::current_path();
        std::string includeFullPath = processPath + "/" + readerStream.RelativePath() + "/" + includeFilename;
        
        ReaderStream readerStreamInclude( includeFullPath );
        if (!readerStream) 
            throw std::runtime_error( ErrorString(fileReadError, readerStream) );
        
        // Recursive call
        ParseStream(readerStreamInclude, *ptStack.top());

    } else {
        throw std::runtime_error( ErrorString(invalidDirective, readerStream) );
    }

    // There cannot be anything after the directive
    if ( !readerStream.AtEndOfLine() ) 
        throw std::runtime_error( ErrorString(expectEndOfLineDirective, readerStream) );


}

// Prcess data when parser is expecting key
void ProcessExpectingKeyState(ReaderStream &readerStream, std::stack<pt::ptree *> &ptStack, pt::ptree *&ptLast, ParserState &parserState) 
{       
    std::string key;
    switch (readerStream.CurrentCharacter()) 
    {
        case BLOCK_OPEN_CHAR:
            if (!ptLast) {
                throw std::runtime_error( ErrorString(unmatchedBlockError, readerStream) );
            }
            ptStack.push(ptLast);
            ptLast = nullptr;
            readerStream.GoToNextCharacter();
            break;

        case BLOCK_CLOSE_CHAR:
            if (ptStack.size() <= 1) {
                throw std::runtime_error( ErrorString(unbalancedBraceError, readerStream) );
            }
            ptStack.pop();
            ptLast = nullptr;
            readerStream.GoToNextCharacter();
            break;

        case ASSIGNMENT_CHAR:
            throw std::runtime_error( ErrorString(noKeyGivenError, readerStream) );
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
                throw std::runtime_error( ErrorString(unbalancedBraceError, readerStream) );
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
                throw std::runtime_error( ErrorString(noDataGivenError, readerStream) );
            }
            ptLast->data() = data;
            parserState = expectingKey;
            break;

        // default:
     
    }
}


// Return error message
std::string ErrorString(const ErrorType &error, const ReaderStream &readerStream) 
{
    std::string errorMessage;

    switch (error)
    {
        case fileReadError:
            errorMessage = "Unable to open file '" + readerStream.Filename() + "'.";
            break;

        case unmatchedBlockError:
            errorMessage = "Line " + std::to_string(readerStream.LineNumber())
                         + ": Block opened with '{'  without matching block name.";
            break;

        case unbalancedBraceError:
            errorMessage = "Line " + std::to_string(readerStream.LineNumber()) 
                         + ": Unbalanced block braces '{' and '}'.";
            break;

        case noKeyGivenError:
            errorMessage = "Line " + std::to_string(readerStream.LineNumber())
                         + ": Assignment operator '=' without corresponding key.";
            break;

        case noDataGivenError:
            errorMessage = "Line " + std::to_string(readerStream.LineNumber())
                         + ": Key assignment '=' specified without data.";
            break;

        case noClosingQuotesError:
            errorMessage = "Line " + std::to_string(readerStream.LineNumber())
                         + ": Unclosed string.";
            break;

        case invalidDirective:
            errorMessage = "Line " + std::to_string(readerStream.LineNumber())
                         + "Invalid directive.";
            break;

        case expectEndOfLineDirective:
            errorMessage = "Line " + std::to_string(readerStream.LineNumber())
                         + ": Expected end of line after directive.";
            break;

        case expectingString:
            errorMessage = "Line " + std::to_string(readerStream.LineNumber())
                         + "Expected string.";
            break;
        
        default:
            errorMessage = "undefined error!";
        
    }
    
    return errorMessage;
}


}


/*-------------------------------------------------------------------------------------*\
                                  User Functions
\*-------------------------------------------------------------------------------------*/

// Read input from file and store in returned property tree
pt::ptree ParseFile(const std::string &inputFileName) 
{

    ReaderStream readerStream(inputFileName);
    pt::ptree pt;

    if (!readerStream)
        throw std::runtime_error( ErrorString(fileReadError, readerStream) );
    
    // Start reading from the root space of the input file
    ParseStream(readerStream, pt);

    // Success
    return pt;
    
}

}   // end namespace INP