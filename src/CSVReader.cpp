#include <string>
#include <vector>
#include <fstream>


enum class ReaderState {
    UnquotedField,
    QuotedField,
    QuotedQuote
};



void ProcessUnquotedField( std::string::const_iterator &rowIterator,
                           ReaderState &state,
                           std::vector< std::string > &rowFields )
{
    switch ( *rowIterator )
    {  
        case ',':
            rowFields.push_back("");
            break;

        case '"':
            state = ReaderState::QuotedField;
            break;

        default:
            rowFields.back() += *rowIterator;
            break;
    }
}


void ProcessQuotedField( std::string::const_iterator &rowIterator,
                         ReaderState &state,
                         std::vector< std::string > &rowFields )
{
    switch ( *rowIterator )
    {  
        case '"':
            state = ReaderState::QuotedQuote;
            break;

        default:
            rowFields.back() += *rowIterator;
            break;
    }
}


void ProcessQuotedQuote( std::string::const_iterator &rowIterator,
                         ReaderState &state,
                         std::vector< std::string > &rowFields )
{
    switch ( *rowIterator )
    {  
        case ',':
            rowFields.push_back("");
            state = ReaderState::UnquotedField;
            break;

        case '"':
            rowFields.back() += "\"";
            state = ReaderState::QuotedField;
            break;

        default:
            rowFields.back() += *rowIterator;
            break;
    }
}


std::vector< std::string > ReadRow( const std::string &row )
{
    std::vector< std::string > rowFields{""};
    ReaderState state = ReaderState::UnquotedField;

    for ( std::string::const_iterator rowIterator = row.begin(); rowIterator != row.end(); rowIterator++ ) {
        
        switch (state) 
        {
            case ReaderState::UnquotedField:
                ProcessUnquotedField( rowIterator, state, rowFields );
                break;

            case ReaderState::QuotedField:
                ProcessQuotedField( rowIterator, state, rowFields );
                break;


            case ReaderState::QuotedQuote:
                ProcessQuotedQuote( rowIterator, state, rowFields );
                break;

        }

    }

    return rowFields;
}



std::vector< std::vector< std::string > > ReadCSV( const std::string &filename )
{

    std::ifstream fStream( filename );
    std::vector< std::vector< std::string > > data;
    std::string rowString;

    while( !fStream.eof() ) {
        std::getline( fStream, rowString );
        if ( fStream.bad() || fStream.fail() ) {
            break;
        }

        std::vector< std::string > rowFields = ReadRow( rowString );
        data.push_back( rowFields );
    }

    return data;
}
