/*
*	- (c++) Yiming Zhu
*	- 2019/09/10
*	- Algorithm : Parse the trait strings of {vertex, edge, halfedge, face} base "c" language (fast version)
*/

#ifndef _MyLib_Parser_H_
#define _MyLib_Parser_H_

#include "My_IOFuncDef.h"
#include<cstring>

namespace MyLib 
{
	//Return type
	enum ReturnType { StringOver = 1, TokenOut = 2, Success = 3, Failure = 4 };

	/*!
	 *	\brief My_Token  class, key=(value), e.g. uv=(x y)
	*/

	class My_Token
	{
	public:
		/*! key of the token */
		char mKey[MAX_TOKEN_STRING_SIZE];
		/*! value of the token */
		char mValue[MAX_TOKEN_STRING_SIZE];
	};

	/*!
	 *	\brief My_Parser class
	*/

	class My_Parser
	{
	public:
		/*!
		 *	\brief My_Parser constructor
		 *  \param str input string
		 */
		My_Parser(const char * str)
		{
			unsigned int len = strlen(str);
			if (len >= MAX_LINE_SIZE) {
				printf("Error: Current string's length is over MAX LINE SIZE!\n");
				return;
			}
			/*Copy string*/
			unsigned int i;
			for (i = 0; i < len; i++)
			{
				mLine[i] = str[i];
			}
			mLine[i] = 0;
			currrentPtr = mLine;
			int returnValue = nextToken();
			switch (returnValue)
			{
			case MyLib::StringOver:
				//printf("Error: Can't find Value! In line :  %s", mLine);
				break;
			case MyLib::Failure:
				//printf("Error: Can't find Value! In line :  %s", mLine);
				break;
			case MyLib::TokenOut:
				//printf("Error: Token's length is over! In line :  %s", mLine);
				break;
			default:
				break;
			}
		};

		/*!
		 *	My_Parser Destructor
		 */
		~My_Parser() {};

		/*!
		 *	List of tokens extracted from the string
		 */
		My_Token* tokens() { return &mToken; };

		/*Jump to next token*/
		int nextToken() {
			/*Skip the blank*/
			SkipBlank();
			/*If string over, return false*/
			if (End()) return StringOver;
			/*Get the key from current line*/
			char * pKey = mKey;
			char   currentChar = NextChar();
			while (currentChar != ' ' && currentChar != '=' && !End())
			{
				*pKey++ = currentChar;
				currentChar = NextChar();
			}
			if (currentChar != '=' && currentChar != ' ')
			{
				*pKey++ = currentChar;
			}
			*pKey = 0;
			while (currentChar == ' ' && !End())
			{
				currentChar = NextChar();
			}
			if (currentChar != '=')
			{
				int keySize = strlen(mKey);
				if (keySize >= MAX_TOKEN_STRING_SIZE) {
					return TokenOut;
				}
				SAFE_STRCPY(mToken.mKey, mKey);
				mToken.mKey[keySize] = '\0';
				return Failure;
			}
			/*If string over, return false*/
			if (End()) return StringOver;

			currentChar = NextChar();

			while (currentChar != '(' && !End()) currentChar = NextChar();

			char * pvalue = mValue;

			while (currentChar != ')' && !End())
			{
				*pvalue++ = currentChar;
				currentChar = NextChar();
			}
			*pvalue++ = currentChar;
			*pvalue = 0;

			int keySize = strlen(mKey);
			if (keySize >= MAX_TOKEN_STRING_SIZE) {
				return TokenOut;
			}
			SAFE_STRCPY(mToken.mKey, mKey);
			mToken.mKey[keySize] = '\0';

			int valueSize = strlen(mValue);
			if (valueSize >= MAX_TOKEN_STRING_SIZE) {
				return TokenOut;
			}
			SAFE_STRCPY(mToken.mValue, mValue);
			mToken.mValue[valueSize] = '\0';

			return Success;
		}

	private:

		/*!
		 *	get the next char in the current string
		 */
		char NextChar()
		{
			char currentChar = *currrentPtr;
			currrentPtr++;
			return currentChar;
		};
		/*!
		 *	skip blank spaces
		 */
		void SkipBlank()
		{
			while (*currrentPtr == ' ')
			{
				currrentPtr++;
			}
		};
		/*!
		 *	verify if we
		 */
		bool End()
		{
			return ((*currrentPtr) == 0);
		};
		/*!
		 *	token
		 */
		My_Token mToken;
		/*!
		 * The buffer to contain the string
		 */
		char mLine[MAX_LINE_SIZE];
		/*!
		 *	current key
		 */
		char mKey[MAX_TOKEN_STRING_SIZE];
		/*!
		 *	current value
		 */
		char mValue[MAX_TOKEN_STRING_SIZE];
		/*!
		 *	current pointer to the char buffer
		 */
		char * currrentPtr;
	};
}
#endif //_MyLib_Parser_H_ defined