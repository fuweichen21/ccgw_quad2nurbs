/*****************************************************************//**
 * \file   NUBE_printInformation.h
 * \brief  
 * 
 * \author FuweiChen
 * \date   March 2025
 *********************************************************************/

#ifndef _NUBE_PRINTINFORMATION_H_
#define _NUBE_PRINTINFORMATION_H_
#include<string>

namespace NUBELib {
	class NUBE_printInformation
	{
	public:
		NUBE_printInformation() { ; };
		~NUBE_printInformation() { ; };
		void PrintLogInformation(
			const std::string logInformation, const bool printToScreen_Log, FILE* filePointer_Log);

	};

}

#endif // !_NUBE_PRINTINFORMATION_H_

