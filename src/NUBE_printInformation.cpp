#include "NUBE_printInformation.h"

void NUBELib::NUBE_printInformation::PrintLogInformation(const std::string logInformation, const bool printToScreen_Log, FILE* filePointer_Log)
{
	if (printToScreen_Log == true)
	{
		printf(logInformation.c_str());
	}
	if (filePointer_Log != NULL)
	{
		fprintf(filePointer_Log, logInformation.c_str());
	}
}
