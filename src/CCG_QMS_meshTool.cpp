#include"CCG_QMS_meshTool.h"
#include"../3rdParty/MeshLib/core/Parser/parser.h"

/*---------------- V start------------------------------*/
void CCG_QMSLib::CCG_QMS_Vertex::_from_string()
{
	CParser parser(m_string);
	for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
	{
		CToken* token = *iter;
		if (token->m_key == "uv") //CPoint2
		{
			token->m_value >> m_uv;
		}
		if (token->m_key == "feature")
		{
			this->feature() = true;
		}
		if (token->m_key == "originalId")
		{
			//std::cout <<"token->m_value: " << token->m_value << std::endl;
			std::string t = token->m_value;
			t.erase(0, t.find_first_not_of("()"));
			t.erase(t.find_last_not_of("()") + 1);
			m_originalId = strutil::parseString<int>(t);
			//std::cout << " m_originalId: " << m_originalId << std::endl;
		}
	}
}

void CCG_QMSLib::CCG_QMS_Vertex::_to_string()
{
	CParser parser(m_string);
	parser._removeToken("uv");

	parser._toString(m_string);
	std::stringstream iss;

	iss << "uv=(" << m_uv[0] << " " << m_uv[1] << ")";

	if (m_string.length() > 0)
	{
		m_string += " ";
	}
	m_string += iss.str();
}
/*---------------- V end------------------------------*/

/*---------------- E start------------------------------*/
void CCG_QMSLib::CCG_QMS_Edge::_from_string()
{
	CParser parser(m_string);
	//std::cout << m_string << std::endl;
	for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
	{
		CToken* token = *iter;
		if (token->m_key == "sharp")
		{
			this->sharp() = true;
			this->feature() = true;
			this->constrinedBoundary() = true;
			//token->m_value >> m_uv;
		}
		if (token->m_key == "smoothFeature")
		{
			this->constrinedBoundary() = true;
			//token->m_value >> m_uv;
		}
		//if (token->m_key == "rgb") // CPoint
		//{
		//	token->m_value >> m_rgb;
		//}
	}
}

void CCG_QMSLib::CCG_QMS_Edge::_to_string()
{
	CParser parser(m_string);
	parser._removeToken("sharp");
	parser._toString(m_string);
	if (this->sharp())
	{

		std::stringstream iss;
		iss << "sharp";

		if (m_string.length() > 0)
		{
			m_string += " ";
		}
		m_string += iss.str();
	}
}
/*---------------- E end------------------------------*/

/*---------------- F start------------------------------*/
void CCG_QMSLib::CCG_QMS_Face::_from_string()
{
	CParser parser(m_string);
	//std::cout << m_string << std::endl;
	for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
	{
		CToken* token = *iter;
		if (token->m_key == "patchIndex")
		{
			//std::cout <<"token->m_value: " << token->m_value << std::endl;
			std::string t = token->m_value;
			t.erase(0, t.find_first_not_of("()"));
			t.erase(t.find_last_not_of("()") + 1);
			m_patchIndex = strutil::parseString<int>(t);
			//std::cout <<" m_patchIndex: " << m_patchIndex << std::endl;
		}
		if (token->m_key == "subPatchIndex")
		{
			//std::cout << "token->m_value: " << token->m_value << std::endl;
			std::string t = token->m_value;
			t.erase(0, t.find_first_not_of("()"));
			t.erase(t.find_last_not_of("()") + 1);
			m_subPatchIndex = strutil::parseString<int>(t);
			//std::cout << " m_subPatchIndex: " << m_subPatchIndex << std::endl;
		}
	}
}

void CCG_QMSLib::CCG_QMS_Face::_to_string()
{
	CParser parser(m_string);
	parser._removeToken("patchIndex");
	parser._removeToken("subPatchIndex");

	parser._toString(m_string);
	std::stringstream iss;

	iss << " patchIndex=(" << m_patchIndex << ")" << " " << "subPatchIndex=(" << m_subPatchIndex << ") ";

	if (m_string.length() > 0)
	{
		m_string += " ";
	}
	m_string += iss.str();
}
/*---------------- F end------------------------------*/

/*---------------- H start------------------------------*/
void CCG_QMSLib::CCG_QMS_HalfEdge::_from_string()
{
	CParser parser(m_string);
	for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
	{
		CToken* token = *iter;
		if (token->m_key == "uv") //CPoint2
		{
			token->m_value >> m_uv;
			//std::cout << "Corner " << this->target()->id() << " " << this->face()->id() << " { uv=(" <<m_uv[0]<<" "<<m_uv[1]<<") }" << std::endl;
		}
		if (token->m_key == "id_FeatureTrajectory")
		{
			//std::cout << "token->m_value: " << token->m_value << std::endl;
			std::string t = token->m_value;
			t.erase(0, t.find_first_not_of("()"));
			t.erase(t.find_last_not_of("()") + 1);
			m_id_FeatureTrajectory = strutil::parseString<int>(t);
			//std::cout << "Halfedge "<<this->source()->id()<<" "<<this->target()->id() << " { m_id_FeatureTrajectory= " << m_id_FeatureTrajectory<<" } " << std::endl;
		}
	}
}

void CCG_QMSLib::CCG_QMS_HalfEdge::_to_string()
{
	CParser parser(m_string);
	parser._removeToken("uv");

	parser._toString(m_string);
	std::stringstream iss;

	iss << "uv=(" << m_uv[0] << " " << m_uv[1] << ")";

	if (m_string.length() > 0)
	{
		m_string += " ";
	}
	m_string += iss.str();
}
/*---------------- H end------------------------------*/

/*---------------- M start------------------------------*/
/*---------------- M end------------------------------*/
