/* Generated by re2c 3.0 */
#line 1 "tokenizer.re"
#include "tokenizer.h"

#include "parser.tab.hh"

namespace SymEngine
{

void Tokenizer::set_string(const std::string &str)
{
    // The input string must be NULL terminated, otherwise the tokenizer will
    // not detect the end of string. After C++11, the std::string is guaranteed
    // to end with \0, but we check this here just in case.
    SYMENGINE_ASSERT(str[str.size()] == '\0');
    cur = (unsigned char *)(&str[0]);
}

int Tokenizer::lex(yy::parser::semantic_type* yylval)
{
    for (;;) {
        tok = cur;
        
#line 25 "tokenizer.cpp"
{
	unsigned char yych;
	static const unsigned char yybm[] = {
		  0,   0,   0,   0,   0,   0,   0,   0, 
		  0,  32,  32,  32,   0,  32,   0,   0, 
		  0,   0,   0,   0,   0,   0,   0,   0, 
		  0,   0,   0,   0,   0,   0,   0,   0, 
		 32,   0,   0,   0,   0,   0,   0,   0, 
		  0,   0,   0,   0,   0,   0,   0,   0, 
		192, 192, 192, 192, 192, 192, 192, 192, 
		192, 192,   0,   0,   0,   0,   0,   0, 
		  0, 128, 128, 128, 128, 128, 128, 128, 
		128, 128, 128, 128, 128, 128, 128, 128, 
		128, 128, 128, 128, 128, 128, 128, 128, 
		128, 128, 128,   0,   0,   0,   0, 128, 
		  0, 128, 128, 128, 128, 128, 128, 128, 
		128, 128, 128, 128, 128, 128, 128, 128, 
		128, 128, 128, 128, 128, 128, 128, 128, 
		128, 128, 128,   0,   0,   0,   0,   0, 
		128, 128, 128, 128, 128, 128, 128, 128, 
		128, 128, 128, 128, 128, 128, 128, 128, 
		128, 128, 128, 128, 128, 128, 128, 128, 
		128, 128, 128, 128, 128, 128, 128, 128, 
		128, 128, 128, 128, 128, 128, 128, 128, 
		128, 128, 128, 128, 128, 128, 128, 128, 
		128, 128, 128, 128, 128, 128, 128, 128, 
		128, 128, 128, 128, 128, 128, 128, 128, 
		128, 128, 128, 128, 128, 128, 128, 128, 
		128, 128, 128, 128, 128, 128, 128, 128, 
		128, 128, 128, 128, 128, 128, 128, 128, 
		128, 128, 128, 128, 128, 128, 128, 128, 
		128, 128, 128, 128, 128, 128, 128, 128, 
		128, 128, 128, 128, 128, 128, 128, 128, 
		128, 128, 128, 128, 128, 128, 128, 128, 
		128, 128, 128, 128, 128, 128, 128, 128, 
	};
	yych = *cur;
	if (yybm[0+yych] & 32) {
		goto yy4;
	}
	if (yych <= '?') {
		if (yych <= '*') {
			if (yych <= '%') {
				if (yych <= 0x00) goto yy1;
				if (yych <= 0x1F) goto yy2;
				if (yych <= '!') goto yy5;
				goto yy2;
			} else {
				if (yych == '\'') goto yy2;
				if (yych <= ')') goto yy6;
				goto yy8;
			}
		} else {
			if (yych <= '9') {
				if (yych == '.') goto yy9;
				if (yych <= '/') goto yy6;
				goto yy10;
			} else {
				if (yych <= '<') {
					if (yych <= ';') goto yy2;
					goto yy12;
				} else {
					if (yych <= '=') goto yy13;
					if (yych <= '>') goto yy14;
					goto yy2;
				}
			}
		}
	} else {
		if (yych <= '^') {
			if (yych <= 'O') {
				if (yych <= '@') goto yy15;
				if (yych == 'D') goto yy18;
				goto yy16;
			} else {
				if (yych <= 'P') goto yy19;
				if (yych <= 'Z') goto yy16;
				if (yych <= ']') goto yy2;
				goto yy6;
			}
		} else {
			if (yych <= '{') {
				if (yych == '`') goto yy2;
				if (yych <= 'z') goto yy16;
				goto yy2;
			} else {
				if (yych <= '}') {
					if (yych <= '|') goto yy6;
					goto yy2;
				} else {
					if (yych <= '~') goto yy6;
					if (yych <= 0x7F) goto yy2;
					goto yy16;
				}
			}
		}
	}
yy1:
	++cur;
#line 45 "tokenizer.re"
	{ return yy::parser::token::yytokentype::END_OF_FILE; }
#line 127 "tokenizer.cpp"
yy2:
	++cur;
yy3:
#line 44 "tokenizer.re"
	{ throw SymEngine::ParseError("Unknown token: '"+token()+"'"); }
#line 133 "tokenizer.cpp"
yy4:
	yych = *++cur;
	if (yybm[0+yych] & 32) {
		goto yy4;
	}
#line 46 "tokenizer.re"
	{ continue; }
#line 141 "tokenizer.cpp"
yy5:
	yych = *++cur;
	if (yych == '=') goto yy20;
	goto yy3;
yy6:
	++cur;
yy7:
#line 49 "tokenizer.re"
	{ return tok[0]; }
#line 151 "tokenizer.cpp"
yy8:
	yych = *++cur;
	if (yych == '*') goto yy15;
	goto yy7;
yy9:
	yych = *++cur;
	if (yych <= '/') goto yy3;
	if (yych <= '9') goto yy21;
	goto yy3;
yy10:
	yych = *++cur;
	if (yybm[0+yych] & 64) {
		goto yy10;
	}
	if (yych <= '^') {
		if (yych <= '@') {
			if (yych == '.') goto yy23;
		} else {
			if (yych == 'E') goto yy27;
			if (yych <= 'Z') goto yy24;
		}
	} else {
		if (yych <= 'd') {
			if (yych != '`') goto yy24;
		} else {
			if (yych <= 'e') goto yy27;
			if (yych <= 'z') goto yy24;
			if (yych >= 0x80) goto yy24;
		}
	}
yy11:
#line 58 "tokenizer.re"
	{ yylval->emplace<std::string>() = token(); return yy::parser::token::yytokentype::NUMERIC; }
#line 185 "tokenizer.cpp"
yy12:
	yych = *++cur;
	if (yych == '=') goto yy28;
	goto yy7;
yy13:
	yych = *++cur;
	if (yych == '=') goto yy29;
	goto yy3;
yy14:
	yych = *++cur;
	if (yych == '=') goto yy30;
	goto yy7;
yy15:
	++cur;
#line 50 "tokenizer.re"
	{ return yy::parser::token::yytokentype::POW; }
#line 202 "tokenizer.cpp"
yy16:
	yych = *++cur;
yy17:
	if (yybm[0+yych] & 128) {
		goto yy16;
	}
#line 57 "tokenizer.re"
	{ yylval->emplace<std::string>() = token(); return yy::parser::token::yytokentype::IDENTIFIER; }
#line 211 "tokenizer.cpp"
yy18:
	yych = *++cur;
	if (yych == 'e') goto yy31;
	goto yy17;
yy19:
	yych = *++cur;
	if (yych == 'i') goto yy32;
	goto yy17;
yy20:
	++cur;
#line 53 "tokenizer.re"
	{ return yy::parser::token::yytokentype::NE; }
#line 224 "tokenizer.cpp"
yy21:
	yych = *++cur;
yy22:
	if (yych <= '^') {
		if (yych <= '@') {
			if (yych <= '/') goto yy11;
			if (yych <= '9') goto yy21;
			goto yy11;
		} else {
			if (yych == 'E') goto yy27;
			if (yych <= 'Z') goto yy24;
			goto yy11;
		}
	} else {
		if (yych <= 'd') {
			if (yych == '`') goto yy11;
			goto yy24;
		} else {
			if (yych <= 'e') goto yy27;
			if (yych <= 'z') goto yy24;
			if (yych <= 0x7F) goto yy11;
			goto yy24;
		}
	}
yy23:
	yych = *++cur;
	if (yych == 'E') goto yy24;
	if (yych != 'e') goto yy22;
yy24:
	yych = *++cur;
yy25:
	if (yych <= '^') {
		if (yych <= '9') {
			if (yych >= '0') goto yy24;
		} else {
			if (yych <= '@') goto yy26;
			if (yych <= 'Z') goto yy24;
		}
	} else {
		if (yych <= '`') {
			if (yych <= '_') goto yy24;
		} else {
			if (yych <= 'z') goto yy24;
			if (yych >= 0x80) goto yy24;
		}
	}
yy26:
#line 59 "tokenizer.re"
	{ yylval->emplace<std::string>() = token(); return yy::parser::token::yytokentype::IMPLICIT_MUL; }
#line 274 "tokenizer.cpp"
yy27:
	yych = *(mar = ++cur);
	if (yych <= ',') {
		if (yych == '+') goto yy33;
		goto yy25;
	} else {
		if (yych <= '-') goto yy33;
		if (yych <= '/') goto yy25;
		if (yych <= '9') goto yy35;
		goto yy25;
	}
yy28:
	++cur;
#line 51 "tokenizer.re"
	{ return yy::parser::token::yytokentype::LE; }
#line 290 "tokenizer.cpp"
yy29:
	++cur;
#line 54 "tokenizer.re"
	{ return yy::parser::token::yytokentype::EQ; }
#line 295 "tokenizer.cpp"
yy30:
	++cur;
#line 52 "tokenizer.re"
	{ return yy::parser::token::yytokentype::GE; }
#line 300 "tokenizer.cpp"
yy31:
	yych = *++cur;
	if (yych == 'r') goto yy36;
	goto yy17;
yy32:
	yych = *++cur;
	if (yych == 'e') goto yy37;
	goto yy17;
yy33:
	yych = *++cur;
	if (yych <= '/') goto yy34;
	if (yych <= '9') goto yy35;
yy34:
	cur = mar;
	goto yy26;
yy35:
	yych = *++cur;
	if (yych <= '^') {
		if (yych <= '9') {
			if (yych <= '/') goto yy11;
			goto yy35;
		} else {
			if (yych <= '@') goto yy11;
			if (yych <= 'Z') goto yy24;
			goto yy11;
		}
	} else {
		if (yych <= '`') {
			if (yych <= '_') goto yy24;
			goto yy11;
		} else {
			if (yych <= 'z') goto yy24;
			if (yych <= 0x7F) goto yy11;
			goto yy24;
		}
	}
yy36:
	yych = *++cur;
	if (yych == 'i') goto yy38;
	goto yy17;
yy37:
	yych = *++cur;
	if (yych == 'c') goto yy39;
	goto yy17;
yy38:
	yych = *++cur;
	if (yych == 'v') goto yy40;
	goto yy17;
yy39:
	yych = *++cur;
	if (yych == 'e') goto yy41;
	goto yy17;
yy40:
	yych = *++cur;
	if (yych == 'a') goto yy42;
	goto yy17;
yy41:
	yych = *++cur;
	if (yych == 'w') goto yy43;
	goto yy17;
yy42:
	yych = *++cur;
	if (yych == 't') goto yy44;
	goto yy17;
yy43:
	yych = *++cur;
	if (yych == 'i') goto yy45;
	goto yy17;
yy44:
	yych = *++cur;
	if (yych == 'i') goto yy46;
	goto yy17;
yy45:
	yych = *++cur;
	if (yych == 's') goto yy47;
	goto yy17;
yy46:
	yych = *++cur;
	if (yych == 'v') goto yy48;
	goto yy17;
yy47:
	yych = *++cur;
	if (yych == 'e') goto yy49;
	goto yy17;
yy48:
	yych = *++cur;
	if (yych == 'e') goto yy50;
	goto yy17;
yy49:
	yych = *++cur;
	if (yybm[0+yych] & 128) {
		goto yy16;
	}
#line 56 "tokenizer.re"
	{ yylval->emplace<std::string>() = token(); return yy::parser::token::yytokentype::PIECEWISE; }
#line 396 "tokenizer.cpp"
yy50:
	yych = *++cur;
	if (yybm[0+yych] & 128) {
		goto yy16;
	}
#line 55 "tokenizer.re"
	{ yylval->emplace<std::string>() = token(); return yy::parser::token::yytokentype::DERIVATIVE; }
#line 404 "tokenizer.cpp"
}
#line 60 "tokenizer.re"

    }
}

} // namespace SymEngine
