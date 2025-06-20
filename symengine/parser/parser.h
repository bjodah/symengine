#ifndef SYMENGINE_PARSER_PARSER_H
#define SYMENGINE_PARSER_PARSER_H

#include <fstream>
#include <algorithm>
#include <memory>

#include <symengine/add.h>
#include <symengine/pow.h>
#include <symengine/logic.h>
#include <symengine/parser.h>

namespace SymEngine
{

/*
   To Parse (default) constructor is expensive as it creates all the maps and
   tables. If just one expression needs to be parsed, then calling
   SymEngine::parse() does the job. But if multiple expressions are to be
   parsed, then first initialize SymEngine::Parser and after that call
   SymEngine::Parser::parse() repeatedly.

   Example:

   Parser p;
   auto r = p.parse("x**2");

*/

class Tokenizer;

class Parser
{
protected:
    std::string inp;
    std::shared_ptr<const ParserSettings> settings;

public:
    std::unique_ptr<Tokenizer> m_tokenizer;
    RCP<const Basic> res;

    RCP<const Basic> parse(const std::string &input);

    RCP<const Basic> functionify(const std::string &name, vec_basic &params);
    RCP<const Basic> parse_numeric(const std::string &expr);
    RCP<const Basic> parse_identifier(const std::string &expr);
    std::tuple<RCP<const Basic>, RCP<const Basic>>
    parse_implicit_mul(const std::string &expr);
    explicit Parser(std::shared_ptr<const ParserSettings> settings);
    Parser(const std::map<const std::string, const RCP<const Basic>> & constants);
    Parser();
    ~Parser();
};

} // namespace SymEngine

#endif
