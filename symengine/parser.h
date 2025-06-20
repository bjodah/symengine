#ifndef SYMENGINE_PARSER_H
#define SYMENGINE_PARSER_H

#include <symengine/basic.h>

namespace SymEngine
{

struct ParserSettings {
    bool convert_xor {true};
    std::shared_ptr<std::map<const std::string, const RCP<const Basic>>> constants {};
};

std::shared_ptr<std::map<const std::string, const RCP<const Basic>>> get_default_parser_constants();

RCP<const Basic>
parse(const std::string &s, std::shared_ptr<const ParserSettings> settings={});
RCP<const Basic>
parse(const std::string &s, bool convert_xor);
RCP<const Basic>
parse(const std::string &s, bool convert_xor, const std::map<const std::string, const RCP<const Basic>> & constants);


RCP<const Basic> parse_old(const std::string &s, bool convert_xor = true);
RCP<const Basic>
parse_sbml(const std::string &s,
           const std::map<const std::string, const RCP<const Basic>> &constants
           = {});
} // namespace SymEngine

#endif
